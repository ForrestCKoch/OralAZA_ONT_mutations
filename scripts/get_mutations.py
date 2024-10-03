import argparse
from bisect import bisect_left
from collections import Counter
import copy
import itertools
import multiprocessing as mp
import os
import pathlib
import re

import edlib as ed
import networkx as nx
import pandas as pd
import pysam
from tqdm import tqdm

def build_parser():
    """
    Construct an argument parser
    """
    parser = argparse.ArgumentParser(
            prog='get_consensus.py',
            description='Generates a consensus read from unique barcodes',
            epilog='Written by Forrest C. Koch (forrest.c.koch@gmail.com) 2024')

    parser.add_argument('-s', '--sam', required=True,
                        type=str, help='Input sam/bam.')

    parser.add_argument('-m', '--mutations', default='mutations-of-interest-reduced.csv',
                        type=str, help='File detailing the mutations to be searched for.')

    parser.add_argument('-o', '--output', required=True, type=str, 
                        help='Output csv.')

    parser.add_argument('-c', '--cores', required=False, type=int, default=8,
                        help='The number of cores to use')

    #parser.add_argument('-q', '--quiet', default=False, type=bool,
    #                    help='Flag to supress message output')

    return parser

def get_index(x, a):
    """
    Locate the leftmost value exactly equal to x
    https://docs.python.org/3/library/bisect.html
    """
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    return None

def get_mutations(read, mutation_df):
    """
    Work out which mutations are relevant to
    """

    # A dictionary to keep track of whether each read can be assigned to a mutation
    #bp_dict = {mutation:None for mutation in mutation_df['coord_hg38']}
    bp_dict = {}

    for index, row in mutation_df.iterrows():
        # Check that this read has actually been assigned to the gene
        if '|'+row['gene']+'|' not in read.qname:
            # No match, increment NA
            # No need to save this in the dictionary -- should help to keep memory usage down
            continue

        # Check for coverage 
        mutation_chr = row['coord_hg38'].split(':')[0]
        if read.to_dict()['ref_name'] != mutation_chr:
            continue
        
        full_reference_positions = read.get_reference_positions(full_length=True)
        reference_positions = read.get_reference_positions(full_length=False)
        mutation_position = int(row['coord_hg38'].split(':')[1])
        target_idx = get_index(mutation_position, reference_positions)
        if target_idx is None:
            #bp_dict[row['coord_hg38']] = 'NC'
            continue

        # If we reach this point, match has been found, increment that BP
        bp = read.get_forward_sequence()[full_reference_positions.index(mutation_position)]
        bp_dict[row['coord_hg38']] = bp

    return bp_dict



def get_mutation_dict(fetch, header, mutation_df, quiet=False, cores=8):
    mutation_dict = {}
    idx_counter = 0

    # Put together an iterator over all of the reads
    # We need to convert reads to dictionaries which can be pickled and
    # used by the multiprocessing library
    if(not quiet):
        print("Processing reads ...")
    read_list = [read.to_dict() for read in tqdm(fetch)]
    arg_iterator = zip(read_list, 
                       itertools.cycle([header.to_dict()]), 
                       itertools.cycle([mutation_df]))

    # Send the iterator to imap to search the target mutation locations in each read
    if(not quiet):
        print("Searching target positions ...")
    pool = mp.Pool(cores)
    results = tqdm(pool.imap(_get_mutation_dict_iter, arg_iterator))

    # Sort the results into a dictionary grouping common cells and umi's
    for cb, umi, bp_dict in results:
        if cb not in mutation_dict:
            mutation_dict[cb] = {}
        if umi not in mutation_dict[cb]:
            mutation_dict[cb][umi] = \
                    {mutation:{'A':0, 'C':0, 'G':0, 'T':0}\
                    for mutation in mutation_df['coord_hg38']}
        # Add the current index to the list and increment counter
        idx_counter += 1 
        for mutation in  bp_dict:
            # Increment the findings for the unique barcode
            mutation_dict[cb][umi][mutation][bp_dict[mutation]] += 1

    pool.close()

    return mutation_dict

def _get_mutation_dict_iter(arg):
    read_dict = arg[0]
    header_dict = arg[1]
    mutation_df = arg[2]
    header = pysam.AlignmentHeader.from_dict(header_dict)
    read = pysam.AlignedSegment.from_dict(read_dict, header = header)

    tags = dict(read.get_tags())

    cb = tags['CB']
    umi = tags['UB']

    bp_dict = get_mutations(read, mutation_df)

    return (cb, umi, bp_dict)

def cluster_mutation_dict_umis(mutation_dict, quiet=False):
    """
    process the output of get_mutation_dict to cluster common UMI's together
    """ 

    if not quiet:
        print("Clustering UMIs ...")
    for cb in tqdm(mutation_dict):
        # For each cell, grab the umi's and cluster them
        umi_list = list(mutation_dict[cb])
        umi_corrections = cluster_umis(umi_list)
        # Then, aggregate the totals within each cluster
        for umi in umi_list:
            correction = umi_corrections[umi]
            if umi == correction:
                continue # No adjustment necessary .. continue
            # Otherwise, add this entry to the "correction" entry, and delete the original 
            for mutation in mutation_dict[cb][umi]:
                # The mutation may not yet be recoreded in the target ..
                # If this crashes, may need to adjust ...
               for bp in mutation_dict[cb][correction][mutation]:
                    # Add the umi to the target correction
                    mutation_dict[cb][correction][mutation][bp] += \
                    mutation_dict[cb][umi][mutation][bp]
            # Remove the entry as it has been transfered over to the target correction
            del mutation_dict[cb][umi] 
 

    return None

def cluster_umis(umis):
    """
    Return a dictionary providing umi corrections based on louvain partitions
    """
    
    # Build our graph of umi's directly adjacent to one another
    edges = {umis[i]:\
            [umis[j] for j in range(i, len(umis)) \
                     if ed.align(umis[i], umis[j], k=1)['editDistance']==1\
            ] for i in range(len(umis))\
            }
    graph = nx.Graph(edges)
    #cliques = graph.find_cliques()

    # Take a very coarse grain partitioning of the umi's
    partitions = next(nx.community.louvain_partitions(graph))

    # Create a dictionary assigining a common umi to each unique umi
    umi_corrections = {}
    for p in partitions:
        # work out which umi has the greatest degree, this is likely the "true" umi
        max_degree = 0
        max_umi = ''
        for umi in p:
            degree = graph.degree[umi]
            max_degree, max_umi = (max_degree, max_umi) if max_degree > degree else (degree, umi)
        
        # assign that value to all umi's 
        for umi in p:
            umi_corrections[umi] = max_umi
        
    return umi_corrections

def assign_umi_mutations(mutation_dict, quiet=False):
    
    if not quiet:
        print("Assigning UMI mutations ...")

    for cb in tqdm(mutation_dict):
        for umi in mutation_dict[cb]:
            # Find the maximal base pair for each mutation and max that assignment
            for mutation in mutation_dict[cb][umi]:
                max_bp = None
                max_count = 0
                for bp, count in mutation_dict[cb][umi][mutation].items():
                    max_bp, max_count = (bp, count) if count > max_count else (max_bp, max_count)
                # Make the assignment
                mutation_dict[cb][umi][mutation] = max_bp #(max_bp, max_count) 

    return None

def get_cb_mutations(mutation_dict, quiet=False):

    cb_mutations = {}
    

    if not quiet:
        print("Assigning cell mutations ...")

    for cb in tqdm(mutation_dict):
        cb_mutations[cb] = {'Total':0}

        # Go through each of the UMIs and tally up how many base pairs we see at each location
        for umi in mutation_dict[cb]:
            cb_mutations[cb]['Total'] += 1
            max_mutation = None
            max_count = 0
            for mutation, bp in mutation_dict[cb][umi].items():
                if mutation not in cb_mutations[cb]:
                    cb_mutations[cb][mutation] = {}
                if bp not in cb_mutations[cb][mutation]:
                    cb_mutations[cb][mutation][bp] = 0
                cb_mutations[cb][mutation][bp] += 1

        # Go back through and work out the most likely status for each possible mutation
        for mutation in cb_mutations[cb]:
            # One of these entries isn't a dictionary ...
            if mutation == 'Total':
                continue
            max_bp = None
            max_count = 0
            coverage = 0
            for bp, count in cb_mutations[cb][mutation].items():
                # Skip if bp is none ...
                if bp is None:
                    continue
                max_bp, max_count = (bp, count) if count > max_count else (max_bp, max_count)
                coverage += count
            cb_mutations[cb][mutation] = (max_bp, max_count, coverage) 

    return cb_mutations

def write_cb_mutations_to_csv(cb_mutations, outfile, quiet=False):
   
    df_dict = {} 

    if not quiet:
        print("Writing to file ...")

    for cb in tqdm(cb_mutations):
        df_dict[cb] = {}
        df_dict[cb]['Total UMIs'] = cb_mutations[cb]['Total']
        
        for mutation in cb_mutations[cb]:
            # One of these entries isn't a dictionary ...
            if mutation == 'Total':
                continue
            bp, count, coverage = cb_mutations[cb][mutation]
            df_dict[cb][mutation+'_bp'] = bp
            df_dict[cb][mutation+'_count'] = count
            df_dict[cb][mutation+'_coverage'] = coverage

    pd.DataFrame.from_dict(df_dict, orient='index').to_csv(outfile)

def main(args):

    sam = pysam.AlignmentFile(args.sam)
    fetch = sam.fetch(multiple_iterators=True)

    mutation_df = pd.read_csv(args.mutations)#'mutations-of-interest-reduced.csv')

    # Get a data structure which has grouped reads into common umi's and cell barcodes
    # which indicates, for each location of interest, the count of each base pair across reads
    mutation_dict = get_mutation_dict(fetch, sam.header, mutation_df, args.cores)

    # Process the dictionary to cluster similar umi's together to account for read errors
    cluster_mutation_dict_umis(mutation_dict)
   
    # Assign mutations by replacing the dictionary of base pair counts with just the maximal base 
    assign_umi_mutations(mutation_dict)

    cb_mutations = get_cb_mutations(mutation_dict)

    write_cb_mutations_to_csv(cb_mutations, args.output)

    return None
 
if __name__ == '__main__':
    parser = build_parser()

    args = parser.parse_args()

    if not os.path.isfile(args.sam):
        raise FileNotFoundError('Input sam {} not found'.format(args.sam))

    main(args)

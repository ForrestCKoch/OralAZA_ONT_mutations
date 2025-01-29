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
import scipy
from tqdm import tqdm

def build_parser():
    """
    Construct an argument parser
    """
    parser = argparse.ArgumentParser(
            prog='get_mutations.py',
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

def get_mutations(read, mt_df):
    """
    Work out which mutations are relevant to
    """

    # A dictionary to keep track of whether each read can be assigned to a mutation
    mt_status = {}

    for index, row in mt_df.iterrows():
        # Check that this read has actually been assigned to the gene
        if '|'+row['gene']+'|' not in read.qname:
            # No match, increment NA
            # No need to save this in the dictionary -- should help to keep memory usage down
            continue

        # Check for coverage 
        mutation_chr = row['mut2__hg38'].split(':')[0]
        if read.to_dict()['ref_name'] != mutation_chr:
            continue
        
        full_ref_pos = read.get_reference_positions(full_length=True)
        ref_pos = read.get_reference_positions(full_length=False)
        mt_pos = int(row['mut2__hg38'].split(':')[1].split('_')[0]) - 1
        mt_desc = row['mut2__hg38'].split('_')[1]
        target_idx = get_index(mt_pos, ref_pos)
        if target_idx is None:
            # No overlap, skip to next iteration
            continue

        # If we reach this point, match has been found, increment that BP
        mt_status[row['mut2__hg38']] = assign_mutation_status(\
                sequence = read.query_sequence, \
                positions = full_ref_pos, \
                mt_pos = mt_pos, \
                mt_desc = mt_desc)

    return mt_status

def assign_mutation_status(sequence, positions, mt_pos, mt_desc, tol=0.05, err=0.1):
    # NOTE, tol and err are used to allow for slight mismatches in longer insertion 
    # sequences. The number of allowed mistmatches is calculated as
    # scipy.stats.binom(1-tol, len(mt_str), err) where mt_str is the insertion 
    # sequence (including the first, unchanged basepair)

    wt_str, mt_str = mt_desc.split('>')
    mt_status = None
    
    start_pos = positions.index(mt_pos)

    if (len(wt_str)==1) and (len(mt_str)==1):
        # SNP
        # See if the observed read matches either of the "possibilities"
        bp = sequence[start_pos]
        return 'MT' if mt_str==bp else 'WT' if wt_str==bp else 'UN'

    elif (len(wt_str)>1) and (mt_str==wt_str[0]):
        # Deletion
        # Search for a gap with the appropriate size
        if len(positions)==(start_pos+1):
            # Somehow this is the end of the read... tough luck
            return 'UN'

        if (positions[start_pos] is None) or (positions[start_pos+1] is None): 
            # Not sure what really to do in the case ....
            return 'UN'

        gap = positions[start_pos+1] - positions[start_pos]
        # If the gap is the same length as wt_str, then we have a deletion that matches
        # If the gap is 1, there is no deletion
        # Otherwise, we don't know
        return 'MT' if gap==len(wt_str) else 'WT' if gap==1 else 'UN'

    elif (len(mt_str)>1) and (wt_str==mt_str[0]):
        # Insertion
        # Check that the reference alignment is 'None' just after the target position 
        # for the correct number of base pairs

        if len(positions)==(start_pos+1):
            # Somehow this is the end of the read... tough luck
            return 'UN'

        if positions[start_pos+1] is not None:
            # No insertion found, assign WT and break
            return 'WT'

        if len(positions)<(start_pos+len(mut_str)):
            # The read doesn't cover the length of the mutation
            return 'UN'

        # Check that the correct number of inserted base pairs is observed
        for i in range(2, len(mut_str)):
            if positions[start_pos+i] is not None:      
                return 'UN'

        # Finally, check that the sequence matches, within a given tolerence
        max_ed = scipy.stats.binom.ppf(1-tol, len(mut_str), err)
        obs_ed = ed.align(sequence[start_pos:(start_pos+len(mut_str))], mut_str)['editDistance']
        # Return MT if it is within tolerence, otherwise Unknown ...
        return 'UN' if ob_ed > max_ed else 'MT'
    else:
        # Unknown
        raise(ValueError, "Unsupported mutation string: {}".format(mt_desc))
    
    return 'UN'


def get_mutation_dict(fetch, header, mt_df, quiet=False, cores=8):
    mt_dict = {}
    idx_counter = 0

    # Put together an iterator over all of the reads
    # We need to convert reads to dictionaries which can be pickled and
    # used by the multiprocessing library
    if(not quiet):
        print("Processing reads ...")
    read_list = [read.to_dict() for read in tqdm(fetch)]
    arg_iterator = zip(read_list, 
                       itertools.cycle([header.to_dict()]), 
                       itertools.cycle([mt_df]))

    # Send the iterator to imap to search the target mutation locations in each read
    if(not quiet):
        print("Searching target positions ...")
    pool = mp.Pool(cores)
    results = tqdm(pool.imap(_get_mutation_dict_iter, arg_iterator))

    # Sort the results into a dictionary grouping common cells and umi's
    for cb, umi, mt_status in results:
        if cb not in mt_dict:
            mt_dict[cb] = {}
        if umi not in mt_dict[cb]:
            mt_dict[cb][umi] = \
                    {mutation:{'WT':0, 'MT':0, 'UN':0}\
                    for mutation in mt_df['mut2__hg38']}
        # Add the current index to the list and increment counter
        idx_counter += 1 
        for mutation in mt_status:
            # Increment the findings for the unique barcode
            mt_dict[cb][umi][mutation][mt_status[mutation]] += 1

    pool.close()

    return mt_dict

def _get_mutation_dict_iter(arg):
    read_dict = arg[0]
    header_dict = arg[1]
    mt_df = arg[2]
    header = pysam.AlignmentHeader.from_dict(header_dict)
    read = pysam.AlignedSegment.from_dict(read_dict, header = header)

    tags = dict(read.get_tags())

    cb = tags['CB']
    umi = tags['UB']

    mts_dict = get_mutations(read, mt_df)

    return (cb, umi, mts_dict)

def cluster_mutation_dict_umis(mt_dict, quiet=False):
    """
    process the output of get_mutation_dict to cluster common UMI's together
    """ 

    if not quiet:
        print("Clustering UMIs ...")
    for cb in tqdm(mt_dict):
        # For each cell, grab the umi's and cluster them
        umi_list = list(mt_dict[cb])
        umi_corrections = cluster_umis(umi_list)
        # Then, aggregate the totals within each cluster
        for umi in umi_list:
            correction = umi_corrections[umi]
            if umi == correction:
                continue # No adjustment necessary .. continue
            # Otherwise, add this entry to the "correction" entry, and delete the original 
            for mutation in mt_dict[cb][umi]:
                # The mutation may not yet be recoreded in the target ..
                # If this crashes, may need to adjust ...
               for mts in mt_dict[cb][correction][mutation]:
                    # Add the umi to the target correction
                    mt_dict[cb][correction][mutation][mts] += \
                    mt_dict[cb][umi][mutation][mts]
            # Remove the entry as it has been transfered over to the target correction
            del mt_dict[cb][umi] 
 

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

def assign_umi_mutations(mt_dict, quiet=False):
    
    if not quiet:
        print("Assigning UMI mutations ...")

    for cb in tqdm(mt_dict):
        for umi in mt_dict[cb]:
            # Find the maximal mutation type for each mutation and max that assignment
            for mutation in mt_dict[cb][umi]:
                max_mts = None
                max_count = 0
                for mts, count in mt_dict[cb][umi][mutation].items():
                    max_mts, max_count = (mts, count) if count > max_count else (max_mts, max_count)
                # Make the assignment
                mt_dict[cb][umi][mutation] = max_mts #(max_mts, max_count) 

    return None

def get_cb_mutations(mt_dict, quiet=False):

    cb_mutations = {}
    

    if not quiet:
        print("Assigning cell mutations ...")

    for cb in tqdm(mt_dict):
        cb_mutations[cb] = {'Total':0}

        # Go through each of the UMIs and tally up how many base pairs we see at each location
        for umi in mt_dict[cb]:
            cb_mutations[cb]['Total'] += 1
            max_mutation = None
            max_count = 0
            for mutation, mts in mt_dict[cb][umi].items():
                if mutation not in cb_mutations[cb]:
                    cb_mutations[cb][mutation] = {}
                if mts not in cb_mutations[cb][mutation]:
                    cb_mutations[cb][mutation][mts] = 0
                cb_mutations[cb][mutation][mts] += 1

        # Go back through and work out the most likely status for each possible mutation
        for mutation in cb_mutations[cb]:
            # One of these entries isn't a dictionary ...
            if mutation == 'Total':
                continue
            #max_mts = None
            #max_count = 0
            #coverage = 0
            #for mts, count in cb_mutations[cb][mutation].items():
            #    # Skip if mts is none ...
            #    if mts is None:
            #        continue
            #    max_mts, max_count = (mts, count) if count > max_count else (max_mts, max_count)
            #    coverage += count

            n_wt = cb_mutations[cb][mutation]['WT'] if 'WT' in cb_mutations[cb][mutation] else 0
            n_mt = cb_mutations[cb][mutation]['MT'] if 'MT' in cb_mutations[cb][mutation] else 0
            n_un = cb_mutations[cb][mutation]['UN'] if 'UN' in cb_mutations[cb][mutation] else 0

            cb_mutations[cb][mutation] = (n_wt, n_mt, n_un)#(max_mts, max_count, coverage) 

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

            wt, mt, un = cb_mutations[cb][mutation]
            coverage = wt + mt + un

            df_dict[cb][mutation+'_wt'] = wt
            df_dict[cb][mutation+'_mt'] = mt
            df_dict[cb][mutation+'_un'] = un
            df_dict[cb][mutation+'_coverage'] = coverage

    pd.DataFrame.from_dict(df_dict, orient='index').to_csv(outfile)

def main(args):

    sam = pysam.AlignmentFile(args.sam)
    fetch = sam.fetch(multiple_iterators=True)

    mt_df = pd.read_csv(args.mutations)#'mutations-of-interest-reduced.csv')

    # Get a data structure which has grouped reads into common umi's and cell barcodes
    # which indicates, for each location of interest, the count of each base pair across reads
    mt_dict = get_mutation_dict(fetch, sam.header, mt_df, args.cores)

    # Process the dictionary to cluster similar umi's together to account for read errors
    cluster_mutation_dict_umis(mt_dict)
   
    # Assign mutations by replacing the dictionary of base pair counts with just the maximal base 
    assign_umi_mutations(mt_dict)

    cb_mutations = get_cb_mutations(mt_dict)

    write_cb_mutations_to_csv(cb_mutations, args.output)

    return None
 
if __name__ == '__main__':
    parser = build_parser()

    args = parser.parse_args()

    if not os.path.isfile(args.sam):
        raise FileNotFoundError('Input sam {} not found'.format(args.sam))

    main(args)

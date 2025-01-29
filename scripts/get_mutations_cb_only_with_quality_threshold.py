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
import numpy as np
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

    parser.add_argument('-q', '--quality-threshold', required=False, type=int, default=25,
                        help='Exclude reads with PHRED score lower than `quality_threshold` at the mutation site')

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
        mt_pos = int(row['mut2__hg38'].split(':')[1].split('_')[0]) - 1 # The -1 is important if the mutation description is 1 indexed beacuse pysam is 0 indexed
        mt_desc = row['mut2__hg38'].split('_')[1]
        target_idx = get_index(mt_pos, ref_pos)
        if target_idx is None:
            # No overlap, skip to next iteration
            continue

        # If we reach this point, match has been found, increment that BP
        mt_status[row['mut2__hg38']] = assign_mutation_status(\
                sequence = read.query_sequence, \
                positions = full_ref_pos, \
                qualities = read.query_qualities, \
                mt_pos = mt_pos, \
                mt_desc = mt_desc)

    return mt_status

def assign_mutation_status(sequence, positions, qualities, mt_pos, mt_desc, tol=0.05, err=0.1, default_qual=30):
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
        qual = qualities[start_pos]
        return ('MT', qual) if mt_str==bp else ('WT', qual) if wt_str==bp else ('UN', qual)

    elif (len(wt_str)>1) and (mt_str==wt_str[0]):
        # Deletion
        # Search for a gap with the appropriate size
        if len(positions)==(start_pos+1):
            # Somehow this is the end of the read... tough luck
            return ('UN', default_qual)

        if (positions[start_pos] is None) or (positions[start_pos+1] is None): 
            # Not sure what really to do in the case ....
            return ('UN', default_qual)

        gap = positions[start_pos+1] - positions[start_pos]
        # If the gap is the same length as wt_str, then we have a deletion that matches
        # If the gap is 1, there is no deletion
        # Otherwise, we don't know
        qual = (qualities[start_pos]+qualities[start_pos+1])/2
        return ('MT', qual) if gap==len(wt_str) else ('WT', qual) if gap==1 else ('UN', qual)

    elif (len(mt_str)>1) and (wt_str==mt_str[0]):
        # Insertion
        # Check that the reference alignment is 'None' just after the target position 
        # for the correct number of base pairs

        qual = qualities[start_pos]
        qual_n = 1

        if len(positions)==(start_pos+1):
            # Somehow this is the end of the read... tough luck
            return ('UN', default_qual)

        if positions[start_pos+1] is not None:
            # No insertion found, assign WT and break
            qual += qualities[start_pos+1]
            qual_n += 1
            return ('WT', qual/qual_n)

        if len(positions)<(start_pos+len(mut_str)):
            # The read doesn't cover the length of the mutation
            return ('UN', default_qual)

        # Check that the correct number of inserted base pairs is observed
        for i in range(2, len(mut_str)):
            qual += qualities[start_pos+i]
            qual_n += 1
            if positions[start_pos+i] is not None:      
                return ('UN', qualities[start_pos+i])

        # Finally, check that the sequence matches, within a given tolerence
        max_ed = scipy.stats.binom.ppf(1-tol, len(mut_str), err)
        obs_ed = ed.align(sequence[start_pos:(start_pos+len(mut_str))], mut_str)['editDistance']
        # Return MT if it is within tolerence, otherwise Unknown ...
        return ('UN', qual/qual_n) if ob_ed > max_ed else ('MT', qual/qual_n)
    else:
        # Unknown
        raise(ValueError, "Unsupported mutation string: {}".format(mt_desc))
    
    return 'UN'


def get_mutation_dict(fetch, header, mt_df, quiet=False, cores=8, quality_threshold=25):
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
    for cb, mt_status in results:
        if cb not in mt_dict:
            mt_dict[cb] = \
                    {mutation:{'WT':np.array((0.0,0.0)), 'MT':np.array((0.0,0.0)), 'UN':np.array((0.0,0.0))}\
                    for mutation in mt_df['mut2__hg38']}
        # Add the current index to the list and increment counter
        idx_counter += 1 
        for mutation in mt_status:
            # Increment the findings for the unique barcode
            status, qual = mt_status[mutation]
            # Skip if the quality is unacceptable
            if qual < quality_threshold:
                continue
            # Otherwise add to the dictionary 
            mt_dict[cb][mutation][status] += np.array((1.0, np.float64(qual).item()))

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

    mts_dict = get_mutations(read, mt_df)

    return (cb, mts_dict)



def get_cb_mutations(mt_dict, quiet=False):

    cb_mutations = {}
    

    if not quiet:
        print("Re-formatting cell mutations ...")

    for cb in tqdm(mt_dict):
        cb_mutations[cb] = {}

        # Go through each of the UMIs and tally up how many base pairs we see at each location
        for mutation, mts in mt_dict[cb].items():
            if mutation not in cb_mutations[cb]:
                cb_mutations[cb][mutation] = {}

            n_wt = mt_dict[cb][mutation]['WT'][0] if 'WT' in mt_dict[cb][mutation] else 0
            n_mt = mt_dict[cb][mutation]['MT'][0] if 'MT' in mt_dict[cb][mutation] else 0
            n_un = mt_dict[cb][mutation]['UN'][0] if 'UN' in mt_dict[cb][mutation] else 0

            q_wt = mt_dict[cb][mutation]['WT'][1] if 'WT' in mt_dict[cb][mutation] else 0
            q_mt = mt_dict[cb][mutation]['MT'][1] if 'MT' in mt_dict[cb][mutation] else 0
            q_un = mt_dict[cb][mutation]['UN'][1] if 'UN' in mt_dict[cb][mutation] else 0

            cb_mutations[cb][mutation] = (n_wt, q_wt, n_mt, q_mt, n_un, q_un)

    return cb_mutations

def write_cb_mutations_to_csv(cb_mutations, outfile, quiet=False):
   
    df_dict = {} 

    if not quiet:
        print("Writing to file ...")

    for cb in tqdm(cb_mutations):
        df_dict[cb] = {}
        
        for mutation in cb_mutations[cb]:

            wt, qwt, mt, qmt, un, qun = cb_mutations[cb][mutation]
            coverage = wt + mt + un

            df_dict[cb][mutation+'_wt'] = wt
            df_dict[cb][mutation+'_mt'] = mt
            df_dict[cb][mutation+'_un'] = un

            df_dict[cb][mutation+'_wt_qual'] = qwt/wt if wt > 0 else None
            df_dict[cb][mutation+'_mt_qual'] = qmt/mt if mt > 0 else None
            df_dict[cb][mutation+'_un_qual'] = qun/un if un > 0 else None

            df_dict[cb][mutation+'_coverage'] = coverage

    pd.DataFrame.from_dict(df_dict, orient='index').to_csv(outfile)

def main(args):

    sam = pysam.AlignmentFile(args.sam)
    fetch = sam.fetch(multiple_iterators=True)

    mt_df = pd.read_csv(args.mutations)#'mutations-of-interest-reduced.csv')

    # Get a data structure which has grouped reads into common umi's and cell barcodes
    # which indicates, for each location of interest, the count of each base pair across reads
    mt_dict = get_mutation_dict(fetch, sam.header, mt_df, args.cores, args.quality_threshold)

    cb_mutations = get_cb_mutations(mt_dict)

    write_cb_mutations_to_csv(cb_mutations, args.output)

    return None
 
if __name__ == '__main__':
    parser = build_parser()

    args = parser.parse_args()

    if not os.path.isfile(args.sam):
        raise FileNotFoundError('Input sam {} not found'.format(args.sam))

    main(args)

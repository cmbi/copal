"""
multiple timewarp module: contains functions that coordinate progressive alignment

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md
"""

# import statements
from . import timewarp
import numpy as np

# functions
def create_sample_entries(totaldata,samplelengths):
    """
    Creates sample entries compatible with timewarp function for each sample

    Args:
        totaldata: 3D list structure with complexome profile sample intensity values
        samplelengths: list of ints, length of each sample
    Returns:
        entries (dict): dictionary with sample entries compatible with timewarp function
                        structure: {sample(int): [sample slice indices (int)](list)}
    """
    entries = {}
    for sample in range(len(totaldata)):
        entries[str(sample+1)] = {str(sample+1):range(samplelengths[sample])}
    return entries

def pairwise_timewarp(totaldata,samplelengths,localdict):
    """
    performs pairwise alignments (first step of progressive alignment)

    Args:
        totaldata: 3D list structure with complexome profile sample intensity values
        samplelengths: list of ints, length of each sample
        localdict: dict with 2D local cost value grid for each sample pair
    Returns:
        pairwisedict (dict): dictionary with alignments results for each pair
                            pair entry: {'samplenum:samplenum'(str): alignment(dict)}
        pairwisecosts: (dict): keys -> alignment costs, values -> 'sample1:sample2' pairs
    """
    entries = create_sample_entries(totaldata,samplelengths)
    pairwisedict = {}
    pairwisecosts = {}
    samples = sorted(entries.keys())[:-1]
    targets = sorted(entries.keys())[1:]
    for sample in samples: # loop that goes over all pairs, without including self-matches
        for target in targets:
            pairwarp = timewarp.timewarp(localdict,entries[sample], entries[target])
            pairwisedict[sample + ':' + target] =  pairwarp[0]
            pairwisecosts[pairwarp[1]] = sample + ':' + target
        del targets[0]
    return (pairwisedict,pairwisecosts)

def multiple_timewarp(pairwisecosts,samplelengths,localdict,pairwisedict):
    """
    determines progressive aligment order and performs progressive alignment

    Args:
        pairwisecosts (dict): pairwise alignment costs(keys), sample pairs (values)
        samplelengths: list of ints, length of each sample
        localdict: dict with 2D local cost value grid for each sample pair
        pairwisedict (dict): dictionary with alignments results for each pair
                            pair entry: {'samplenum:samplenum'(str): alignment(dict)}
    Returns:
        final alignment(dict): samples as keys (int), aligned slice indices values (list)
            multiple_alignment_order (listof tuples): multiple alignment order
                                      format:[((align1 samples),(align2 samples), cost)]
    """
    #sort pairwise alignments from pairs with lowest cost to pairs with highest cost
    keys = sorted(pairwisecosts.keys())
    msa_order = []                      # list with alignments have been performed with[(al1,al2,cost), ...]
    alignments = {}
    for key in keys:                    # loop through pairwise alignments. pair value looks like: '1:2'

    #determine for each pair in sorted list if the 2 samples have already been aligned, either seperately or in the same alignment
        pair = pairwisecosts[key]
        sample_in_alignments = []       # shows for both samples in pair whether and if so, in which alignment it is using None/slicenr of alignments list
        for sample in [0,1]:            # loops through 2 samples in each pair
            sample_loc = None
            for al_key in alignments.keys():   #for each alignment, checks whether sample is included in it. if so, stores location in sample_loc
                if int(pair.split(':')[sample]) in al_key:
                    sample_loc = al_key
            sample_in_alignments.append(sample_loc)     # for each sample, appends sample_value to sample_in_alignments

    #perform time warping alignment for all possible scenarios of each pair
        if sample_in_alignments[0] == None and sample_in_alignments[1] == None:  # both samples haven't been aligned --> take pairwise alignment
            alignments[(int(pair.split(':')[0]),int(pair.split(':')[1]))] = pairwisedict[pair]       # get alignment from pairwisedict, store in alignments
            msa_order.append(((int(pair.split(':')[0])),(int(pair.split(':')[1])), key))                 # add alignment info to msa_order

        elif sample_in_alignments[0] == sample_in_alignments[1]:           # samples already in same alignment --> skip this pair
            pass  # skip this alignment

        elif sample_in_alignments[0] != None and sample_in_alignments[1] != None:    # samples are in different alignments --> align these
            align1key = sample_in_alignments[0]                   # get keys from the 2 alignments
            align2key = sample_in_alignments[1]
            new_key = tuple(list(align1key)+list(align2key))                         # create new key by combining old ones
            alignments[new_key] = timewarp.timewarp(localdict,alignments[align1key], alignments[align2key])[0]   # perform alignment, store in alignments
            msa_order.append((align1key,align2key, key))
            del alignments[align1key]                               # delete old alignments
            del alignments[align2key]

        elif sample_in_alignments[0] == None:                             #sample 1 not yet aligned, 2 has been --> align 1 with 2's alignment
            new_align = {pair.split(':')[0]:range(samplelengths[int(pair.split(':')[0])-1])}       # create new entry for unaligned sample
            align2key = sample_in_alignments[1]        # get key for 2's alignment
            new_key = tuple(list(align2key) + [int(pair.split(':')[0])])             # generate new key by adding sample 1
            alignments[new_key] = timewarp.timewarp(localdict,new_align,alignments[align2key])[0]
            msa_order.append(((int(pair[0])), align2key, key))
            del alignments[align2key]

        elif sample_in_alignments[1] == None:                             # sample 2 not yet aligned, 1 has been --> align 2 with 1's alignment
            new_align = {pair.split(':')[1]:range(samplelengths[int(pair.split(':')[1])-1])}
            align1key = sample_in_alignments[0]
            new_key = tuple(list(align1key) + [int(pair.split(':')[1])])
            alignments[new_key]= timewarp.timewarp(localdict,new_align,alignments[align1key])[0]
            msa_order.append((align1key, (int(pair.split(':')[1])), key))
            del alignments[align1key]

    return (next(iter(alignments.values())), msa_order)

def datawarp(data, alignment):
    """
    uses final alignemnt to warp complexome profiling data along the 'slice axis'

    warps protein migration patterns by repeating intensity values at alignment gaps
    Args:
        data: 3D list structure with complexome profile sample intensity values
                   of unaligned samples
        final alignment(dict): samples as keys (int), aligned slice indices values (list)
    Returns:
        warpeddata: 3D list structure with aligned intensity values for each sample
    """
    warpeddata = []
    samplenum = 0
    for sample in data:
        newsample = []
        sample_alignment = alignment[str(samplenum+1)]       #get alignment 'x values' for this sample from alignment dictionary
        for x in sample_alignment:                          # loop through all alignment x values for this sample
            newsample.append(data[samplenum][x])            # append slice corresponding to alignment x value to new sample
        warpeddata.append(newsample)                          # add new sample to warpeddata
        samplenum += 1
    return warpeddata

def interpolate_sample(sample,alignment):
    """
    warp sample, using interpolation to fill gaps
    """
    dimensions = (sample.shape[0],len(alignment))
    warped_data = np.full(dimensions, np.nan)
    alignment_iterator = enumerate(alignment)
    for ix,value in alignment_iterator:
        if value != None:
            warped_data[:,ix] = sample[:,value]
        else:
            gaps = 1
            # look into the future for any extra Nones
            for value in alignment[ix+1:]:
                if value == None:
                    gaps += 1
                else:
                    break

            # error if sample starts with a gap
            if ix == 0:
                raise ValueError("sample starts with gap. "
                "Option not implemented as this normally can't happen")

            # get adjacent values
            left_adj = sample[:,alignment[ix-1]]
            try:
                right_adj = sample[:,alignment[ix+gaps]]
        
            # catch and handle case where gaps last till end of sample
            except IndexError:
                fill_values = left_adj.repeat(gaps).reshape(dimensions[0],gaps)
                warped_data[:,ix:] = fill_values
                break
            except Exception:
                raise ValueError('something else went wrong with right_adj..')

            # determine step sizes
            diff = right_adj - left_adj
            stepsizes = diff/(gaps + 1)

            # get values to fill gaps
            gap_array = np.arange(1,gaps+1)
            to_add = np.outer(gap_array, stepsizes)
            fill_values = (to_add + left_adj).transpose()
            warped_data[:,ix:ix+gaps] = fill_values

            # skip iterations over gap indexes that have been handled            
            for i in range(gaps-1):
                next(alignment_iterator)
    
    return warped_data

def find_gaps(alignment):
    """
    """
    previous = None
    with_gaps = []
    for value in alignment:
        if value == previous:
            with_gaps.append(None)
        else:
            with_gaps.append(value)
        previous = value
    return with_gaps

def interpolation_warp(data, alignment):
    """
    """
    warpeddata = []
    for ix,sample in enumerate(data):
        sample_array = np.array(sample).transpose()
        sample_alignment = alignment[str(ix+1)]
        sample_gaps = find_gaps(sample_alignment)
        warped_sample_array =interpolate_sample(sample_array,sample_gaps)
        warped_sample = warped_sample_array.transpose().tolist()
        warpeddata.append(warped_sample)
    return warpeddata

if __name__ == "__main__":
    # get test sample
    sample1 = [[0,3,8,1],[1,3,7,1],[2,1,6,1],[6,1,5,1],[7,1,4,1]]
    sample2 = [[3,4,6,0],[8,8,6,0],[4,3,8,0],[3,7,6,0],[7,2,0,0]]
    test_data = [sample1,sample2]
    test_alignment = {
        '1':[0,1,2,2,3,3,4,4],
        '2':[0,1,1,2,2,3,3,4],
                }

    warped_data = interpolation_warp(test_data, test_alignment)
    
    print(warped_data[0])
    print(warped_data[1])

    print(np.array(warped_data[0]).transpose())
    print(np.array(warped_data[1]).transpose())
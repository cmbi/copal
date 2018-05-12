""" 
shiftscore module, uses data with hausdorff distances to compute hausdorff effect size

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md
"""

# import statements
import pandas as pd
import numpy as np


#functions

# determine comparisons to be made
def comp_between(groups):
    """
    determines the sample combinations to be evaluated to calculate scores between groups 
    
    Args:
        groups (list): contains two lists with sample numbers belonging to each group
    Returns:
        comparisons (list): list of sample pairs to be evaluated: ['samplenum:samplenum'] 
    """ 
    comparisons = []
    for sample in groups[0]:
        for target in groups[1]:
            if int(sample) < int(target):
                comparisons.append(str(sample) + ':' + str(target))
            else:
                comparisons.append(str(target) + ':' + str(sample))
    return comparisons

def comp_within(groups):
    """
    determines the sample combinations to be evaluated to calculate scores within groups 
    
    Args:
        groups (list): contains two lists with sample numbers belonging to each group
    Returns:
        comparisons (list): list of sample pairs to be evaluated: ['samplenum:samplenum'] 
    """ 
    comparisons = []
    for group in groups:
        samples = group[:-1]
        targets = group[1:]
        for sample in samples:
            for target in targets:
                if int(sample) < int(target):
                    comparisons.append(str(sample) + ':' + str(target))
                else:
                    comparisons.append(str(target) + ':' + str(sample))
            del targets[0]
    return comparisons
    
def get_sample_name_headers(sample_names):
    """
    
    Args:
        sample_names (list): contains name of each sample
    Returns:
        headers (list): list of headers for every pairwise comparison between samples
                        ['samplename:samplename']
    """
    samples = sample_names[:-1]
    targets = sample_names[1:]
    headers = []
    for sample in samples:
        for target in targets:
            headers.append("%s:%s"%(sample, target))
        del targets[0]
    return headers


#loop through all proteins, perform calculations
def shift_score(dataframe, groups, sample_names):
    """
    determines hausdorff scores for given sample groups based on hausdorff distance df
    
    Args:
        dataframe (pd df): dataframe containing hausdorff distances between sample pairs
        groups (list): contains two lists with sample numbers belonging to each group
        sample_names (list): contains name of each sample
    Returns:
        final_frame (pd df): input dataframe with added column headers, and combined score
                             columns.
    """
    print "determining hausdorff scores..."
    between_comp = comp_between(groups)
    within_comp = comp_within(groups)
    subframe = dataframe

    # between score
    sub_between_comp = comp_between(groups)
    sub_between_scores = np.sum(abs(subframe[between_comp]), axis = 1)/len(between_comp)
    
    # within score
    sub_within_comp = comp_within(groups)
    sub_within_scores = np.sum(abs(subframe[within_comp]), axis = 1)/len(within_comp)
        
    # total score
    sub_total_scores = sub_between_scores/sub_within_scores
    
    #merge scores with shift frame
    new_subdata = [subframe, sub_between_scores, sub_within_scores, sub_total_scores]
    new_subframe =  pd.concat(new_subdata, axis = 1)
    
    scoresframe = new_subframe.ix[:, -3:]
    scoresframe.columns = ['Between', 'Within', 'Combined']

    final_frame = pd.merge(dataframe, scoresframe, left_index = True, right_index = True, how = 'outer')
    final_frame.sort_values('Combined', inplace = True, ascending = False)
    
    #replace 1:2 style comparison headers with sample names, after scores have been computed. (for user convenience)
    comparison_headers = get_sample_name_headers(sample_names)    # create headers list with sample name headers
    all_headers = list(final_frame.columns)                       # take headers from final frame
    all_headers[:len(comparison_headers)] = comparison_headers      
    final_frame.columns = all_headers                   # assign new sample name headers to corresponding comparison columns 

    return final_frame
        





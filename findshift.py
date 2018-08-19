"""
findshift module: calculates hausdorff scores between samples

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md
"""
# import statements
import pandas as pd
import numpy as np
import hausdorff

    
# functions
def strip_frame(dataframe, dataloc):
    """
    strips dataframe containing complexome profiling data from any extra columns
    
    Keyword arguments:
        dataframe -- pandas dataframe containing complexome profiling data
        dataloc   -- tuple containing first and last column header of complexome
                     profiling data
    Returns:
        stripped_frame -- pandas dataframe containing just columns with complexome profile
                          data
    """
    stripped_frame = dataframe.ix[:, dataloc[0]:dataloc[1]]
    return stripped_frame
    
def slicer(row, samplelengths):
    """
    splits dataframe row into seperate series per sample
    
    Keyword arguments:
        row -- list, row from df containing complexome profile data
        samplelengths -- list of numbers, length of each sample
    Returns:
        slices -- list, containing lists of values for each sample
    """
    slices = []
    start = 0
    for samplel in samplelengths:
        slices.append(row[start:start+samplel])
        start += samplel
    return slices
        
def get_headers(samplenum):
    """
    gets header numbers for all pairwise sample combinations for hausdorff distance data
    
    Keyword argument: samplenum -- number of samples
    Returns headers -- list of pairwise combinations of sample numbers in '1:2' format  
    """
    samples = range(1, samplenum+1)[:-1]
    targets = range(1, samplenum+1)[1:]
    headers = []
    for sample in samples:
        for target in targets:
            headers.append("%s:%s"%(sample, target))
        del targets[0]
    return headers
    
def get_normlist(dataframe, target_frame, normcolumn):
    """
    takes normcolumn from complexome profile frame, adds to hausdorff score frame
    
    Keyword arguments:
        dataframe -- pandas dataframe containing complexome profile data and normcol
        target_frame -- pd df, hausdorff score frame to which normcol will be added
        normcolumn -- string, header of normalisation column in complexome profile df
    Returns:
        new_frame -- pandas dataframe, hausdorff score frame with added norm column
    """
    # gets norm column from main dataframe, adds it to new dataframe, created by find_hausdorff function. 
    normlist = dataframe[normcolumn]
    new_frame = pd.concat([target_frame, normlist], axis = 1)
    return new_frame
    
def find_hausdorff(dataframe, samplelengths, mainframe, alignment, groups, normcolumn, norm_factor):
    """
    calculates hausdorff scores, stores results in new dataframe
    
    Keyword arguments:
        dataframe -- stripped pandas dataframe containing just complexome profiling data
        samplelengths -- list of numbers, number of slices for each sample
        mainframe -- pandas df, complexome data containing extra columns (normcol etc)
        alignment -- final alignment result, dictionary with aligned slice index list for
                     each sample
        groups -- list of 2 lists containing sample numbers for groups to be used for
                      calculation of intensity difference between sample groups
        normcolumn -- string, header of column from mainframe containing normalisation 
        norm_factor (float): factor that determines ratio between plain height and width
                             height = width (align length) * norm_factor
    Returns:
        dist_frame -- pandas dataframe containing:
                        - hausdorff distances between sample pairs for each protein
                        - normalisation column
                        - intesity differences: differences in protein abundance between 
                                                sample groups.
                                               (group2_intensity/group1_intensity)
    """
    print("determining hausdorff distances...")
    distdata = []
    intensity_diffs = []
    counter = 1
    for row in dataframe.itertuples():                     
        print("starting on new protein: #", counter)
        index = row[0]                                             
        data = list(row)[1:]
        slices = slicer(data, samplelengths)   # slice rows into seperate series per sample
        intensity_diffs.append(intensity_diff(slices, groups, alignment))  
        slices = hausdorff.square_series(slices, norm_factor)    # create square 2D plane
        hausdorffs = pairwise_hausdorff(slices) # calculate hausdorff distances
        distdata.append((index, hausdorffs))                               
        counter += 1
    samplenum = len(samplelengths)                       
    headers = get_headers(samplenum)                       
    dist_frame = pd.DataFrame.from_items(distdata, orient = 'index', columns = headers)
    if normcolumn != None:
        dist_frame = get_normlist(mainframe, dist_frame, normcolumn)
    dist_frame['intensity diffs'] = intensity_diffs                 
    dist_frame.index.name = dataframe.index.name                   
    return dist_frame 
            
def pairwise_hausdorff(slices):
    """
    determines hausdorff distance between pairwise combinations of given sample series
    
    Keyword arguments:
        slices -- list, containing list with intensity values for each protein
    Returns:
        hausdorffs -- list of hausdorff distances for each pairwise combination
    """
    samples = slices[:-1]
    targets = slices[1:]
    hausdorffs = []
    for sample in samples:
        for target in targets:
            hausdorffs.append(hausdorff.hausdorff(sample, target))
        del targets[0]
    return hausdorffs   


def gap_correct(alignments, slices):
    """
    replaces repeated intensity values through alignment with 0's
    
    repeated vals are replaced to calculate original protein intensities before alignment
    Keyword arguments:
        alignments -- final alignment result, dictionary with aligned slice index list for
                     each sample
        slices -- list, containing list with intensity values for each protein
    """
    # corrects for gaps that have been filled by alignment, replaces filled gap values with 0 (to avoid artificially increased intensity)
    for slice in range(len(slices)):            # loop through sample slices
        alignment = alignments[str(slice+1)]     # get alignment corresponding to slice from alignments
        for i in range(len(alignment)):          # loop through alignment values
            if i != 0:                                     
                if alignment[i] == alignment[i-1]:          # find gaps, if gap assign 0 to corresponding value in slice
                    slices[slice][i] = 0
    return slices

def intensity_diff(slices, groups, alignment):
    """
    calculates original intensity difference between sample groups
    
    Keyword arguments:
        alignments -- final alignment result, dictionary with aligned slice index list for
                     each sample
        slices -- list, containing list with intensity values for each protein
        groups -- list of 2 lists containing sample numbers for groups to be used for
                      calculation of intensity difference between sample groups
    Returns
        average_group1_intensity/average_group_2_intensity
        None if average intensity of both groups == 0
    """
    corrected_slices = gap_correct(alignment, slices)   # correct for extra intensity caused by repeated values through gaps
    avg_group1_intensity = np.mean([sum(corrected_slices[i-1]) for i in groups[0]])    # determine average sample intensity for group 1
    avg_group2_intensity = np.mean([sum(corrected_slices[i-1]) for i in groups[1]])    # determine average sample intensity for group 2

    if avg_group1_intensity != 0 and avg_group2_intensity != 0:            # check averages to avoid dividing by 0
        return avg_group1_intensity/avg_group2_intensity         
    else:
        return None 
    


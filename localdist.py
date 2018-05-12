"""
localdist module: calculates local distances between all complexome profile sample pairs

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md
"""
# import statements
import itertools

def gridprint(grid):
    """prints a grid, for debugging purposes"""
    row = 1
    for i in grid:
        print i
        row += 1

def distcalc(series1, series2):
    """
    calculates local distance at slice between samples based on all proteins.
    
    takes two lists of intensity values for all proteins at a certain slice
    calculates local distance by summing all differences in intensity per protein    
    Keyword arguments:
        series1/2 -- list of values. one value for each protein at given slice of sample
    Returns
        distance between samples at certain slice
    """
    dist = 0
    for first, second in itertools.izip(series1, series2):          # iterates trough both lists in paralel
        dist += abs(first - second)                                 # calculates absolute distance between corresponding values from 2 lists
    return dist    

def localdist(sample1, sample2):
    """
    takes protein series data from 2 samples, and returns a 2D list grid with local costs
    
    Keyword arguments:
        sample1/2 -- 2D list structure with all intensity values for each protein at each
                     slice for one sample
    Returns:
        localcost -- 2D list grid with local costs between two samples
    """
    # initialise grid
    localcost = []
    for column in range(len(sample1)): 
        row = []
        for point in range(len(sample2)):
            row.append(None)
        localcost.append(row)
    
    # loop that passes through all grid points 
    for y in range(len(sample1)):
        for x in range(len(sample2)):
            localcost[y][x] = distcalc(sample1[y], sample2[x])     # calculates local cost with two protein series using distcalc function
    
    return localcost

def pairwise_localdist(data):
    """
    takes complexome profile data of all samples, creates local cost grids for each pair
    
    Keyword argument:
        data -- 3D list structure with complexome profile data for all samples
    Returns:
        localdict -- dictionary with local distance grid for each sample pair
                   - keys 'samplenum:samplenum'  (ie: '1:2')
                   - values 2D list structure
    """
    localdict = {}
    samples = range(1,len(data)+1)             # the samples the main loop iterates through
    targets = range(1,len(data)+1)             # the target samples that get aligned with samples
    for sample in samples:                     #loops through each sample
        for target in targets:                 # per sample, loops through all targets
            localdict[str(sample) + ":" + str(target)] = localdist(data[sample-1], data[target-1])
        del targets[0]                         # removes first target after it has been aligned with all samples already
    return localdict    

"""
hausdorff module: computes hausdorff distances

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md
"""

# import statements
import math

# functions
def convert_set(set):
    """converts 1dimensional series to two dimensional data using list index"""
    new_set = []
    for i in range(len(set)):
        new_set.append((i,set[i]))
    return new_set
    
def distance(point, target):
    """ returns absolute distance between 2 points using (x,y) coordinates"""
    xdist = abs(point[0]-target[0])
    ydist = abs(point[1]-target[1])
    squareddist = xdist**2 + ydist**2
    dist = math.sqrt(squareddist)
    return dist

def hausdorff(set1, set2):
    """
    computes hausdorff distance between two sets of points
    
    Keyword arguments:
        set1/2 -- list of numeric values
    Returns:
        hausdorff distance between sets, based on euclidian distance between sets.
    """
    set1 = convert_set(set1)                    # convert 1D series into 2D points
    set2 = convert_set(set2)
    hausdorffs = []  
    for set in [(set1, set2), (set2, set1)]:                      # loop through 2 hausdorf directions: h(1,2) and h(2,1), store dist in hausdorffs
        nearest_distances = []
        for point in set[0]:    # loop through all points of set1
            distances = []  
            #determine distance from point to all targets, store in distances
            for target in set[1]:       
                distances.append(distance(point, target))
            #select lowest distance --> closest distance from point to set2, store in nearest
            nearest_distances.append(min(distances))
        hausdorffs.append(max(nearest_distances))       # select max hausdorff distance (which is the hausdorff distance between set1 and set2)
    return max(hausdorffs)

def square_series(slices, norm_factor):
    """
    normalises list of protein intensity values per sample
    
    sets maximum value across all samples to the length of the series, creating a square 
    2D plane in which to calculate hausdorff distances
    Keyword arguments:
        slices -- list, containing list of numeric values per sample
        norm_factor (float): factor that determines ratio between plain height and width
                             height = width (align length) * norm_factor
                                
    Returns:
        slices -- list, containing lists of normalised numeric values per sample
    """
    samplelength = len(slices[0])            # determine sample length (assuming same length for all slices)
    plain_height = samplelength * float(norm_factor)
    maximum = max([max(slice) for slice in slices])          # check whether row of data is not 0 everywhere
    if maximum != 0:
        factor = plain_height/float(maximum)                    # determining factor to multiply all values with
        new_slices= []              
        for slice in slices:                                    # applying factor to all values in slices
            new_slices.append([x*factor for x in slice])
        return new_slices
    else:
        return slices

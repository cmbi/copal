"""
hausdorff module: computes hausdorff distances

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md
"""

# import statements
import math
import numpy as np

# functions
def convert_set(series):
    """converts 1dimensional series to two dimensional np array """
    two_D_series = [[ix,val] for ix,val in enumerate(series)]
    two_D_array = np.array(two_D_series)
    return two_D_array

def min_euclid(point, series):
    """
    determine min euclid dist between point and every point of series

    Args:
        point: np array of length 2 with x,y coords
        series: 2D np array of shape(n,2) with x,y coords for each point
    Returns:
        closest: minimum euclid dist between point and series
    """
    diff = np.square(np.abs(series - point))
    summed = np.sum(diff, axis = 1)
    result = np.sqrt(summed)
    closest = np.min(result)
    return closest

def hausdorff_one_sided(template, target):
    """
    determine one-sided hausdorff distance from template to target
    """
    minima = np.apply_along_axis(min_euclid,1,template,target)
    max_min = minima.max()
    return max_min

def hausdorff(left, right):
    """
    determine hausdorff distance between two (list-like) series 
    """
    left = convert_set(left)
    right = convert_set(right)
    first_side = hausdorff_one_sided(left, right)
    second_side = hausdorff_one_sided(right,left)
    hausdorff = max(first_side,second_side)
    return hausdorff

def square_series(slices, norm_factor):
    """
    normalises list of protein intensity values per sample
    
    sets maximum value across samples to a factor of the series length,
    creating a 2D plane in which to calculate hausdorff distances
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

def new_square_series():
    pass

if __name__ == "__main__":
    pass
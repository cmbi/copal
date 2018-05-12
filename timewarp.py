"""
timewarp module: performs dynamic time warping between two samples or alignments

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md
"""
        
def distcalc(localdict, align1, align2,y,x):
    """
    calculates local distances between local alignments at certain slice coordinate
    
    Args:
        localdict: dict with 2D local cost value grid for each sample pair
        align1/2 (dict): samples as keys (int), aligned slice indices as values (list)
        x/y: indices of slice pair to calculate local distance between.
    Return: 
        realdist: mean local distance between two alignemnts at slice pair: x,y
    """
    dist = 0
    for key1 in align1.keys():
        for key2 in align2.keys():
            if int(key1) <= int(key2):
                dist += localdict[key1 + ':' + key2][align1[key1][y]][align2[key2][x]]
                
            else:
                dist += localdict[key2 + ':' + key1][align2[key2][x]][align1[key1][y]]
    comparisons = len(align1.keys())*len(align2.keys())
    realdist = dist/comparisons
    return realdist


def timewarp(localdict,align1,align2):
    """
    performs timewarping between two alignments, using local distances stored in localdict
    
    Args:
        localdict: dict with 2D local cost value grid for each sample pair
        align1/2 (dict): input alignments samples as keys (int), aligned slice indices as
                         values (list)
    Returns:
        new_align: resulting alignment, samples as keys (int), aligned slice indices 
                    as values (list)
        alignment_cost: global cost of alignment
    """
    # function that performs dynamic time warping, using local distance grids
    # creates grid datastructures for local dists and backtrace arrows
    # loops through all grid points
    # perform backtrace
    # warp starting series, store in new alignment dictionary 
    
    al1len = len(align1[align1.keys()[0]])                   # length of alignments
    al2len = len(align2[align2.keys()[0]])

    # initialise data structures. 2 grids, one with distance cost values, one with arrows for backtracing.
    distcost = []
    arrows = []
    for i in range(al1len):          #add rows to the 2 grids  
        distrow = []                    
        arrowrow = []
        for col in range(al2len):          #fill rows of both grids
            distrow.append(None)
            arrowrow.append((None,None))
        distcost.append(distrow)
        arrows.append(arrowrow)
    
    # loop that passes trough all grid points
    for y in range(al1len):
        for x in range(al2len):
            # calculate distance costs and arrows for all grid points, store in respective data structures
            # if at start point
            if y == 0 and x == 0:
                distcost[y][x] = distcalc(localdict,align1, align2,y,x)
                arrows[y][x] = (None,None)
            # elif at top edge
            elif y == 0:
                distcost[y][x] = distcost[y][(x-1)] + distcalc(localdict, align1, align2,y,x)
                arrows[y][x] = (y,x-1)
            # elif at left edge
            elif x == 0:
                distcost[y][x] = distcost[y-1][x] + distcalc(localdict, align1, align2,y,x)
                arrows[y][x] = (y-1,x)
            # else:  normal cases. dynamic programming takes place here. least costly point of origin is chosen
            else:
                lowestcost = distcost[y-1][x]             #
                arrows[y][x] = (y-1,x)
                if lowestcost > distcost[y][x-1]:
                    lowestcost = distcost[y][x-1]
                    arrows[y][x] = (y, x-1)
                if lowestcost > distcost[y-1][x-1]:
                    lowestcost = distcost[y-1][x-1]
                    arrows[y][x] = (y-1,x-1)
                distcost[y][x] = lowestcost + distcalc(localdict, align1, align2,y,x)


    # backtrace alignment, store in list    --> MOVE TO NEW FUNCTION
    alignment = []
    y = al1len-1                           #y,x --> outer corner of arrows and distcost grid. location of alignment cost and first arrow
    x = al2len-1                          
    alignment_cost = distcost[y][x]                # storing alignment cost
    alignment.append((y,x))                        # store first location in alignment
    while x != None or y != None:
        alignment.append(arrows[y][x])             # store next location in alignment
        y,x = arrows[y][x][0], arrows[y][x][1]     # change x,y values according to arrow tuple at former x,y coordinates in arrows
    del alignment[-1]                              # delete last entry of alignment, which is (None,None)
    alignment.reverse()                            # reverse order of backtrace, to make it a forward trace
    
    # warp old alignment, store in new alignment  --> MOVE TO NEW FUNCTION
    new_align = {}
    for key in align1.keys():                          # loop through all sequences in align1
        new_align[key] = []
        for i in alignment:                            # loop through all values of the forward trace alignment
            new_align[key].append(align1[key][i[0]])   # modify old sequences, acording to newly made alignment
    for key in align2.keys():                          # do the same for align 2 sequences
        new_align[key] = []
        for i in alignment:
            new_align[key].append(align2[key][i[1]])
    
    return (new_align, alignment_cost)
    


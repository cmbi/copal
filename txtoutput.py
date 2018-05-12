"""
txtoutput module: performs dynamic time warping between two samples or alignments

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md 
"""

# import statements
import pandas as pd
import numpy as np
import os

# functions
def show_gap(alignment):
    """
    visualizes alignment slice indices of one sample. X's indicate gaps.
    
    Args:
        alignment (list): list of aligned index values of 1 sample
    Returns: 
        strline (str): string of X's and _'s, displaying gap locations
    """
    line = []
    for i in range(len(alignment)):
        if i != 0:
            if alignment[i] == alignment[i-1]:
                line.append('X')
            else:
                line.append('_')
        else:
            line.append('_')
    strline = "".join(line)
    return strline

                                
def show_align(alignment, file):
    """
    outputs string visualisation of gap locations of given alignment
    
    outputs visualisation to stdout and writes to file object 
    Args:
        alignment (dict): given alignment, samples as keys (int), aligned slice indices 
                          as values (list)
        file: file object of output text file with align info
    Returns None
    """
    samplekeys = alignment.keys()
    samplekeys.sort()
    for key in samplekeys:
        put("sample {:4s}: ".format(key) + str(show_gap(alignment[key])), file)


def put(string, file):
    """
    prints string to stdout and writes to given file object
    
    Args:
        string(str): line to write and print
        file: file object to write to
    Returns None
    """
    print string
    file.write(string+ "\n")


def create_ranked_file(rnk_filename, dataframe, columns):
    #dataframe['Fasta headers'] = dataframe.index.values
    newframe = pd.concat([dataframe[columns[0]], dataframe[columns[1]]], axis = 1)
    newframe.dropna(inplace = True)
    np.savetxt(rnk_filename, newframe.values, fmt = '%s', delimiter = ' ')
    print "outputfile: ", rnk_filename
    print "outputfile length: ", len(newframe.index.values)
        
def text_output(input, output, align_info, output_path):
    """
    ouputs information on COPAL alignment analysis to stdout and writes to txt file
    
    outputs information on: 
        - input and output dataset size
        - normalisation process
        - length of samples before and after alignment
        - clustering step of progressive alignment
        - gap locations for each sample
        - groups used for hausdorff calculations
        - location where output is stored on drive
    
    Args:
        input (dict): containing all provided input for analysis
        output (dict): containing results produced during COPAL analysis
        align_info (file_object): file object to write output to
        output_path (str): file path where output folder is stored
    Returns None
    """     
    
    txtfile = open(align_info, 'w')
    print "-" * 60
    put("ALIGNMENT INFO\n", txtfile)
    put("source file(s): " + str(input['filename']), txtfile)
    
    # dataset size info
    matched_length = output['newdataframe'].shape[0]
    put("Length of original dataframe(s): ", txtfile)
    for i in output['original_frame_lengths']:
        put(str(i), txtfile)
    if len(output['original_frame_lengths']) > 1:
        put("length of matched datasets: " + str(matched_length), txtfile)

    # normalisation info
    if output['mod_factors'] != None:
        mod_list = []
        for i in range(len(output['mod_factors'])):
            mod_list.append(input['samplenames'][i] + ": " + str(output['mod_factors'][i]))
        put("normalisation factors:", txtfile)      
        for i in mod_list:
            put(i, txtfile)
    else:
        put("no normalisation performed.", txtfile)
    
    # sample and alignment length info
    put("samplelengths: " +  str(output['samplelengths']), txtfile)
    put("Length of final alignment: " + str(output['final_align_length']), txtfile)

    put("Samples used:\n#  sample name", txtfile)
    for i in range(len(input['samplenames'])):
        put("%d  %s"%(i+1,input['samplenames'][i]), txtfile)

    # clustering info
    paircostkeys = output['pairwisecosts'].keys()
    paircostkeys.sort()
    put("\npairwise clustering results (from low to high cost):", txtfile)
    for key in paircostkeys:
        put("pair: " + str(output['pairwisecosts'][key]) + "   cost: " + str(key), txtfile) 

    put("\npairwise matches ordered by cost scores from low to high:", txtfile)
    sortedpairs = []
    for i in paircostkeys:
        sortedpairs.append(output['pairwisecosts'][i])
    put(str(sortedpairs), txtfile)  
    put("", txtfile)
    put("Alignment order: " +  str(output['multiple_alignment_order']), txtfile)
    put("", txtfile)
    
    # gap location info
    put("visualision of gap locations of final alignment: ", txtfile)
    put("-----------------------------------------------------------\n", txtfile)
    show_align(output['final_alignment'], txtfile)
    put("\n-----------------------------------------------------------", txtfile)
    if input['hausdorff_scoring']:
        put("groups used for combined scores: " +str(input['groups']), txtfile)

    # close file and return to starting directory after storing file in results folder
    txtfile.close()
    os.chdir('..')
    print "Complexome profile alignment Complete! Your output can be found in: ", output_path
    
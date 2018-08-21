""" 
module that complements complexome profiling datasets so that proteins match 

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md

"""

# import statements
import pandas as pd

# functions
def complement_dataframes(dataframes, normcol, gseacol, samplecolumns):
    """
    takes list of dataframes, complements dataframes that are missing proteins
    
    Adds rows for missing proteins with profiling data and columns required for analysis
    Keyword arguments:
        dataframes -- list of pandas dataframes with complexome profiling data
        normcol    -- string, name of column header specifying proteins to be used for
                      normalisation
        gseacol    -- string, name of column header containing identifiers for .rnk output
        samplecolumns -- list of strings, containing first and last column header for
                         each sample 
    Returns:
        list of complemented pandas dataframes
    """
    if len(dataframes) > 1:
        print( "complementing samples that are missing proteins...")
        matcheddataframes = []
        datamatches = dataframes[1:]
        target_frame = dataframes[0]
        # get important columns to be transferred to target_frame for each new protein
        identcolumns = []
        if normcol != None:
            identcolumns.append(normcol)
        if gseacol != None:
            identcolumns.append(gseacol)
            
        #target_frame (first dataframe) samplecolumns
        target_samplecolumns = samplecolumns[0]
        
        # loop through all datamatches( all frames except first.)
        # add any proteins from match frames missing in target_frame to target_frame
        for match in datamatches:
            target_frame = match_two_frames(target_frame, match, identcolumns, target_samplecolumns)
        matcheddataframes.append(target_frame)
        
        # loop through all datamatches again.
        # add proteins not in match frames to these frames, so all dataframes have all proteins.
        for i in range(len(datamatches)):
            match_columns = samplecolumns[i+1]
            new_match_frame = match_two_frames(datamatches[i], matcheddataframes[0], identcolumns, match_columns)
            new_match_frame = new_match_frame.loc[matcheddataframes[0].index.values]
            matcheddataframes.append(new_match_frame)
        return matcheddataframes
    else:
        return dataframes

def find_missing(target_frame, match_frame):
    """
    finds proteins missing from target frame, present in match_frame
    
    Keyword arguments:
        target_frame -- pandas dataframe of which missing proteins will be determined 
        match_frame -- pandas dataframe with proteins to be checked for
    Returns:
        list of strings, identifiers of missing proteins
    """
    series1 = list(target_frame.index.values)
    series2 = list(match_frame.index.values)
    return [i for i in series2 if i not in series1]
        
def match_two_frames(target_frame, match_frame, identcolumns, target_samplecolumns):
    """
    complements protein rows of two target frame to contain all proteins from match_frame
    
    Keyword arguments:
        target frame -- pandas dataframe with complexome profiling data
        match_frame  -- pandas dataframe with complexome profiling data
        identcolumns -- list of strings, names of relevant column headers
        target_samplecolumns -- samplecolumns of target_frame
    Returns:
        target_frame -- pandas dataframe, complemented with proteins from match_frame
    """
    # get list of missing proteins
    missing_list = find_missing(target_frame, match_frame)
    # for each missing protein: 
    for protein in missing_list:
        # get identifiers with appropriate column headers
        protein_data = []
        protein_indexes = []
        for i in identcolumns:
            protein_data.append(match_frame.ix[protein, i])
            protein_indexes.append(i)
        
        # get sample data headers
        #target frame headers
        target_headers = list(target_frame.columns)        # get all headers from target frame
        data_headers = []
        for sample in target_samplecolumns:      # go trough all samples in target frame samplecolumns
            first_index = target_headers.index(sample[0])
            last_index = target_headers.index(sample[1])+1
            
            # get the headers belonging to this sample, add to data headers
            data_headers += target_headers[first_index:last_index]           
        
        data_zeros = [0]*len(data_headers)  #create list of 0's to fill protein data rows
        protein_data += data_zeros        # add to protein_data list
        protein_indexes += data_headers   # add data headers to indexes list
        
        # create pandas series with new data and headers as indexes
        new_row = pd.Series(protein_data, index = protein_indexes, name = protein)
        # add series to target_frame as a new row
        target_frame = target_frame.append(new_row)
    print( "length target_frame: ", len(target_frame))
        
    return target_frame

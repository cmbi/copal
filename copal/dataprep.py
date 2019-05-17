"""
dataprep module: contains functions that prepare data for COPAL alignment

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md
"""

# import statements
from . import complement_frames as complement
import pandas as pd

def input_checker(dataframes, input):
    """
    Checks if provided input matches loaded data, raises ValueError if otherwise

    checks:
        - norm column
        - sample columns
        - gsea column

    Args:
        dataframes (list): list of pd dfs with loaded input datasets
        input (dict): input dictionary containing all user provided input
    Returns None
    """
    # check if scoring is checked without alignment
    if input['align_check'] == False and input['hausdorff_scoring'] == True:
        raise ValueError("Cannot do hausdorff scoring without alignment!")

    for ix,df in enumerate(dataframes):
        datacols = df.columns
        # test if all input samplecolumns match with data
        for samplecol in input['samplecolumns'][ix]:
            left_scol, right_scol = samplecol
            if not left_scol in datacols:
                raise ValueError("Sample column not found: '{}'".format(left_scol))
            if not right_scol in datacols:
                raise ValueError("Sample column not found: '{}'".format(right_scol))

        # test if input normcol matches with data
        if input['normcol']:
            if not input['normcol'] in datacols:
                raise ValueError("normalisation column not found: '{}'".format(input['normcol']))

        # test if gsea column matches with data
        if input['gsea_output'] and input['GSEA_rank_column']:
            if not input['GSEA_rank_column'] in datacols:
                raise ValueError("GSEA output column not found: '{}'".format(input['GSEA_rank_column']))

# functions
def multi_dataload(identifier, filename, sheet, skip_rows, file_type):
    """
    parses excel or text-based dataset, loads into pandas dataframe

    Keyword arguments:
        identifier -- string, header of column containing protein identifiers
        filename   -- string, path/filename of data file
        sheet      -- string, name of sheet containing data (in case of excel file)
        skip_rows  -- int, number of header rows to skip (should still have 1 header row)
        file_type  -- info containing info on type of file.
                    - 'excel' if excel file.
                    - ('column seperator','decimal seperator') if text based file
        Returns:
            pandas dataframe, with protein identifiers as index column
    """
    if file_type == 'excel':
        data = pd.read_excel(filename, sheet_name = sheet, skiprows = skip_rows)
    else:
        data = pd.read_csv(filename, skiprows = skip_rows, sep = file_type[0], decimal = file_type[1])
    if not identifier in data.columns:
        raise ValueError('given protein ID not a column header! check input or file format')
    data.dropna(subset = [identifier], inplace = True)
    data[identifier] = data[identifier].astype(str)
    data = data.set_index(identifier, drop = True)
    return data

def get_norm_list(filen):
    """
    gets list of protein IDs from plain text file with one ID on each line.

    Args:
        filen -- path/filename of text file containing protein IDs to be used for
                 normalisation
    Returns:
        list of strings, containing protein IDs to be used for normalisation
    """
    with open(filen, 'r') as file_object:
        return [line.strip() for line in file_object]

def is_norm(id,normlist):
    """
    for row from pandas dataframe, determines if identifier is in normlist

    Keyword arguments:
        row -- pandas series, row from pandas df with complecome profiling data
        column -- string, header of dataframe column containing protein identifier
        normlist -- list of strings,
    """
    if id in normlist:
        return 1
    else:
        return 0

def get_norm_data(normfile, compare_column, dataframe):
    """
    loads list of protein IDs to be used for normalisation, adds bool column to profile df

    Keyword arguments:
        normlist -- list of strings, protein IDs to be used for normalisation
        compare_column -- string, header of column with protein IDs to compare to

    Returns:
        input pandas df, with extra boolean column specifying presence in normlist
    """
    norm_list = get_norm_list(normfile)
    dataframe['norm'] = dataframe.apply\
                        (lambda x: is_norm(x.name, norm_list), axis = 1)
    return dataframe

def varied_length_extract(data, samplecolumns):
    """
    gets profilng data from pandas dataframe structure, loads into 3D structure

    Keyword arguments:
        data -- pandas dataframe containing complexome profiling data
        samplecolumns -- list of tuples, with header names of first and last column for
                         each sample
    Returns:
        totaldata -- 3D list strucure with numeric values. profiling data ready for
                     alignment
                   - outer list with one list object for each sample
                   - each sample is 2D array of lists containing slices and protein rows
    """
    totaldata = []
    for s in samplecolumns:
        if not s[0] in data.columns or not s[1] in data.columns:
            raise ValueError("Samplecolumn not found: '{}' or '{}'".format(s[0],s[1]))
        sample = []
        sampledata = data.ix[:, s[0]:s[1]]
        for col in sampledata.keys():
            slice = [float(i) for i in data.ix[:,col]]
            sample.append(slice)
        totaldata.append(sample)
    return totaldata

def match_frames(dataframes):
    """
    matches list of pd dfs, to contain only shared proteins in the same order

    Keyword arguments:
        dataframes -- list of pandas dataframes
    Returns:
        list of pd dfs, matched and ordered based on index column
    """
    if len(dataframes) > 1:
        print( "creating subsets of data with matched and ordered proteins...")
        matcheddataframes = []
        datamatches = dataframes[1:]
        matchframe = dataframes[0]
        for i in datamatches:
            matchframe = matchframe.loc[i.index.values]
            matchframe.dropna(inplace = True, how = 'all')
        matcheddataframes.append(matchframe)
        for i in datamatches:
            matcheddataframes.append(i.loc[matcheddataframes[0].index.values])
        return matcheddataframes
    else:
        return dataframes

def gel_mod(totaldata, normalisationdata):
    """
    normalises total data, correcting for intensity differences in normalisation data

    Keyword arguments:
        totaldata -- 3D list structure containing numeric profiling data
        normalisationdata -- 3D list structure containing subset of totaldata

    Returns (modifieddata, mod_factors)
        modifieddata -- 3D list structure, normalised data based on intensity differences
                        in normalisation data
                       - normalised so that samples in normalisation data have equal
                         intensity
        mod_factors -- list of numeric values, normalisation factor used for each sample
    """
    #calculate total intensity of each gel using normalisationdata
    intensities = []
    for sample in normalisationdata:
        summed_sample = 0
        for slice in sample:
            summed_sample += sum(slice)
        intensities.append(summed_sample)
        if summed_sample == 0:
            raise ValueError('one sample did not detect any of the normalisation proteins')
    #determine modifications
    mod_factors = []
    for sample in intensities:
        mod_factors.append(intensities[0]/sample)
    #apply modifications do totaldata
    modifieddata = []
    samplenum = 0
    for sample in totaldata:
        new_sample = []
        for slice in sample:
            new_slice = []
            for value in slice:
                new_value = value * mod_factors[samplenum]
                new_slice.append(new_value)
            new_sample.append(new_slice)
        modifieddata.append(new_sample)
        samplenum += 1

    return (modifieddata, mod_factors)

def prep_data(dataframes, samplecolumns, normcol, normalisation, gseacol):
    """
    preforms data preperation required before complexome alignment

    matches datasets so they have matching proteins in the same order
    normalises data based on list of proteins to be used for normalisation
    Keyword arguments:
        dataframes -- list of pandas dataframes for each provided dataset
        samplecolumns -- list of tuples, with column headers for first and last column
                         containing data for each sample
        normcol -- string, header of boolean column specyfing proteins to be used for
                   normalisation
        normalisation -- boolean, whether normalisation is to be performed
        gseacol -- string, header of column containing IDs for GSEA .rnk output
    Returns (matcheddataframes, samplelengths, mofifieddata, mod_factors):
        matcheddataframes -- list of dataframes with matched and ordered datasets
        samplelengths -- list of numbers, number of slices of each sample
        modifieddata -- 3D list structure containing normalised profiling data numbers
        mod_factors -- list of numeric values, containing normalisation factors used per
                       sample
    """
    # in case of multiple dfs: complement dfs with missing proteins, and reorder frames
    matcheddataframes = complement.complement_dataframes\
                                        (dataframes, normcol, gseacol, samplecolumns)

    # extract totaldata  from dataframe(s)
    totaldata = []
    for i in range(len(matcheddataframes)):
        totaldata += varied_length_extract(matcheddataframes[i], samplecolumns[i])

    # create subsets of dataframes for normalisation using normcol variable
    if normcol != None:
        print( "creating subsets for normalisation...")
        normalisationdataframes = []
        for i in matcheddataframes:
            if normcol not in i.columns:
                raise ValueError("normalisation column name not found: '{}'".format(normcol))
            normalisationdataframes.append(i.loc[i[normcol] == 1.0])

        # only take common norm proteins to prevent skewed sample intensities.
        normalisationdataframes = match_frames(normalisationdataframes)

        # extract normalisation data from normalisation dataframes
        normalisationdata = []
        for i in range(len(normalisationdataframes)):
            normalisationdata += varied_length_extract(normalisationdataframes[i],
                                                       samplecolumns[i])
    else:
        normalisationdata = totaldata

    # calculate length of samples, store lengths in list
    samplelengths = []
    for i in totaldata:
        samplelengths.append(len(i))

    # perform normalisation
    if normalisation:
        # correct for varying gel intensities
        print( "normalising data to correct for varying gel intensities...")
        # select subset of normalisation proteins
        normalisation_results = gel_mod(totaldata, normalisationdata)
        modifieddata = normalisation_results[0]
        mod_factors = normalisation_results[1]
    else:
        modifieddata = totaldata
        mod_factors = None

    return (matcheddataframes, samplelengths, modifieddata, mod_factors)

if __name__ == "__main__":
    pass

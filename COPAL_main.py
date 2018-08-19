""" 
main COPAL module that runs analysis without GUI wrapper

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md
"""
# import statements
import pandas as pd
import dataprep, os
import localdist  
import multipletimewarp as msa_warp
import findshift
import shiftscore
import datatoexcel
import txtoutput as out

# functions
# ---------------------- input Data processing ---------------------------#
def input_processing(input):
    """
    performs input processing. loading data, combining datasets, normalisation.
    
    Keyword arguments:
        input -- dict, containing all provided input for analysis
    Returns:
        tuple -- containing:
                  - matcheddataframes: list of pandas dataframes after matching proteins
                  - samplelengths: list  of ints, number of slices for each sample
                  - mofifieddata: normalised iBAQ data in list in list in list structure
                  - mod_factors: list of numerics, factors for normalisation per sample
                  - original_frame_lengths: nr of rows/proteins before combining datasets
    """
    samplenum = len(input['samplenames'])

    #LOAD DATA FROM EXCEL
    print( "loading data from input files...")
    dataframes = []
    for i in range(len(input['filename'])):
        dataframes.append(dataprep.multi_dataload(input['identifier'][i],
                          input['filename'][i], input['sheetname'][i],
                          input['skiprows'][i], input['input_type'][i]))
    
    #CHECK IF INPUT MATCHES LOADED DATA
    dataprep.input_checker(dataframes, input)
                          
    #IF NESCESSARY: ADD NORMALISATION COLUMN
    if input['normfile'] != None:          # check if normlist value is not None.
        print( "loading normalisation IDs from file..."  )
        new_dataframes = []            # use normlist variable to add normcolumn to dfs
        for i in range(len(dataframes)):
            new_dataframes.append(dataprep.get_norm_data(input['normfile'],
                                  input['identifier'][i], dataframes[i]))
        dataframes = new_dataframes
        # change normcolumn variable to name of column that was just created.
        input['normcol'] = 'norm'                  
    
    # PREPARE DATA FOR TIME WARPING
    prepareddata = dataprep.prep_data(dataframes, input['samplecolumns'],
                                      input['normcol'], input['norm_check'],
                                      input['GSEA_rank_column'])
                                      
    matcheddataframes = prepareddata[0]
    samplelengths = prepareddata[1]
    modifieddata = prepareddata[2]
    mod_factors = prepareddata[3]
    original_frame_lengths = []
    print( "length of original dataframes:")
    for i in dataframes:
        print( len(i.index.values))
        original_frame_lengths.append(len(i.index.values))
    print( "length of matched dataframes:")
    for i in matcheddataframes:
        print( len(i.index.values))
    return (matcheddataframes, samplelengths, modifieddata, mod_factors, original_frame_lengths) 


# ---------------------- Complexome profile alignment ---------------------------#
def complexome_alignment(template_df, samplelengths, modifieddata, input):
    """
    performs complexome alignment
    
    Keyword arguments:
        template_df   -- pandas dataframe, contains all data provided in first dataset.
                         aligned data is inserted into this template dataframe.
        samplelengths -- list of ints, number of slices for each sample
        modifieddata  -- normalised iBAQ data, list in list in list structure
        input         -- dict, containing all provided input for analysis
    Returns:
        tuple containing:
            final_alignment -- dictionary with aligned slice index list for each sample
                            -- format: {'samplenumber':[slice indices]}
            multiple_alignment_order  -- list of tuples with alignment order
                                      -- format:[((align1 samples),(align2 samples), cost)]
            final_align_length -- integer, number of slices after alignment
            dataloc -- tuple with 2 strings, name of first and last columns containing data
            newdataframe -- pandas df, containing aligned data and columns from template
            pairwisecosts -- dictionary containing costs of all pairwise alignments
                          -- format: {'sample1:sample2':cost(type:float)}
    """
    # CALCULATE LOCAL DISTANCE GRIDS FOR ALL PAIRWISE ALIGNMENTS
    print( "calculating local distances...")
    localdict = localdist.pairwise_localdist(modifieddata)
    #print localdict.keys()

    # PERFORM PAIRWISE ALIGNMENTS
    print( "performing pairwise alignments...")
    # iterate through all pairwise combinations of samples
    pairwise_warp = msa_warp.pairwise_timewarp(modifieddata,samplelengths,localdict)
    pairwisedict = pairwise_warp[0]
    pairwisecosts = pairwise_warp[1]
    paircostkeys = pairwisecosts.keys()
    
    # DETERMINE ALINGMENT ORDER FOR MSA, PERFORM ALIGNMENTS
    print( "performing multiple alignments..."   )
    multiple_alignment = msa_warp.multiple_timewarp(pairwisecosts,samplelengths,localdict,pairwisedict)
    final_alignment = multiple_alignment[0]
    multiple_alignment_order = multiple_alignment[1]
    final_align_length = len(final_alignment[final_alignment.keys()[0]])
    samplenum = len(samplelengths)
    align_lengths = [final_align_length] * samplenum    # required as input for some functions, in order for the functions to also be able to handle varying sample lengths

    # WARP DATA USING FINAL ALIGNMENT
    print( "warping data using final multiple alignment...")
    warpeddata = msa_warp.datawarp(modifieddata, final_alignment)
    #transform warped data to pandas dataframe
    warpeddataframe = datatoexcel.datatoframe(warpeddata, input['samplenames'], align_lengths)

    # REPLACE ORIGINAL DATA WITH WARPED DATA IN TEMPLATE DATAFRAME
    newdataframe = datatoexcel.integratedata(template_df, warpeddataframe, input['samplecolumns'][0], input['identifier'][0])
    dataloc = (input['samplenames'][0] + "_1", input['samplenames'][-1] + "_" + str(align_lengths[-1]))   # headers of first and last data containing columns

    return (final_alignment, multiple_alignment_order, final_align_length, dataloc, newdataframe, pairwisecosts)
    
    
# ---------------------- Hausdorff score calculation ---------------------------#
def hausdorff_scoring(input, final_alignment, final_align_length, dataloc, newdataframe):
    """
    performs hausdorff scoring analysis on aligned dataframe
    
    Keyword arguments:
        input -- dict, containing all provided input for analysis
        final_alignment -- dictionary with aligned slice index list for each sample
                        -- format: {'samplenumber':[slice indices]}
        final_align_length -- integer, number of slices after alignment
        dataloc -- tuple with 2 strings, name of first and last columns containing data
        newdataframe -- pandas df, containing aligned data and columns from template
    Returns:
        score_frame -- pandas df containing columns with hausdorff scores
        newdataframe -- pandas df containing aligned data and hausdorff effect sizes
    """
    samplenum = len(input['samplenames'])
    align_lengths = [final_align_length] * samplenum
    
    # determine hausdorff distances and combined hausdorff scores
    nakedframe = findshift.strip_frame(newdataframe, dataloc)
    shift_frame = findshift.find_hausdorff(nakedframe, align_lengths, newdataframe, final_alignment, input['groups'], input['normcol'], input['hausd_factor'])
    score_frame = shiftscore.shift_score(shift_frame, input['groups'],input['samplenames'])

    # ADD COMBINED SCORE AND INT DIFFS TO MAIN DATAFRAME
    combined_score = score_frame['Combined']
    intensity_diffs = score_frame['intensity diffs']
    newdataframe = pd.concat([newdataframe, intensity_diffs, combined_score], axis = 1)
    newdataframe.index.name = score_frame.index.name
    
    return (score_frame, newdataframe)      

# ---------------------- Output results  ---------------------------#
def output_results(input, output):
    """
    produces output files of results and prints information to stdout
    
    produces results folder containing
        - excel file with aligned profiles, hausdorff scores, migration pattern graphs 
        - csv file of aligned profiles
        - csv file of hausdorff scores
        - text file with info on alignment process
        - .rnk file: ranked list with hausdorff effect size values for GSEA 
        - info on alignment process to stdout
    Keyword arguments:
        input -- dict, containing all provided input for analysis
        output -- dict, containing results produced during COPAL analysis
    Returns:
        None
    """ 
    samplenum = len(input['samplenames'])
    align_lengths = [output['final_align_length']] * samplenum
    
    # PREPARE FOR OUTPUT 
    # get file names for output files
    excel_output = input['analysis_name'] + ".xlsx"
    csv_score_frame = input['analysis_name'] + "_score_frame.csv"
    csv_dataframe = input['analysis_name'] + "_dataframe.csv"
    align_info = input['analysis_name'] + "_align_info.txt"

    # create output folder, change current directory to output folder
    cwd = os.getcwd()
    output_path = cwd + '/' + input['analysis_name'] + "_results"
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    os.chdir(output_path)

    # CSV OUTPUT OF DATAFRAMES (FOR POST-ANALYSIS)
    print( "exporting data to csv...")
    if input['hausdorff_scoring']:
        output['score_frame'].to_csv(csv_score_frame)
    output['newdataframe'].to_csv(csv_dataframe)

    # EXCEL OUTPUT OF DATAFRAMES, GRAPHS OF PROTEINS
    print( "exporting data to excel...")
    datatoexcel.exceldump(output['newdataframe'], excel_output, align_lengths, input['samplenames'], output['score_frame'], output['dataloc'])

    # .RNK OUTPUT FOR GENE SET ENRICHMENT ANALYSIS
    if input['hausdorff_scoring'] and input['gsea_output']:
        print( "generating .rnk files...")
        # add identifier column for .rnk file to score frame
        rank_ident = output['newdataframe'][input['GSEA_rank_column']]
        rank_frame = pd.concat([output['score_frame'], rank_ident], axis = 1)

        # create ranked file with combined hausdorff score 
        out.create_ranked_file(input['analysis_name'] +'_comb_hausd_ranked_file.rnk', rank_frame, (input['GSEA_rank_column'], 'Combined'))
        
    # PROVIDE ALIGNMENT INFO TEXT OUTPUT (to stdout and .txt file)
    out.text_output(input, output, align_info, output_path)
    
# ---------------------- main function ---------------------------#
def main(input = input):
    """
    main function that runs COPAL analysis
    
    Keyword arguments:
        input -- dictionary containing required input for COPAL analysis
            'analysis name': string, name for current analysis, used in output files
            'filename': list of strings, path/name of files containing complexome profiles  
            'skiprows': list of numbers, number of header rows to skip for each input file
            'input_type': list, specifying format of each input file.
                        - if excel: 'excel'
                        - if text: ('column seperator','decimal seperator') 
            'samplecolumns': list of tuples with headers of first and last columns
                             containing profiling data for each sample
            'identifier': list of strings, header of column identifying proteins for each
                          file
            'sheetname': list of strings, names of excel sheets containing data for each
                         file. empty strings if not excel files
            'samplenames': list of strings, name for each sample
            'norm_check': boolean, specifying if normalisation is to be performed
            'normcol': string, header of column indicating proteins to be used for
                       normalisation
            'normfile': string: path/file of text file containing identifiers of proteins 
                        to be used for normalisation
            'hausdorff_scoring': boolean or 1/0, specifying whether to perform hausdorff
                                 scoring
            'gsea_output': boolean, specifying whether to produce .rnk file output for
                           GSEA analysis
            'GSEA_rank_column': string, header of column to be used as identifier in .rnk 
                                output for GSEA analysis
            'groups': list of 2 lists containing sample numbers for groups to be used for
                      hausdorff effect size calculation
            'hausd_factor' (float): ratio of height to width of plain in which hausdorff
                                    distances are calculated
    Returns None
    """
    output = {}
    
    processed_input = input_processing(input)
    matched_dataframes = processed_input[0]
    output['samplelengths'] = processed_input[1]
    modifieddata = processed_input[2]
    output['mod_factors']  = processed_input[3]
    output['original_frame_lengths'] = processed_input[4]

    align_result = complexome_alignment(matched_dataframes[0], output['samplelengths'], modifieddata, input)
    output['final_alignment'] = align_result[0]
    output['multiple_alignment_order'] = align_result[1]
    output['final_align_length'] = align_result[2]
    output['dataloc'] = align_result[3]
    output['newdataframe'] = align_result[4]
    output['pairwisecosts'] = align_result[5]
        
    if input['hausdorff_scoring']:
        hausdorff_result = hausdorff_scoring(input, output['final_alignment'], output['final_align_length'],
                                             output['dataloc'], output['newdataframe'])
        output['score_frame'] = hausdorff_result[0]
        output['newdataframe'] = hausdorff_result[1]
    else:
        output['score_frame'] = None
        
    output_results(input, output)
    

if __name__ == "__main__":


    # input dict values for analysis run with test files
    input = {}
    
    input['analysis_name'] = 'test_run'
    input['filename'] = ['Gel1_unaligned_samples_500_rows.csv', 'Gel2_unaligned_samples_500_rows.csv', 'Gel3_unaligned_samples_500_rows.csv']
    input['skiprows'] = [0, 0, 0]
    input['input_type'] = [(',', '.'), (',', '.'), (',', '.')]
    input['samplecolumns'] = [[(u'S1_G1_01', u'S1_G1_60'), (u'S2_G1_01', u'S2_G1_60'), (u'S3_G1_01', u'S3_G1_60'), (u'S4_G1_01', u'S4_G1_60'), (u'S5_G1_01', u'S5_G1_60')], [(u'S1_G2_01', u'S1_G2_60'), (u'S2_G2_01', u'S2_G2_60'), (u'S3_G2_01', u'S3_G2_60')], [(u'S1_G3_01', u'S1_G3_60')]]
    input['identifier'] = ['Gi', 'Gi', 'Gi']
    input['sheetname'] = ['', '', '']
    input['samplenames'] = [u'S1_G1', u'S2_G1', u'S3_G1', u'S4_G1', u'S5_G1', u'S1_G2', u'S2_G2', u'S3_G2', u'S1_G3']
    input['norm_check'] = True
    input['normcol'] = None
    input['normfile'] = None
    input['hausdorff_scoring'] = 1
    input['gsea_output'] = 1
    input['GSEA_rank_column'] = 'Gene Symbol' 
    input['groups'] = [[1, 2, 3, 4, 5], [6, 7, 8, 9]]
    input['hausd_factor'] = 1.0

    # main function call, starting COPAL analysis
    main(input = input)
    
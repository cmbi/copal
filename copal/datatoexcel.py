"""
datatoexcel module: contains function that export data to excel, output helper functions
author: Joeri van Strien

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md
"""
# import statements
import numpy as np
import pandas as pd

# functions
def get_headers(align_lengths, sample_names):
    """
    gets column headers for aligned complexome profiles
    
    Keyword arguments:
        align_lengths -- list of numbers, nr of slices for each sample (after alignment) 
        sample_names  -- list of strings, name of all samples
    Returns
        headers -- list of strings, column headers of aligned data columns
    """
    headers = []
    for i in range(len(align_lengths)):
        sample_headers = []
        for j in range(1, align_lengths[i]+1):
            sample_headers.append(sample_names[i] + "_" + str(j))
        headers = headers + sample_headers
    return headers

def datatoframe(totaldata, sample_names, align_lengths):
    """
    converts data from 3D list structure to pandas dataframe
    
    Keyword arguments:
        totaldata -- numeric data in 3D list structure
        sample_names -- list of strings, names of samples
        align_lenghts -- list of numbers, lengths of samples
    Returns:
        pandas dataframe, containing data with column headers based on sample names
    """
    headers = get_headers(align_lengths, sample_names)
    dataframes = []
    for i in totaldata:
        df = pd.DataFrame(data = i)
        df = df.transpose()
        dataframes.append(df)
    dataframe = pd.concat(dataframes, axis = 1)
    dataframe.columns = headers
    return dataframe
    
def integratedata(totaldataframe, dataframe, samplecolumns, identifier):
    """
    inserts df with aligned data into template dataframe with protein information columns
    
    Keyword arguments:
        totaldataframe -- dataframe containing aligned profiling data
        dataframe -- template dataframe with original extra data columns
        samplecolumns -- list of tuples, first and last column headers for each sample
        identifier -- string, header of column containing protein identifiers
    Returns
        pandas dataframe, containing aligned data and original extra protein info columns
    """
    datastartcol = totaldataframe.columns.get_loc(samplecolumns[0][0])
    dataendcol = totaldataframe.columns.get_loc(samplecolumns[-1][-1])
    beforeframe = totaldataframe.ix[:,:datastartcol]
    beforeframe.reset_index(drop = False, inplace = True)
    afterframe = totaldataframe.ix[:,dataendcol+1:]
    afterframe.reset_index(drop = True, inplace = True)
    newdataframe = pd.concat([beforeframe, dataframe, afterframe] ,axis = 1)
    newdataframe.set_index(identifier, inplace = True, drop = True)
    return newdataframe

def exceldump(dataframe,filename, samplelengths, sample_names, shift_frame, dataloc,
              normframe = None, normloc = None):
    """
    takes aligned complexome profile data and score data, outputs to excel
    
    first sheet will contain aligned complexome profiles
    second sheet will contain hausdorff scoring results
    third sheet will contain migration pattern graphs of the first 200 rows
    Keyword arguments:
        dataframe -- pandas df, aligned complexome profiles
        filename -- string, name of output excel file
        samplelengths -- list of numbers, number of slices per sample
        sample_names -- list of strings, names of samples
        shift_frame -- pandas df, containing hausdorff score results
        dataloc -- tuple, containing header of (first,last) col of complexome profile data
    Returns None
    """
    writer = pd.ExcelWriter(filename, engine = 'xlsxwriter')
    workbook = writer.book
    # dump aligned dataframe into excel
    dataframe.to_excel(writer, sheet_name = 'data')
    # apply conditional formatting to data columns
    datastartcol = dataframe.columns.get_loc(dataloc[0]) + 1
    dataendcol = dataframe.columns.get_loc(dataloc[1]) + 1
    worksheet = writer.sheets['data']
    apply_cond_format(dataframe,datastartcol,dataendcol,writer,worksheet,workbook)

    # set column width of data columns
    worksheet.set_column(datastartcol,dataendcol, 0.4)
    # set headers row as freeze pane
    worksheet.freeze_panes(1, 0)
    
    if normframe is not None and normloc is not None: 
        # dump normalised data into excel
        normframe.to_excel(writer, sheet_name='normed_not_aligned')
        
        # add formatting
        normstartcol = normframe.columns.get_loc(normloc[0]) + 1
        normendcol = normframe.columns.get_loc(normloc[1]) + 1
        normsheet = writer.sheets['normed_not_aligned']
        apply_cond_format(normframe,normstartcol,normendcol,writer,normsheet,workbook)
        normsheet.set_column(normstartcol,normendcol,0.4)
        normsheet.freeze_panes(1,0)

    #add shift_frame to seperate sheet if score analysis was performed
    if isinstance(shift_frame, pd.DataFrame):
        shift_frame.to_excel(writer, sheet_name = 'scores')
        
        #add graphs to excel file if score analysis was performed
        print("adding graphs...")
        graphsheet = workbook.add_worksheet('graphs')
        dumprow = 1              #row in which graph gets placed
        count = 1                                           
        for row in range(1,201):                            # goes through protein_index                                
            startcol = datastartcol
            samplecol = startcol
            chart = workbook.add_chart({'type': 'line'})
            chart.set_title({'name': ['data', row , 0, row, 0]}) 
            chart.set_size({'x_scale' : 1.5, 'y_scale': 1.5})
            chart.set_x_axis({'name': 'Molecular Mass (logarithmic)', 'label_position' : 'none'})
            chart.set_y_axis({'name': 'Relative Abundance', 'label_position' : 'none'})
            for i in range(len(samplelengths)):
                endcol = samplecol + samplelengths[i] -1 
                chart.add_series({'values': ['data', row, samplecol, row, endcol],
                                  'name'  : sample_names[i]})
                samplecol += samplelengths[i]
            graphsheet.write("A" + str(dumprow), count)
            graphsheet.insert_chart('B' + str(dumprow), chart)
            dumprow += 23
            count += 1

    writer.save()
    workbook.close()

def apply_cond_format(dataframe,startcol,endcol,writer,worksheet,workbook):
    """apply conditional formatting to relevant columns"""
    rownum = len(dataframe.index.values)
    worksheet.conditional_format(1,startcol,rownum, endcol,
                                              {'type'      : '3_color_scale',
                                               'min_color' : "#000000",
                                               'mid_type'  :  'percentile',
                                               'mid_value' : 95,
                                               'mid_color' : "#FFFF00",
                                               'max_color' : "#FF0000"})
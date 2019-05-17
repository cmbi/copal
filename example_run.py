import copal
import os

# test data paths
file_1 = os.path.join('test_data','Gel1_unaligned_samples_500_rows.csv')
file_2 = os.path.join('test_data','Gel2_unaligned_samples_500_rows.csv')
file_3 = os.path.join('test_data','Gel3_unaligned_samples_500_rows.csv')

# input dict values for analysis run with test files
input = {}
input['analysis_name'] = 'test_run'
input['filename'] = [file_1, file_2, file_3]
input['skiprows'] = [0, 0, 0]
input['input_type'] = [(',', '.'), (',', '.'), (',', '.')]
input['samplecolumns'] = [[(u'S1_G1_01', u'S1_G1_60'), (u'S2_G1_01', u'S2_G1_60'), (u'S3_G1_01', u'S3_G1_60'), (u'S4_G1_01', u'S4_G1_60'), (u'S5_G1_01', u'S5_G1_60')], [(u'S1_G2_01', u'S1_G2_60'), (u'S2_G2_01', u'S2_G2_60'), (u'S3_G2_01', u'S3_G2_60')], [(u'S1_G3_01', u'S1_G3_60')]]
input['identifier'] = ['Gi', 'Gi', 'Gi']
input['sheetname'] = ['', '', '']
input['samplenames'] = [u'S1_G1', u'S2_G1', u'S3_G1', u'S4_G1', u'S5_G1', u'S1_G2', u'S2_G2', u'S3_G2', u'S1_G3']
input['norm_check'] = True
input['normcol'] = None
input['normfile'] = None
input['align_check'] = True
input['hausdorff_scoring'] = False
input['gsea_output'] = False
input['GSEA_rank_column'] = 'Gene Symbol'
input['groups'] = [[1, 2, 3, 4, 5], [6, 7, 8, 9]]
input['hausd_factor'] = 1.0
input['warp_method'] = 'interpolate'

copal.main(input = input)

"""
----input dictionary: containing required input for COPAL analysis----

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
"""
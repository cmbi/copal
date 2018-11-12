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
input['hausdorff_scoring'] = 1
input['gsea_output'] = 1
input['GSEA_rank_column'] = 'Gene Symbol'
input['groups'] = [[1, 2, 3, 4, 5], [6, 7, 8, 9]]
input['hausd_factor'] = 1.0

copal.main(input = input)

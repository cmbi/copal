"""
script that renames lower scoring duplicates from a .rnk file

required command line arguments:
[1] -- path/filename of input .rnk file
[2] -- preferred path/filename of output .rnk file
"""

# import statements
import numpy as np
import pandas as pd
from sys import argv

# functions
def rename_duplicates(data):
    """
    renames duplicate gene symbols, highest scoring isoform keeps original name
    
    Args:
        data (pd_df): table with rank ordered list with index and following columns:
                                                - 'Symbol': Gene Symbol
                                                - 'Score' : Hausdorff effect size score
    Returns:
        data (pd_df): lower scoring isoforms have been renamed
    """
    symbols = list(set(data['Symbol']))
    for symbol in symbols:
        sub_rows = data.loc[data['Symbol'] == symbol]
        if len(sub_rows) > 1:
            sub_rows.sort_values('Score', ascending = False, inplace = True)
            for i in range(1,len(sub_rows)):
                idx =  sub_rows.iloc[i].name
                data.loc[idx,'Symbol'] += '_{}'.format(i+1)
    return data         
    

if __name__ == "__main__":
    
    # command line arguments:
    in_fname = argv[1]
    out_fname = argv[2]
    
    # load .rnk file
    df = pd.read_csv(in_fname, sep = '\t', header = None)
    df.columns = ['Symbol','Score']
    
    # rename duplicates
    result = rename_duplicates(df)
    
    # write updated .rnk file
    result.to_csv(out_fname, sep = '\t', header = False, index = False) 

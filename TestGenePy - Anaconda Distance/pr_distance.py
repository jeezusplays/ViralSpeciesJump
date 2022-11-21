"""PR_Distance.ipynb"""

# Start ===== Import Packages =====

import pandas as pd
pd.options.mode.chained_assignment = None #default= warn
#import numpy as np
import hashlib
from pylint import epylint as lint
(pylint_stdout, pylint_stderr) = lint.py_run('module_name.py', return_std=True)

# START ===== Read and Clean Data =====

original_df = pd.read_csv("Data/CSV_PR_original.csv") #import dataset

original_df.drop(['RefID', 'PtID', 'IsolateName', 'Region', 'Year', 'Subtype', 'PIList', 'AccessionID', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11','P12', 'P13', 'P14', 
        'P15', 'P16', 'P17', 'P18', 'P19', 'P20', 'P21', 'P22','P23', 'P24', 'P25', 'P26', 'P27', 'P28', 'P29', 'P30', 'P31', 'P32', 'P33','P34', 'P35', 'P36', 'P37', 'P38', 
        'P39', 'P40', 'P41', 'P42', 'P43', 'P44','P45', 'P46', 'P47', 'P48', 'P49', 'P50', 'P51', 'P52', 'P53', 'P54', 'P55','P56', 'P57', 'P58', 'P59', 'P60', 'P61', 'P62', 
        'P63', 'P64', 'P65', 'P66' ,'P67', 'P68', 'P69', 'P70', 'P71', 'P72', 'P73', 'P74', 'P75', 'P76', 'P77','P78', 'P79', 'P80', 'P81', 'P82', 'P83', 'P84', 'P85', 'P86', 
        'P87', 'P88','P89', 'P90', 'P91', 'P92', 'P93', 'P94', 'P95', 'P96', 'P97', 'P98', 'P99'], axis=1, inplace=True) #drop columns

original_df.columns = ['end'] #change 'NASeq' to 'end'
#original_df.head()

original_df.drop_duplicates(subset='end') #drop duplicate rows in 'end' column
#start = protease no mutation sequence
start = original_df[original_df['end'] == 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT' ].index

original_df.drop(start , inplace=True) #delete these row exact to 'start' from df
#original_df.shape #print shape of new df

#remove rows where sequences are not 297 string length
original_df['end'] = original_df['end'].astype('str')
mask = (original_df['end'].str.len() == 297) #df = df.drop(df[str.len(df.end) < 297].index)
original_df = original_df.loc[mask] #df.drop(df[mask].index)
#original_df.head()

original_df = pd.DataFrame(original_df) # Make edited dataframe as df to save it to a csv_file
original_df.to_csv('PR_Sequences.csv') 
df = pd.read_csv('PR_Sequences.csv') #import new dataset
df.shape #new shape of df with string length of 297

#check name of other column 
for col in df.columns: 
    print(col)
    
del df['Unnamed: 0'] #delete extra column

df_test = df.copy() #copy 'df' dataset
df_test = df_test.head(50) #extract 50 rows to use as testing
print(df_test)

#insert a column where all rows are protease start sequence
df_test['start'] = 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT'
df_test.insert(1, "hamming", "") #insert column 1 named hamming with no values
df_test.insert(2, "levenshtein", "") #insert column 2 named levenshtein with no values

df_test = df_test[['start', 'end', 'hamming', 'levenshtein']] #rearrange dataframe for better visualization
#df_test.head()

# Start ===== Hamming Distance =====
count = 0
row = df_test.shape[0]

for i in range(0, row): #for all i in all rows of df
    
    row_num = count
    start = df_test.iloc[row_num, 0]
    end = df_test.iloc[row_num, 1]
    hamming = df_test.iloc[row_num, 2]
    
    def hamming(start, end):
        return sum(start != end for start, end in zip(start, end))

    if __name__=="__main__":    
        start = hashlib.md5("start".encode()).hexdigest()
        end = hashlib.md5("end".encode()).hexdigest()

        start = df_test.iloc[row_num, 0]
        end = df_test.iloc[row_num, 1]

        assert len(start) == len(end)
        
        df_test.at[i, 'hamming'] = (hamming(start, end))
        
    count += 1

# Start ===== Iterative Levenshtein Distance =====
count = 0
row = df_test.shape[0]

for i in range(0, row):
    
    row_num = count
    start = df_test.iloc[row_num, 0]
    end = df_test.iloc[row_num, 1]
    levenshtein = df_test.iloc[row_num, 3]

    def levenshtein(start, end, costs=(1, 1, 1)):

        rows = len(start)+1
        cols = len(end)+1
        deletes, inserts, substitutes = costs
    
        dist = [[0 for x in range(cols)] for x in range(rows)]
        
        for col in range(1, cols):
            dist[0][col] = col * inserts
            for row in range(0, rows):
                if start[row-1] == end[col-1]:
                    cost = 0
                else:
                    cost = substitutes
                dist[row][col] = min(dist[row-1][col] + deletes,
                                     dist[row][col-1] + inserts,
                                     dist[row-1][col-1] + cost) # substitution
    
        return dist[row][col]

    df_test.at[i, 'levenshtein'] = (levenshtein(start, end))
    count += 1

#save new df into a csv file
df_test = pd.DataFrame(df_test)
df_test.to_csv('PR_Test.csv')
print("Calculated Hamming and Levenshtein Distance")
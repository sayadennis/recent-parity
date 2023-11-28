import os
import sys
import re
import numpy as np
import pandas as pd

sys.path.append('recent_parity/recent-parity/02_process_notes')
import ParseNotesGenetic

dn = '/share/fsmresfiles/breast_cancer_pregnancy/data/02_isolated_sections/genetic'

data = pd.read_csv(os.path.join(dn, 'isolated_sections_genetic_testing.csv'), sep='!')

# genelist = ['brca', 'brca1', 'brca2', 'palb2', 'palb', 'tp53', 'pten', 'cdh1', 'cdh', 'stk11', 'stk', 'chek2', 'atm']
genelist = ['brca', 'brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm']

#######################
#### Clean up text ####
#######################

for i in range(data.shape[0]):
    data.loc[i,'test_results'] = ParseNotesGenetic.clean_string(data.loc[i,'test_results'])

data.drop_duplicates(inplace=True, ignore_index=True)

#################################
#### Read and record results ####
#################################

result_table = pd.DataFrame(None, columns=['ir_id', 'general'] + genelist)
region_table = pd.DataFrame(None, columns=['ir_id', 'general'] + genelist)
for i in range(data.shape[0]):
    resultdict, regiondict = ParseNotesGenetic.string_to_dict(data.loc[i, 'ir_id'], data.loc[i,'test_results'])
    result_table = result_table.append(resultdict, ignore_index=True)
    region_table = region_table.append(regiondict, ignore_index=True)

# remove duplicates
result_table.drop_duplicates(inplace=True, ignore_index=True)
region_table.drop_duplicates(inplace=True, ignore_index=True)
# remove rows where everything except ir_id is NaN 
result_table = result_table.iloc[~np.all(pd.isnull(result_table.iloc[:,1:]), axis=1).values,:]
region_table = region_table.iloc[~np.all(pd.isnull(region_table.iloc[:,1:]), axis=1).values,:]

# result_table.to_csv('/share/fsmresfiles/breast_cancer_pregnancy/formatted_genetic_testing_results.csv', index=False)

##############################################################################################
#### Consolidate rows to create one-to-one relationship between test results and patients ####
##############################################################################################

dup = result_table.iloc[[sum(result_table['ir_id']==x)>1 for x in result_table['ir_id']],:]
dup_ids = dup['ir_id'].unique()

def consolidate(df_slice): # slices are single rows of Pandas DataFrane 
    consolidated = {}
    # get column names
    colnames = df_slice.columns
    # iterate through columns
    isnan = pd.isnull(df_slice)
    for colname in colnames:
        if np.all(isnan.loc[:,colname]): # all NaN
            continue
        elif np.sum(~isnan.loc[:,colname])==1: # only one is not NaN:
            consolidated[colname] = list(df_slice.loc[:,colname].values)[np.where(~isnan.loc[:,colname])[0][0]]
        else: # multiple non-NaN values
            if np.all([df_slice.loc[i,colname]==df_slice.loc[df_slice.index[0],colname] for i in df_slice.index]): # if values match for that specific columns
                consolidated[colname] = df_slice.loc[df_slice.index[0],colname] # simply append that value
            elif np.any([df_slice.loc[i,colname]=='positive' for i in df_slice.index]):
                consolidated[colname] = 'positive'
            elif np.any([df_slice.loc[i,colname]=='variant' for i in df_slice.index]):
                consolidated[colname] = 'variant'
            else: # if none of them are positive or variant... this should be caught by the first "if" so print if this case is met
                print('Uncertain which values to take for %s:' % colname)
                print(df_slice)
    return consolidated

singles = result_table.iloc[[sum(result_table['ir_id']==x)==1 for x in result_table['ir_id']],:]

for dup_id in dup_ids:
    singles = singles.append(consolidate(dup.iloc[dup['ir_id'].values==dup_id,:]), ignore_index=True)

singles.to_csv('/share/fsmresfiles/breast_cancer_pregnancy/data/formatted_genetic_testing_results.csv', index=False, header=True)


#### Do the same for regions ####

def consolidate_region(df_slice):
    consolidated = {}
    # get column names
    colnames = df_slice.columns
    # iterate through columns
    isnan = pd.isnull(df_slice)
    for colname in colnames:
        if np.all(isnan.loc[:,colname]): # all NaN
            continue
        elif np.sum(~isnan.loc[:,colname])==1: # only one is not NaN:
            consolidated[colname] = list(df_slice.loc[:,colname].values)[np.where(~isnan.loc[:,colname])[0][0]]
        elif len(df_slice[colname].values[0]) >= len(df_slice[colname].values[1]):
            consolidated[colname] = df_slice[colname].values[0]
        elif len(df_slice[colname].values[1]) > len(df_slice[colname].values[0]):
            consolidated[colname] = df_slice[colname].values[1]
        else:
            print('Uncertain which values to take for %s:' % colname)
            print(df_slice)
    return consolidated

dup = region_table.iloc[[sum(region_table['ir_id']==x)>1 for x in region_table['ir_id']],:]
dup_ids = dup['ir_id'].unique()

singles = region_table.iloc[[sum(region_table['ir_id']==x)==1 for x in region_table['ir_id']],:]

for dup_id in dup_ids:
    singles = singles.append(consolidate(dup.iloc[dup['ir_id'].values==dup_id,:]), ignore_index=True)

singles.to_csv('/share/fsmresfiles/breast_cancer_pregnancy/data/formatted_genetic_testing_results_regions.csv', index=False, header=True)




# for section in data['test_results']:
#     if ~np.any([
#         # ('Hypermethylation' in section), # discard this 
#         # ('Positive' in section), 
#         # ('POSITIVE:' in section),
#         ('mutation found-' in section),
#         ('mutation identified-' in section),
#         # ('BRCA1 (c.68_69del, c.5266dupC) and BRCA2 (c.5946del)' in section), # positive for BRCA1/2
#         # ('Likely Pathogenic' in section),
#         # ('Negative' in section), 
#         # ('Variant of uncertain significance' in section), 
#         # ('Variant of Uncertain Significance' in section),
#         # ('Variants of uncertain significance' in section), 
#         # ('Variant(s) of uncertain significance' in section),
#         # ('of uncertain significance-' in section),
#         # ('No mutation found' in section),
#         # ('No mutations found:' in section),
#         # ('No Mutation Detected' in section),
#         # ('No pathogenic mutation found' in section),
#         # ('No pathogenic mutations found ' in section),
#         # ('No BRCA mutation found' in section)
#         # ('Negative for the three common Ashkenazi Jewish' in section)
#     ]):
#         print(section + '\n')

#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import re
import sys
import time



#chromosome = 'S9'
#ncycles = 200
chromosome = sys.argv[1]
ncycles = sys.argv[2]

print("Starting on " + chromosome + " ncycles " + str(ncycles))



data = pd.read_csv('/media/drive_6tb/projects/kakapo-genomics/output/21_alphapeel/raw_alphapeel_files/'+ 
                   chromosome +                   
                   '_output_ncycles' +
                   str(ncycles) +
                   '.seg',
                         sep="\s+", header=None)

print("Done step 1/7: Read in data")


long_format = data.T
long_format.columns = long_format.iloc[0,]
del data
print("Done step 2/7: Transposed dataframe")


all_seq = pd.read_csv("~/projects/kakapo-genomics/data/all_seq_birdnames.txt", 
                       sep = "\t", header = None)[0].tolist()
unseq = ["Pegasus", "John-girl"]

all_the_names = all_seq + unseq
# all_the_names = ['Zephyr', 'Hoki']



blah = long_format.columns.unique()

all_colnames = []
for name in blah:
    for i in range(1,5):
        new_names = name + '_' + str(i)
        all_colnames.append(new_names)


long_format.columns = all_colnames


new_df = long_format.drop([0])
new_df['pos'] = range(1, len(new_df)+1)
del long_format
longer_df = pd.melt(new_df, id_vars = 'pos', value_vars= all_colnames, value_name = 'prob', var_name = 'anc_allele')
del new_df
print("Done step 3/7: Melted dataframe")

blah2 = longer_df.assign(pat_only = lambda dataframe: 
            dataframe['anc_allele'].map(lambda variable: 
                                      'pat' if re.search('_1$', variable) or re.search('_2$', variable) else 'mat'),
                         mat_only = lambda dataframe: 
            dataframe['anc_allele'].map(lambda variable: 
                                      'pat' if re.search('_1$', variable) or re.search('_3$', variable) else 'mat')
                        )
print("Done step 4/7: Recoded values")

blah4 = pd.melt(blah2,id_vars = ['anc_allele', 'prob', 'pos'], 
		value_vars = ['mat_only', 'pat_only'], 
		value_name = 'inherited', 
		var_name = 'focal_parent')
del blah2
blah4['chick'] = blah4['anc_allele'].str.replace('_\d$', '')
blah4.drop(['anc_allele'], axis = 1)
print("Done step 5/7: Melted again. Next step is very slow")

start_time = time.time()
blah4 = blah4.astype({'prob' : 'float'})
blah5 = blah4.groupby(['chick','pos', 'focal_parent', 'inherited'], as_index = False)[['prob']].sum()
del blah4
print("Slow step done, elapsed time: " + str(time.time() - start_time) + " seconds.")
print("Step 6/7: Added probabilities")
blah5.rename(columns= {'prob':'sum_val'}, inplace = True)



blah5[['pos', 'chick', 'focal_parent', 'inherited', 'sum_val']].to_csv('added_haplos_' + 
                          chromosome +
                          '_ncycles' +
                          str(ncycles) +
                          '.txt', 
             index = False,
             sep = '\t')
print("Step 7/7: Wrote datafile")
print("Done on " + chromosome + " ncycles " + str(ncycles))


# ** Glossary: pd ; df is the standar way of importing pandas package into a DataFrame **

import pandas as pd#To install pandas package, enter the command "pip install pandas" in the Terminal of editor
#need to import the panda package to work with this dataset

df = pd.read_csv('DRG_61122.csv')
print(df.to_string())# ** Using to.string() prints the entire DataFrame **
# prints the matrix of data ** Loads the CSV into a DataFrame **

av_row = df.mean(axis=1)
df['meanCT'] = av_row
print(df)
# computes the average of all technical replicates and bio replicates - the whole row

df['combined'] = df['gene_id'] + df['group']
#concatenation of columns gene id + group A/B

dCT = []
#defines list

# function to show all combined unique values i.e. TRPV1A/ TRPV1B
for n in set(df['combined']):
    dfU = df[df['combined'] == n]

# takes the for loop through each gene_id/ group combo

    goi = dfU[dfU['gene_type'] == 'GOI']['meanCT'].iloc[0]
    gapdh = dfU[dfU['gene_type'] == 'GAPDH']['meanCT'].iloc[0]
    dCT.append([n[0:-1], n[-1], goi - gapdh])
#computes dCT = CT_GOI-CT_GAPDH
dCT = pd.DataFrame(dCT, columns=['gene_id', 'group', 'dCT'])
#converts dCT from a list to a Pandas data frame


print(df)
# shows converted data set

#same steps below for the ddCT as for dCT
ddCT_fc = []

# set function to show all unique values i.e. all genes in a list duplicates removed
for n in set(dCT['gene_id']):
    ddfU = dCT[dCT['gene_id'] == n]

    # shows rows within group A
    A = ddfU[ddfU['group'] == 'A']['dCT'].iloc[0]
    B = ddfU[ddfU['group'] == 'B']['dCT'].iloc[0]
    # iloc= gives the value in the indexed location, first index in this case
    ddCT_fc.append([n, A - B, round(2 ** (-(A - B)), 2)])

ddCT_fc = pd.DataFrame(ddCT_fc, columns=['gene_id', 'ddCT', 'fc'])

print(ddCT_fc)

# bar chart script below
import numpy as np
import matplotlib.pyplot as plt

objects = ddCT_fc['gene_id']

y_pos = np.arange(len(objects))

performance = ddCT_fc['fc']
# extract the fc  column from tge ddCT_fc list

plt.bar(y_pos, performance, align='center', alpha=0.5)
plt.xticks(y_pos, objects)
plt.ylabel('fold change gene expression')
plt.title('Gene expression of neurotrophic factors in fractures compared to control')
plt.axhline(y=1,linewidth=1, color='k')
plt.show()
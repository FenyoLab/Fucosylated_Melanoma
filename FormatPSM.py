#! /usr/bin/env python

##Author Samson Jacob
##Format the file: LC_MS1.txt so that Gene_Name is included in a new column

import sys
import pandas as pd
import re

##Find the GN* in the description column from sequest output
def split_it(val):
    return re.findall('GN=([^\s]+)',val)

##read in the source file and create a new column with the Gene_Name Adjusted.
def Correct_Gene(datafr):
	vv = []
	df1 = pd.read_table(datafr)
	desc = df1.Description.tolist()
	for i in desc:
		m = split_it(i)
		if len(m)>0:
			vv.append(''.join(c for c in m if c not in '[]')) ## remove the brackets
		elif len(m)==0:
			vv.append('NaN')
	df1['Gene_Name_Corrected']= vv
	df1.to_csv('FORMATED_PSM.txt',sep='\t',na_rep='NaN')

if __name__ == '__main__':
	Correct_Gene(sys.argv[1])

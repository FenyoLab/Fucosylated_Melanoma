
#!/usr/bin/env python

"""
Author: Samson Jacob
NOTE: this was done verbosely for CLARITY !
"""
import pandas as pd


UNFILT=pd.read_csv("/Users/fenyolab/Desktop/CORRECT_PRAVEEN.txt",sep='\t')

##select the columns that have the word "PSM" in it
PSMCOLS = UNFILT.filter(regex='PSM')
##Add the Accession and Description information
PSMCOLS['Accession']=UNFILT.Accession
PSMCOLS['Description']=UNFILT.Description

##rename the columns
PSMCOLS.columns = ['PSM_4L1','PSM_4L2','PSM_Skmel1','PSM_Skmel2','PSM_Mewo1','PSM_Mewo2','Acession','Description']

##Drop the unused technical replicates
df = PSMCOLS.drop(['PSM_4L2','PSM_Skmel1','PSM_Mewo2'], axis=1)
##re-order the columns
df = df[['Accession','PSM_4L1','PSM_Skmel2','PSM_Mewo1','Description']]

##write_out the file for downstream processsing via R
df.to_csv('FORMATED_PSM.txt',sep='\t')

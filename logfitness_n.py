import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
pd.options.mode.chained_assignment = None  # default='warn'

#Handle the various arguments
parser = argparse.ArgumentParser(description='filter merged fastq reads to those of appropiate length')
parser.add_argument("-s1", "--string1", type=str, required=True)
parser.add_argument("-s2", "--string2", type=str, required=True)
parser.add_argument("-l", "--library", type=str, required=True)
print('start tabulation')

args = parser.parse_args()
#Define variables for either library a or b
if args.library is 'A':
    WT = 'MFKRLMMVALLVIAPLSAATAADQTNPYKLMDEAAQKTFDRLKNEQPQIRANPDYLRTIVDQELLPYVQVKYAGALVLGQYYKSATPAQREAYFAAFREYLKQAYGQ'
    codon = 0
else:
    WT = 'ALAMYHGQTYQIAPEQPLGDKTIVPIRVTIIDPNGRPPVRLDFQWRKNSQTGNWQAYDMIAEGVSMITTKQNEWGTLLRTKGIDGLTAQLKSISQQKITLEEKK'
    codon = 107

#Functions
#formater ensures the incoming dataframe is read as intergers and not strings. It also applies the RCPM normalization to all the values
def formater(table, length):
    table = table[0:21]
    for x in list(table.columns.values[0::]):
        table[x] = pd.to_numeric(table[x], errors='coerce').fillna(0).astype(np.int64)
        table[x] = (table[x]/length)*1000000
    return table

#ratioizer takes the ratio between RCPM of the selection/no-selection datasets and takes then takes the log10 value.
def ratioizer(table1, table2):
    output_table = table1 / table2
    output_table_log = output_table.apply(lambda x: np.log10(x) if np.issubdtype(x.dtype, np.number) else x)
    return output_table_log


#subtractor will subtract the mutant log10 ratio by the WT log10 ratio for each position on the protein.
def Substractor(table, codon):
    for x in WT:
        codon += 1
        row = table.loc[table["Amino Acid"] == x].index[0]
        WT_value = table.loc[row, "codon{}".format(codon)]
        table["codon{}".format(codon)] = table["codon{}".format(codon)].apply(lambda x: x - WT_value)
    return(table)

Length1 =0
Length2 =0


#read in the two data sets (table 1 = selection condition; table 2 = no-selection condition)
table1 = pd.read_csv("{}.raw_counts.csv".format(args.string1), keep_default_na=False)
table1.rename( columns={'Unnamed: 0':'Amino Acid'}, inplace=True )
table1_r = table1.drop('Amino Acid', 1)
table2 = pd.read_csv("{}.raw_counts.csv".format(args.string2), keep_default_na=False)
table2.rename( columns={'Unnamed: 0':'Amino Acid'}, inplace=True )
table2_r = table2.drop('Amino Acid', 1)

#Two loops to count the total sequences in each table for normalization (slowest part of the program)
with open("{}.translated.fasta".format(args.string1), "w") as handle:
    for seq_record in SeqIO.parse("{}.fasta".format(args.string1), "fasta"):
        Length1 = Length1+1

with open("{}.translated.fasta".format(args.string2), "w") as handle:
    for seq_record in SeqIO.parse("{}.fasta".format(args.string2), "fasta"):
        Length2 = Length2+1

#Call the formater function for table 1 and table 2
table1_num = formater(table1_r, Length1)
table2_num = formater(table2_r, Length2)



#Call the ratioizer function for the two tables
table_log = ratioizer(table1_num, table2_num)

#Insert a column at position 0 in the table
AminoAcid =('A', 'R', 'N', 'D', 'C', 'E', 'Q', \
          'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T','W', 'Y','V', "*")

table_log.insert(0, 'Amino Acid', AminoAcid)

#call the subtractor function on the log10 ratio table
table_log_dif = Substractor(table_log, codon)

#export final table
table_log_dif.to_csv('diflog_{}_{}.csv'.format(args.string1,args.string2), index=False)

import pandas as pd
from Bio import SeqIO 
import sys

tsv_file = sys.argv[1]
fasta_file = sys.argv[2]
output_file = sys.argv[3]
Maximum_Ns = sys.argv[4]

df = pd.read_csv(tsv_file, sep="\t")
df["MAF"] = df["AF"]
df = df[["ID", "coordinates", "MAF", "REF", "ALT"]]
# Add a new column with the part after the ':' from the ID column
df['position'] = df['ID'].str.split(':').str[1]
# Add a new column with the ratio of REF to ALT
df['Value'] = '[' + df['REF'] + '/' + df['ALT'] + ']'
df.set_index('coordinates', inplace=True)
df["Seq"] = ""

# Parse the FASTA file and fill the Seq column

for record in SeqIO.parse(fasta_file, "fasta"):
    header = record.id
    if header in df.index:
        sequence = str(record.seq)
        Value = df.loc[header, 'Value']
        sequence_list = list(sequence)
        sequence_list[30] = Value
        sequence = ''.join(sequence_list)
        df.loc[header, 'Seq'] = sequence

# Now add N filter
df["Ns"] = df['Seq'].str.count('N') 

df_filtered = df[df['Ns'] <= Maximum_Ns]

df_filtered = df_filtered.sort_values("MAF", ascending=False).head(n=1000)
# Export data
df_filtered.to_csv(path_or_buf= output_file)



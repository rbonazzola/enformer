import pandas as pd
from pyfaidx import Fasta
from tqdm import tqdm

# Load the FASTA file

regions_df = pd.read_csv("sequences.bed", sep="\t", header=None)
regions_df.columns = ['chr', 'start', 'end', 'partition']

regions_df.start -= 131_072
regions_df.end += 131_072
# regions_df = df.query("chr == @chromosome").head(20)

# Loop through regions and write each to a separate file
for _, row in tqdm(regions_df.iterrows()):
  chrom, start, end = row.chr, row.start, row.end
  fasta = Fasta(f'data/{chrom}.fa')
  region_seq = fasta[chrom][start:end]  # 0-based indexing in pyfaidx
  output_file = f'regions/{chrom}_{start+1}_{end}.fa'
                
  with open(output_file, 'w') as f:
    f.write(f'>{chrom}:{start+1}-{end}\n')
    f.write(str(region_seq) + '\n')

print("Regions extracted and saved.")

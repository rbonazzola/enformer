import argparse
import pandas as pd
import numpy as np
from pyfaidx import Fasta
from tqdm import tqdm
from IPython import embed

ONE_HOT_ENCODING = {
    'A': [1, 0, 0, 0],
    'C': [0, 1, 0, 0],
    'G': [0, 0, 1, 0],
    'T': [0, 0, 0, 1],
    'a': [1, 0, 0, 0],
    'c': [0, 1, 0, 0],
    'g': [0, 0, 1, 0],
    't': [0, 0, 0, 1],
    'N': [0, 0, 0, 0],
    'n': [0, 0, 0, 0]
}


def one_hot_encode_sequence(sequence):    
    
    one_hot_encoded_sequence = []
    
    for nucleotide in sequence:
        if nucleotide in ONE_HOT_ENCODING:
            one_hot_encoded_sequence.append(ONE_HOT_ENCODING[nucleotide])
        else:
            raise ValueError(f"Unexpected nucleotide: {nucleotide}")
    
    return np.array(one_hot_encoded_sequence)


def process_chunk(regions_df, chunk_index, num_chunks):

    chunk_size = len(regions_df) // num_chunks
    start_index = (chunk_index - 1) * chunk_size  # Adjust for 1-based indexing    
    end_index = len(regions_df) if chunk_index == num_chunks else start_index + chunk_size

    chunk_df = regions_df.iloc[start_index:end_index]

    for _, row in tqdm(chunk_df.iterrows(), total=len(chunk_df)):
        chrom, start, end = row.chr, row.start, row.end
        fasta = Fasta(f'data/{chrom}.fa')
        region_seq = fasta[chrom][start:end]
        encoded_seq = one_hot_encode_sequence(str(region_seq))
        output_file = f'regions_mouse_npy/{chrom}_{start+1}_{end}.npy'
        np.save(output_file, encoded_seq)

    print(f"Chunk {chunk_index}/{num_chunks} processed and saved.")


def main():

    parser = argparse.ArgumentParser(description="Split regions and one-hot encode sequences.")
    
    parser.add_argument('--species', type=str, required=True, default="human", help="Species to process (human or mouse).")
    parser.add_argument('--reference_genome_file_format', type=str, required=True, default="data/datasets/reference_{species}_genome/{chrom}.fa")
    parser.add_argument('--output_hdf5', type=str, required=True, default="data/transforms/{species}_genome/regions_{species}_{chunk_index}.hdf5")    
    parser.add_argument('--regions-file', type=str, required=True, help="Path to the regions file (BED format).")
    parser.add_argument('--num-chunks', type=int, required=True, help="Total number of chunks to split the regions into.")
    parser.add_argument('--chunk-index', type=int, required=True, help="Index of the chunk to process (1-based).")
    
    args = parser.parse_args()

    # Ensure that the chunk index is within the valid range
    assert 1 <= args.chunk_index <= args.num_chunks, '''Chunk index must be between 1 and the total number of chunks.'''

    regions_df = pd.read_csv(args.regions_file, sep="\t", header=None)
    # regions_df.columns = ['chr', 'start', 'end', 'partition']
    regions_df.columns = ['chr', 'start', 'end']
    print(regions_df.end-regions_df.start)
    
    #regions_df.start -= 131_072
    #regions_df.end += 131_072
    
    process_chunk(regions_df, args.chunk_index, args.num_chunks)

if __name__ == "__main__":
    main()

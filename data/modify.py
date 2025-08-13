import os
import pandas as pd
import pickle as pkl
import urllib.request
from tqdm import tqdm

import cyvcf2

from personalized import (
    find_variants_in_vcf_file,
    extract_reference_sequence,
    replace_variants_in_reference_sequence,
    create_mapping_dictionary,
    get_fastaExtractor,
    compute_avg_sequence,
)

ENFORMER_DATA_DIR = "/home/home01/scrb/nobackup/enformer/"
# ENFORMER_DATA_DIR = "/mnt/data/data/enformer/"
GEUVADIS_PATH = f"{ENFORMER_DATA_DIR}/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt"
FASTA_PATH    = f"{ENFORMER_DATA_DIR}/genomes/reference_human_genome"
VCF_PATH      = f"{ENFORMER_DATA_DIR}/genomes/1000G/vcfs/1000G_GRCh38_vcfs/"
VCF_FILE_PATTERN = f"{VCF_PATH}/ALL.chr{{chromosome}}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

for path in [GEUVADIS_PATH, VCF_PATH]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Path {path} does not exist. Please check the path.")


def download_geuvadis_sample_list(output_path: str = "datasets/E-GEUV-1.sdrf.txt") -> str:
    """
    Downloads the GEUVADIS SDRF file from ArrayExpress and returns its path.
    Skips download if the file already exists.
    """
    url = "https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt"

    if os.path.exists(output_path):
        print(f"📎 GEUVADIS SDRF already exists at {output_path}, skipping download.")
    else:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        print(f"⬇️ Downloading GEUVADIS SDRF from {url}")
        urllib.request.urlretrieve(url, output_path)
        print(f"💾 Saved to {output_path}")

    return output_path


def extract_geuvadis_european_ids(sdrf_path: str) -> list[str]:
    """
    Extracts European GEUVADIS sample IDs (CEU, FIN, GBR, TSI) from the SDRF file.
    """
    print(f"📋 Extracting European GEUVADIS sample IDs from {sdrf_path}")
    european_pops = {"Utah", "Finnish", "British", "Tuscan"}
    sample_ids = set()

    with open(sdrf_path) as f:
        header = f.readline().strip().split("\t")
        print(header)
        pop_idx = header.index("Characteristics[ancestry category]")
        sample_idx = header.index("Source Name")

        for line in f:
            fields = line.strip().split("\t")
            population = fields[pop_idx].replace("population: ", "").strip()
            sample = fields[sample_idx].strip()
            if population in european_pops:
                sample_ids.add(sample)

    print(f"✅ Found {len(sample_ids)} European GEUVADIS sample IDs")
    return sorted(sample_ids)


def pipeline(interval, vcf_file, samples, fasta_extractor):
    
    # 1. Extract reference sequence from FASTA
    ref = extract_reference_sequence(region=interval, fasta_func=fasta_extractor, resize_for_enformer=True)    
    onehot_ref = ref["sequence"]           # one-hot encoded reference sequence
    interval = ref["interval"]             # dict with chrom, start, end

    # 2. Load VCF and extract variants per sample within the interval
    vcf = cyvcf2.VCF(vcf_file, samples=samples)
    variant_calls = find_variants_in_vcf_file(cyvcf2_object=vcf, interval=interval, samples=samples, chromosome_length=ref['max_position'])
    
    # 3. Map variants to relative positions and encode haplotypes
    mapping = create_mapping_dictionary(variants_array=variant_calls, interval_start=interval["start"], haplotype="both", samples=samples)

    # 4. Replace variants in the reference one-hot sequence
    encoded_sequences = replace_variants_in_reference_sequence(query_sequences_encoded=onehot_ref, mapping_dict=mapping, samples=samples)
    avg_sequence = compute_avg_sequence(encoded_sequences)

    return encoded_sequences, avg_sequence


def process_chunk(regions_df, chunk_index, num_chunks):

    chunk_size = len(regions_df) // num_chunks
    start_index = (chunk_index - 1) * chunk_size  # Adjust for 1-based indexing    
    end_index = len(regions_df) if chunk_index == num_chunks else start_index + chunk_size

    chunk_df = regions_df.iloc[start_index:end_index]
    return chunk_df


def get_interval_from_region(region):
    """
    Get the interval from the region.
    """
    chromosome, start, end = region.chr.replace("chr", ""), region.start, region.end
    if chromosome != 'X':
        chromosome = int(chromosome)
    start, end = int(start), int(end)
    interval = f"chr{chromosome}_{start}_{end}"
    return chromosome, start, end, interval


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="Geuvadis reference generator")
    
    parser.add_argument('--regions-file', "--regions_file", type=str, required=False, help="Path to the regions file (BED format)", default="sequences_human.bed")
    parser.add_argument('--num-chunks', "--num_chunks", type=int, required=True, help="Total number of chunks to split the regions into.")
    parser.add_argument('--chunk-index', "--chunk_index", type=int, required=True, help="Index of the chunk to process (1-based).")
    parser.add_argument('--output-folder', "--output_folder",  type=str, required=False, help="Output folder", default="EUR_reference")

    args = parser.parse_args()
    args.chunk_index = int(args.chunk_index)
    args.num_chunks  = int(args.num_chunks)

    regions_df         = pd.read_csv(args.regions_file, sep="\t", header=None)
    regions_df.columns = ['chr', 'start', 'end', 'partition']
    regions_df.start  -= 131_072
    regions_df.end    += 131_072

    # geuvadis_lcl = pd.read_csv(GEUVADIS_PATH, sep="\t")
    # geuvadis_samples = geuvadis_lcl.columns[4:].to_list()
    
    geuvadis_file = download_geuvadis_sample_list()
    eur_geuvadis_samples = extract_geuvadis_european_ids(geuvadis_file)

    region_subdf = process_chunk(regions_df, args.chunk_index, args.num_chunks)

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    for _, region in tqdm(region_subdf.iterrows(), total=len(region_subdf)):
                
        chromosome, start, end, interval = get_interval_from_region(region)
        output_file = f"{args.output_folder}/chr{chromosome}_{start}-{end}.npy"
        
        if os.path.exists(output_file):
            print(f"File {output_file} already exists, skipping.")
            continue

        fasta_file = os.path.join(FASTA_PATH, f"chr{chromosome}.fa")
        fasta_extractor = get_fastaExtractor(fasta_file)
        
        vcf_file = VCF_FILE_PATTERN.format(chromosome=chromosome)
        
        print(vcf_file)
        encoded_sequences, avg_sequence = pipeline(interval, vcf_file, eur_geuvadis_samples, fasta_extractor)
        
        pkl.dump(avg_sequence, open(output_file, "wb"))

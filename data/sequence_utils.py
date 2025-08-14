import os, sys
import subprocess
from pathlib import Path
from tqdm import tqdm

import numpy as np
import tensorflow as tf

from Bio import SeqIO

repo_root = Path("..")
os.chdir(repo_root)

HUMAN_REF_FOLDER = os.getenv("HUMAN_REF_FOLDER", "./data/datasets/ref/human")
HUMAN_TFR_FOLDER = os.getenv("HUMAN_TFR_FOLDER", "./data/datasets/basenji/human")

MOUSE_REF_FOLDER = os.getenv("MOUSE_REF_FOLDER", "./data/datasets/ref/mouse")
MOUSE_TFR_FOLDER = os.getenv("MOUSE_TFR_FOLDER", "./data/datasets/basenji/mouse")

HUMAN_REF_FOLDER = Path(HUMAN_REF_FOLDER)
HUMAN_TFR_FOLDER = Path(HUMAN_TFR_FOLDER)

MOUSE_REF_FOLDER = Path(MOUSE_REF_FOLDER)
MOUSE_TFR_FOLDER = Path(MOUSE_TFR_FOLDER)


def deserialize(serialized_example, metadata):

    feature_map = {
        'sequence': tf.io.FixedLenFeature([], tf.string),
        'target': tf.io.FixedLenFeature([], tf.string),
    }
    example = tf.io.parse_example(serialized_example, feature_map)
    sequence = tf.io.decode_raw(example['sequence'], tf.bool)
    sequence = tf.reshape(sequence, (metadata['seq_length'], 4))
    sequence = tf.cast(sequence, tf.float32)

    target = tf.io.decode_raw(example['target'], tf.float16)
    target = tf.reshape(target, (metadata['target_length'], metadata['num_targets']))
    target = tf.cast(target, tf.float32)

    return {'sequence': sequence, 'target': target}


def load_tfrecord_to_numpy(tfrecord_path, metadata):
    dataset = tf.data.TFRecordDataset([tfrecord_path], compression_type='ZLIB')
    dataset = dataset.map(lambda x: deserialize(x, metadata))
    sequences = []
    targets = []
    for example in dataset:
        sequences.append(example['sequence'].numpy())
        targets.append(example['target'].numpy())
    sequences = np.stack(sequences)
    targets = np.stack(targets)
    return {'sequence': sequences, 'target': targets}


def onehot_to_seq(arr):
    # arr: (length, 4), float32 values
    indices = np.argmax(arr, axis=1)
    return ''.join(np.array(['a', 'c', 'g', 't'])[indices])


def extract_region(region, expected_length=131_072):
    if not str(region.chromosome).startswith("chr"):
        chromosome = f"chr{region.chromosome}"
    else:
        chromosome = region.chromosome
    start = int(region.start)
    end   = int(region.end)
    sequence = mm10_per_chr[chromosome][start:end]

    assert (len(sequence)) == expected_length, f"Difference between start and end is not {expected_length} but {len(sequence)}."
    return sequence


def get_all_mouse_sequences(path=MOUSE_TFR_FOLDER):
    
    metadata = {
      'seq_length': 131072,
      'target_length': 896,
      'num_targets': 1643
    }

    all_seqs = []
    for basename in tqdm(os.listdir(path)):
        tfrecord_path = path + basename
        record_data = load_tfrecord_to_numpy(tfrecord_path, metadata)
        seqs = record_data['sequence'].astype(np.int8)
        all_seqs.append(seqs)
    
    all_seqs = np.concatenate(all_seqs, axis=0)
    return all_seqs


def get_hg38():
    genome_fasta = HUMAN_REF_FOLDER / "hg38.fa"
    seqs_per_chr = { record.id: str(record.seq).lower() for record in tqdm(SeqIO.parse(genome_fasta, "fasta")) }
    return seqs_per_chr


def get_mm10():
    genome_fasta = MOUSE_REF_FOLDER / "mm10.fa"
    seqs_per_chr = { record.id: str(record.seq).lower() for record in tqdm(SeqIO.parse(genome_fasta, "fasta")) }
    return seqs_per_chr
    
get_data_subset = lambda tfr_filename: tfr_filename.replace(".tfr", "").split('-')[0]   
get_record_number = lambda tfr_filename: int(tfr_filename.replace(".tfr", "").split('-')[2])


def get_human_record_ids():
        
    return sorted(
        [ f.replace(".tfr", "") for f in os.listdir(HUMAN_TFR_FOLDER / "tfrecords") ], 
        key=lambda f: (get_data_subset(f), get_record_number(f))
    )


def get_mouse_record_ids():
    
    return sorted(
        [ f.replace(".tfr", "") for f in os.listdir(MOUSE_TFR_FOLDER  / "tfrecords") ], 
        key=lambda f: (get_data_subset(f), get_record_number(f))
    )


def get_sequences_for_record(record_id):
    seq_for_record = { 
        record.id: str(record.seq).lower() 
        for record in tqdm(SeqIO.parse(f"fasta_seq/{record_id}.fasta", "fasta"))
    }    
    return seq_for_record


def get_human_basenji_regions():
    human_seqs = pd.read_csv(HUMAN_TFR_FOLDER / "sequences.bed", sep='\t', header=None)
    human_seqs = human_seqs.set_axis(["chromosome", "start", "end", "subset"], axis=1)
    return human_seqs


def get_mouse_basenji_regions():
    mouse_seqs = pd.read_csv( MOUSE_TFR_FOLDER / "sequences.bed", sep='\t', header=None)
    mouse_seqs = mouse_seqs.set_axis(["chromosome", "start", "end", "subset"], axis=1)
    return mouse_seqs


def expand_regions(regions_df, to_left:int=131_072, to_right:int=131_072):
    regions_df.start -= to_left
    regions_df.end += to_right
    return regions_df


def get_region_from_npy_filename(filename):
    return {
      'chr': filename.split("_")[0], 
      'start': int(filename.split("_")[1]), 
      'end': int(filename.split("_")[2].split(".")[0])
    }


def write_fasta(sequences, output_path, cut_every=80):
    with open(output_path, "w") as f:
        for name, seq in sequences:
            f.write(f">{name}\n")
            # cut lines every cut_every characters
            for i in range(0, len(seq), cut_every):
                f.write(seq[i:i+cut_every] + "\n")


def get_hash(seq):
    return hashlib.sha256(seq.encode()).hexdigest()


def find_matches(genome_fasta, tfrecord_hashes, window_size=131072):
    for record in SeqIO.parse(genome_fasta, "fasta"):
        chrom = record.id
        seq = str(record.seq).upper()
        matches = []
        for i in range(len(seq) - window_size + 1):
            subseq = seq[i:i+window_size]
            h = get_hash(subseq)
            if h in tfrecord_hashes:
                matches.append((chrom, i, i+window_size, h))
                print(f"Match at {chrom}:{i}-{i+window_size}")
        return matches


def get_tfr_files(tfr_folder):

    train_files = [ tfr_folder / f"train-1-{i}.tfr" for i in range(200) ]
    valid_files = [ tfr_folder / f"valid-1-{i}.tfr" for i in range(200) ]
    test_files  = [ tfr_folder / f"test-1-{i}.tfr" for i in range(200) ]
    train_files = [ f for f in train_files if os.path.exists(f) ]
    valid_files = [ f for f in valid_files if os.path.exists(f) ]
    test_files  = [ f for f in test_files  if os.path.exists(f) ]

    return train_files, valid_files, test_files
    

def chunk_by_subset(regions_df, chunk_size=256):
    result = {}

    for subset in ["train", "valid", "test"]:
        df_subset = regions_df[regions_df["subset"] == subset]
        chunks = [df_subset.iloc[i:i+chunk_size] for i in range(0, len(df_subset), chunk_size)]
        for i, chunk in enumerate(chunks):
            name = f"{subset}-1-{i}"
            result[name] = chunk

    return result


def get_hash_from_seq(seq):
    return hashlib.sha256(seq.encode()).hexdigest()


def find_matches(genome_fasta, tfrecord_hashes, window_size=131072):
    for record in SeqIO.parse(genome_fasta, "fasta"):
        chrom = record.id
        seq = str(record.seq).upper()
        matches = []
        for i in range(len(seq) - window_size + 1):
            subseq = seq[i:i+window_size]
            h = get_hash_from_seq(subseq)
            if h in tfrecord_hashes:
                matches.append((chrom, i, i+window_size, h))
                print(f"Match at {chrom}:{i}-{i+window_size}")
        return matches


def get_all_hashes():
    hashes = [ open("sequence_hashes/"+x, "rt").readlines() for x in os.listdir("sequence_hashes") ]
    return hashes


def flatten_list(lst):
    return[x for y in lst for x in y]
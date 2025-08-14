import os
import numpy as np
import pandas as pd
import tensorflow
import json

from Bio import SeqIO
from tqdm import tqdm
from pathlib import Path

from tfr_to_hash import *

from sequence_utils import (
    deserialize,
    load_tfrecord_to_numpy,
    onehot_to_seq,
    write_fasta
)

metadata = dict(seq_length=131072, target_length=896, num_targets=1643)

filename = sys.argv[1]

tfrecord_path = "tfrecords" + filename
record_data = load_tfrecord_to_numpy(tfrecord_path, metadata)
seqs = record_data['sequence'].astype(np.int8)
sequences = [ onehot_to_seq(seq) for seq in tqdm(seqs) ]
write_fasta([(f"seq{i}", seq) for i, seq in enumerate(sequences)], f"{filename.replace('.tfr','')}.fasta")
import os
import json
from pathlib import Path
from tqdm import tqdm
import gc

import numpy as np
import h5py
import pandas as pd
import tensorflow as tf

from Bio import SeqIO

from sequence_utils import (
    deserialize,
    load_tfrecord_to_numpy,
    get_tfr_files,
    MOUSE_TFR_FOLDER, MOUSE_REF_FOLDER,
    HUMAN_TFR_FOLDER, HUMAN_REF_FOLDER
)

# Change this by 'human' if you want to add targets to the human dataset
SPECIES = "mouse"

assert SPECIES in ["mouse", "human"], f"SPECIES should be 'mouse' or 'human' but {SPECIES=}"

TFR_FOLDER = MOUSE_TFR_FOLDER if SPECIES == "mouse" else HUMAN_TFR_FOLDER

train_files, valid_files, test_files = get_tfr_files(tfr_folder=TFR_FOLDER)
train_files = train_files + valid_files

NUM_TARGETS = 1643 if SPECIES == 'mouse' else 5803

metadata = dict(seq_length=131072, target_length=896, num_targets=NUM_TARGETS )

h5_filename, tfr_files = f"test_{SPECIES}_copy.h5" , test_files
h5_filename, tfr_files = f"train_{SPECIES}_copy.h5", train_files

with h5py.File(h5_filename, "r+") as f:
    dset_target = f["target"]
    num_samples = dset_target.shape[0]

    start = 0
    for i, tfr_file in enumerate(tfr_files):
        with tf.device('/CPU:0'):
            tgt_chunk = load_tfrecord_to_numpy(tfr_file, metadata=metadata)['target']
        end = start + len(tgt_chunk)
        print(i, start, end)
        dset_target[start:end] = tgt_chunk
        start = end
        del tgt_chunk
        gc.collect()
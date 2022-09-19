#!/usr/bin/env python3
from joblib import Parallel, delayed
import numpy as np
import glob
import argparse
from shutil import copyfileobj
import bz2
import gzip
import os
from pathlib import Path

def compress_function(file, delete, compress):
    if file.endswith('.bz2') or file.endswith('.gz'):
        return
    print (f"Compress file: {file} is started.")
    with open(file, 'rb') as input:
        if compress=='bz2':
            with bz2.BZ2File(f'{file}.bz2', 'wb', compresslevel=9) as output:
                copyfileobj(input, output)
        else:
            with gzip.open(f'{file}.gz', 'wb') as output:
                copyfileobj(input, output)
    print (f"Compress file: {file} is finished.")
    os.remove(file)
    if delete:
        print (f"Uncompressed file: {file} is deleted.")

parser = argparse.ArgumentParser()
parser.add_argument("-pattern", type=str, help="Pattern for compress file: pattern=dump-* is by default.", nargs='?', default="dump-*")
parser.add_argument("-recursive", type=int, help="Recursively file compression: recursive=0 is by default.", nargs='?', default=0)
parser.add_argument("-delete", type=int, help="Delete uncompressed file 1, otherwise 0. delete=0 is by default.", nargs='?', default=0)
parser.add_argument("-compress", type=str, help="Specify compress format: [gz, bz2]. compress=.gz is by default.", nargs='?', default="gz")
parser.add_argument("-jobs", type=int, help="Specify number of jobs. jobs=-1 is by default: all cores will be run.", nargs='?', default=-1)
args = parser.parse_args()
print(args)
pattern = args.pattern
delete = args.delete
compress = args.compress
n_jobs = args.jobs
if args.recursive:
    recursive = True
else:
    recursive = False
if recursive:
    files = [str(p) for p in Path('.').rglob(f'{pattern}')]
else:
    files = glob.glob(f'{pattern}', recursive=False)
print(f'files: {files}')
element_run = Parallel(n_jobs=n_jobs)(delayed(compress_function)(file, delete, compress) for file in files)

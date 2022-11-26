#!/usr/bin/env python3
import gzip
import os
import sys

input_file = sys.argv[1]
output_file = "tmp_file.gz" #tempfile.NamedTemporaryFile(mode="w", delete=False)
with gzip.open(input_file, "rb") as in_file, gzip.open(output_file, "wb") as out_file:
    for i, l in enumerate(in_file):
        line = l.decode('utf-8')
        if i % 4 == 0:
            out_file.write(f"{'_'.join(line.split('_')[:-1])}\n".encode('utf-8'))
        else:
            out_file.write(l)

# os.rename(output_file, input_file)

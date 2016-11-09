#!/bin/bash

tar -xzf python3.tar.gz
export PATH=miniconda3/bin:$PATH
python calc_L57_chunk_vi.py $1

#!/bin/bash

tar -xzf python3.tar.gz
export PATH=miniconda3/bin:$PATH
python chunk_L57_bands.py $1 $2

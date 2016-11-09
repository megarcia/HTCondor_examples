#!/bin/bash

tar -xzf python3.tar.gz
tar -xzf $1_$2.tar.gz
export PATH=miniconda3/bin:$PATH
python stitch_vi.py $1 $2

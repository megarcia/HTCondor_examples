#!/bin/bash

tar -xzf python3.tar.gz
export PATH=miniconda3/bin:$PATH
python diff_vi.py $1 $2 $3

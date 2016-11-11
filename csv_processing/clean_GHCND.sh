#!/bin/bash

tar -xzf python3.tar.gz
export PATH=miniconda3/bin:$PATH
python clean_GHCND.py $1

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Tests for acrg_tools
"""

import pytest
import gzip
import bz2
from acrg_tools import file_decompression as decompress


from acrg_config.paths import paths
acrg_path = paths.acrg

def test_decompress_files(tmp_path):
    string_one = b"this is a string of words"
    string_two = b"this is another string of words"
    string_three = b"the earth is flat and sits on a flying turtle"

    strings_to_compress = [string_one, string_two, string_three]

    # Test compressed and decompression by the function with bz2
    bzip_paths = []
    for i, s in enumerate(strings_to_compress):
        tmp_filepath_bz2 = tmp_path / f"{i}.bz2"
        
        with bz2.open(tmp_filepath_bz2, "wb") as f:
            f.write(s)
        
        bzip_paths.append(tmp_filepath_bz2)

    decompressed_bz = decompress.decompress_files(bzip_paths)

    for fpath, test_str in zip(decompressed_bz, strings_to_compress):
        with open(fpath, "r") as f:
            read_str = f.read()

        assert test_str.decode("utf-8") == read_str

    # Now do the same but with gzip
    gzip_paths = []
    for i, s in enumerate(strings_to_compress):
        tmp_filepath_gz = tmp_path / f"{i}.gz"

        with gzip.open(tmp_filepath_gz, "wb") as f:
                f.write(s)
        
        gzip_paths.append(tmp_filepath_gz)

    decompressed_gz = decompress.decompress_files(gzip_paths)
    
    for fpath, test_str in zip(decompressed_gz, strings_to_compress):
        with open(fpath, "r") as f:
            read_str = f.read()

        assert test_str.decode("utf-8") == read_str
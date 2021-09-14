#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Tests for acrg_tools
"""
import gzip
import bz2
import os 
import shutil

from acrg_config.paths import paths
acrg_path = paths.acrg

def decompress_files(filepaths):
    """ Decompress files at paths in filepaths

        Args:
            filepaths (list): List of paths to files to decompress
        Returns:
            list: List of paths to decompressed files
    """
    _, extension = os.path.splitext(filepaths[0])

    if extension == ".bz2":
        complib = bz2
    elif extension == ".gz":
        complib = gzip
    else:
        raise ValueError("Unable to decompress files other than bz2 and gz")
        
    decompressed_files = []
   
    for compressed_path in filepaths:
        # Here decompressed_path will be the path minus the extension
        decompressed_path, extension = os.path.splitext(compressed_path)

        with complib.open(compressed_path, "rb") as f_in, open(decompressed_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        decompressed_files.append(decompressed_path)

    return decompressed_files
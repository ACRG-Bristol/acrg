# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 10:28:38 2015

@author: chxmr
"""

from cPickle import dump as cpdump
from cPickle import load as cpload
import gzip

def save(filename, list_of_objects):
    '''
    Save and compress Python objects.
    Example:
    
    save("/my/file/path/file.p.gz", (obja, objb, objc))
    '''
    
    with gzip.open(filename, "wb") as f:
        cpdump(list_of_objects, f, protocol = 2)

def load(filename):
    '''
    Load a list of Python objects stored in a compressed file
    Example:
    
    obja, objb, objc = load("/my/file/path/file.p.gz")
    '''

    with gzip.open(filename, "rb") as f:
        list_of_objects = cpload(f)
    return list_of_objects

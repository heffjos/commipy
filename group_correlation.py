#!/usr/bin/env python3

import os
import argparse
import sys
import time

import nibabel as ni
import numpy as np

import cifti_helpers as CH

class groupData():
    def __init__(self):
        self.in_files = {}

    def parse_data(input_csv):
        with open(input_csv, "r") as f:
            lines = f.readlines()
        
    
def init_parser(argv=sys.argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--roi", required=True,
        help="path to mask as dscalar.nii")
    parser.add_argument("--prefix", default="./output",
        help="output file prefix")
    parser.add_argument("--individual", action="store_true",
        help="save connectivity maps for each individual")
    parser.add_argument("--input", nargs="+",
        help="input csv file c1:file path, c2:file identifier")
    return parser

def group_correlation(scans, roi, 

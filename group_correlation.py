#!/usr/bin/env python3

import os
import argparse
import sys
import time

import nibabel as ni
import numpy as np
import pandas as pd

import cifti_helpers as CH

def parse_data(input_csv):
    with open(input_csv, "r") as f:
        lines = f.readlines()

    results = []
    for i, line in enumerate(lines):
        tokenized = line.split(",")
        if len(tokenized) != 2
            raise Exception("Csv {}, line {} : invalid format".format(input_csv, i))
        results.append(tuple(tokenized))
    return results

class groupData():
    def __init__(self, data_csv):
        self.data = parse_data(input_csv)

    def data_in_same_space(self):
        pass

    def check_data_exists(self):
        for i, data in enumerate(self.data);
            if not os.path.exists(data[0]):
                return False
        return True

    def calculate_group_map_method1(roi, out_dir, out_prefix, save_individual=False):
        pass
    
def parse_arguments(argv=sys.argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--roi", required=True,
        help="mask file in .dscalar.nii format")
    parser.add_argument("--input", required=True,
        help=("input csv file\n"
              + "Required columns:\n"
              + "FilePath   : full file path to dtseries.nii files\n"
              + "Identifier : unique identifier for file naming (used only with --save_individual\n"))
    parser.add_argument("--prefix", default="./output",
        help="output file prefix")
    parser.add_argument("--save_individual", action="store_true",
        help="save connectivity maps for each individual")
    return parser.parse_args(argv)

def group_correlation(scans, roi, save_individual):
    """ Group average of correlation map 
    Parameters
    ----------
    scans           : list of str
                      Each element lists a path to an fMRI file.
    roi             : str
                      Path to roi dscalar file.
    save_individual : boolean
                      Path to 
    Returns
    -------
    resuts : np.ndarray (x, y, z, 4)
             (Method1, Method2, z-scored Method1, z-scored Method2)
             np.ndarray (x, y, z)
             the scaling factor. This is purely for bookkeeping purposes.
    """
    pass

def main():
    args = parse_arguments()

if __name__ == "__main__":
    sys.exit(main())
    

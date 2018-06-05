#!/usr/bin/env python3

import os
import argparse
import sys
import time

import nibabel as ni
import numpy as np
import pandas as pd

from nibabel import cifti2 as ci

import cifti_helpers as CH
import math_helpers as MH

def parse_data(input_csv):
    with open(input_csv, "r") as f:
        lines = f.readlines()

    results = []
    for i, line in enumerate(lines):
        tokenized = line.split(",")
        if len(tokenized) != 2:
            raise Exception("Csv {}, line {} : invalid format".format(input_csv, i))
        results.append(tuple(tokenized))
    return results

class groupData():
    def __init__(self, data_csv):
        self.data = parse_data(input_csv)
        for i, data in enumerate(self.data):
            if not os.path.exists(data[0]):
                raise Exception("Line {}, file does not exist: {}".format(i, data[0]))

    def data_in_same_space(self):
        pass

    def calculate_group_map_method1(roi, out_dir, out_prefix, save_individual=False):
        pass
    
def parse_arguments(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument("--t", action="store_true",
        help="print input arugments and exit.")
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

def calc_group_correlation(df, roi, save_individual, out_dir):
    """ Group average of correlation map 
    Parameters
    ----------
    df              : pandas dataframe
                      REQUIRED COLUMNS
                      FilePath   - lists file path to dtseries files
                      Identifier - unique identifier (used for indiual 
                      connectivity map naming)
    roi             : str
                      Path to roi dscalar file.
    save_individual : boolean
                      flag whether to save individual connectivity maps 
    out_dir         : str
                      Path to output directory. This is only used if
                      save_individual is True.
    Returns
    -------
    resuts : np.ndarray (x, y, z, 4)
             (Method1, Method2, z-scored Method1, z-scored Method2)
             np.ndarray (x, y, z)
             the scaling factor. This is purely for bookkeeping purposes.
    """
    rcii = ci.load(roi)
    rdata = np.squeeze(rcii.get_fdata()) != 0

    N = df.shape[0]
    m = np.zeros((4, len(rdata)))
    for i, row in df.iterrows():
        print("Working on image {} of {}: {}".format(i + 1, N, row.FilePath))
        dcii = ci.load(row.FilePath)
        ts = dcii.get_fdata()
        ts_roi = ts[:, rdata]

        cm1 = np.squeeze(MH.pearson_2d(ts.T, ts_roi.T.mean(0, keepdims=True)))
        cm2 = MH.pearson_2d(ts.T, ts_roi.T).mean(1)

        if save_individual:
            for iMethod, val in enumerate((cm1, cm2), start=1): 
                base_file = "{}_r_m{}".format(row.Identifier, iMethod)
                out_file = os.path.join(out_dir, base_file + ".dscalar.nii")
                dscalar = CH.create_dscalar_from_template(rcii, val[np.newaxis], base_file)
                ci.save(dscalar, out_file)
        m[0, :] += cm1
        m[1, :] += cm2
        m[2, :] += np.arctanh(cm1)
        m[3, :] += np.arctanh(cm2)
    m /= N
    return m

def main():
    args = parse_arguments()
    if args.t:
        print(
            "roi             : {}\n".format(args.roi) +
            "input           : {}\n".format(args.input) +
            "prefix          : {}\n".format(args.prefix) +
            "save_individual : {}\n".format(args.save_individual)
        )
        return 0

    # do file checking here
    if not os.path.isfile(args.roi):
        print("roi file does not exist: {}".format(args.roi))
        return -1

    df = pd.read_csv(args.input)
    for i in df.FilePath:
        if not os.path.isfile(i):
            print("input file does not exist: {}".format(i))
            return -1

    out_dir = os.path.dirname(args.prefix)
    out_base = os.path.basename(args.prefix)
    if not os.path.isdir(out_dir):
        print("out directory does not exist: {}".format(out_dir))
        return -1

    conn_out = os.path.join(out_dir, "conn")
    if args.save_individual and not os.path.isdir(conn_out):
        os.mkdir(conn_out)

    m = calc_group_correlation(df, args.roi, args.save_individual, conn_out)
    rcii = ci.load(args.roi)
    map_names = [os.path.join(out_dir, "{}_{}".format(out_base, x))
        for x in ["1r", "2r", "1z", "2z"]]
    out_cifti = CH.create_dscalar_from_template(rcii, m, map_names)
    ci.save(out_cifti, args.prefix + ".dscalar.nii")

    return 0

if __name__ == "__main__":
    sys.exit(main())
    

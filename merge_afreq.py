#!/usr/bin/env python
# coding: utf-8

import os
import glob
import argparse

def get_ID(filename):
    """
    Read tab separated .afreq file and extract ID columns
    """
    ID = ["Pop"]
    with open(filename, "r") as reader:
        lines = reader.readlines()[1:]
        for line in lines:
            line = line.split("\t")
            ID.append(line[1].replace(":", "_"))
    return ID

def get_ALT_FREQS(filename):
    """
    Read tab separated .afreq file and extract ALT_FREQS columns
    """
    pop = filename.split(".")[-2]
    ALT_FREQS = [pop]
    with open(filename, "r") as reader:
        lines = reader.readlines()[1:]
        for line in lines:
            line = line.split("\t")
            ALT_FREQS.append(line[4])
    return ALT_FREQS

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge a folder path.")
    parser.add_argument("dir_path", help="Path to the directory.")
    parser.add_argument("out_path", default="output.txt", help="Output path.")
    args = parser.parse_args()

    filenames = glob.glob(f"{os.path.normpath(args.dir_path)}/*.afreq")

    ID_list = get_ID(filenames[0])
    ALT_FREQS_list = []
    for filename in filenames:
        ALT_FREQS_list.append(get_ALT_FREQS(filename))

    # sort alphabetically
    ALT_FREQS_list = sorted(ALT_FREQS_list, key=lambda x: x[0].upper())

    # write output into file
    with open(args.out_path, "w") as writer:
        writer.write("\t".join(ID_list))
        writer.write("\n")
        for ALT_FREQS in ALT_FREQS_list:
            writer.write("\t".join(ALT_FREQS))
            writer.write("\n")

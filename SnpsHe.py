#!/usr/bin/env python
# coding: utf-8

# read HWE.txt file in read mode
with open("test_hwe.hwe", "r") as f:
    # add chr and pos to results for entry whose obs het > 60
    results = []
    # iterate over each line except first header line
    for line in f.readlines()[1:]:
        # extract chr, pos and obs entry
        Chr, pos, obs = line.strip().split("\t")[:3]
        # extract obs het from obs
        x1, x2, x3 = list(map(int, obs.strip().split("/")))
        het_pcent = (x2*100)/(x1+x2+x3)
        
        # check if obs het percentage > 60
        if het_pcent > 60:
            # write it into results list
            entry = Chr + "\t" + pos
            results.append(entry)

# create a file output.txt and write results into it
with open("SNPShighPropHE.txt", "w") as f:
    # write header CHR and POS
    f.write("CHR\tPOS\n")
    # write each entry of the results list
    for entry in results:
        f.write(entry + "\n")

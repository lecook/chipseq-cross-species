#! /usr/bin/env python

# Author: Adrian Foucal https://gitlab.univ-nantes.fr/foucal-a/full-chipseq/-/blob/master/scripts/count_peaks.py

import json
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--peak_type", help="Required. Peak type for section title. eg broadPeak, narrowpeak...")
parser.add_argument("--sample_name",nargs='+', help="Required. name of your samples")
parser.add_argument("--peaks", nargs = '+', help = "Called peaks, no header")

args = parser.parse_args()


def count_lines(fn):
    return(int(subprocess.check_output("wc -l " + fn, shell=True).split()[0]))

# Go through all sample_name, identify corresponding peak files and get the count in a dict
sample_count = {}
for sample_name in args.sample_name:
    sample_count[sample_name] = {"Peak number": [count_lines(fn) for index, fn in enumerate(args.peaks) if sample_name in fn][0]}


# Get the MultiQC output

multiqc_output = {
    "id": "MACS2_" + args.peak_type + "_count" ,
    "section_anchor": "MACS2 peak",
    "section_name": "MACS2 " + args.peak_type + " counts",
    "description": "",
    "plot_type": "bargraph",
    "pconfig": {
        "id": "MACS2 " + args.peak_type,
        "ylab": "# Peaks",
        "cpswitch": False,                       # Show the 'Counts / Percentages' switch?
        "cpswitch_c_active": True
    },
    "data": sample_count
}

print(json.dumps(multiqc_output, indent = 4))

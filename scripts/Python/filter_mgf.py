import argparse
import re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--spectra_array", required=True, nargs="+")
parser.add_argument("--filter", required=True, help="Regex that will be matched to the files in spectra_array so only those matching will be kept")
args = parser.parse_args()
arguments = vars(args)

match=[re.search(arguments["filter"], sp) for sp in arguments["spectra_array"]]
match=[e is not None for e in match]
match=np.where(np.array(match))[0]
spectra_array=np.array(arguments["spectra_array"])
spectra_array=spectra_array[match]
spectra_array=" ".join(spectra_array)
print(spectra_array)

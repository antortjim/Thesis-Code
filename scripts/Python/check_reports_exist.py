import os.path
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--files", required=True, nargs="+")
parser.add_argument("--root_dir", required=True)
parser.add_argument("--exp_name", required=True)
parser.add_argument("--reports_dir", required=True)
parser.add_argument("--suffix", required=True)
arguments=vars(parser.parse_args())

files=[arguments["root_dir"] + "/" + arguments["exp_name"] + "/peptideShaker_out/" + arguments["reports_dir"] + "/" + e + arguments["suffix"] + ".txt" for e in arguments["files"]]
#print(files)
booleans = np.array([os.path.isfile(e) for e in files])
result=int(np.all(booleans))
print(result)


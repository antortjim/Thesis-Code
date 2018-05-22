import argparse
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--updated_list")
args = parser.parse_args()
arguments = vars(args)
updated_list = arguments["updated_list"]

updated_list = " ".join(["-{} 1".format(e) for e in updated_list.split(" ")])
print(updated_list)


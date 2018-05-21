import argparse
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--updated_list")
parser.add_argument("--se")
parser.add_argument("--searchgui_out")
parser.add_argument("--sample_name")
args = parser.parse_args()
arguments = vars(args)
updated_list = arguments["updated_list"]
extensions = {"msgf": "msgf.mzid", "comet", "comet.pep.xml", "xtandem", "t.xml"}

if os.path.isdir(arguments["searchgui_out"]):
    search_engine_output = os.path.join(arguments["searchgui_out"], arguments["sample_name"]+extensions[arguments["se"]])
    if os.path.exists(search_engine_output):
        updated_list.remove(arguments["se"])

print updated_list





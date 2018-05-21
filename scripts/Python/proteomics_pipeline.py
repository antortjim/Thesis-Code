import subprocess
import numpy as np
import argparse
from argparse import RawTextHelpFormatter


parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)

parser.add_argument("--database_names", default = "Homo_sapiens NZ_Contaminants")
parser.add_argument("--spectra", default = "data/mgf_input")
parser.add_argument("--exp_name", required=True)
parser.add_argument("--searchgui_path", default = "/z/home/aoj/opt/SearchGUI-3.3.1/SearchGUI-3.3.1.jar")
parser.add_argument("--peptideshaker_path", default = "/z/home/aoj/opt/PeptideShaker-1.16.23/PeptideShaker-1.16.23.jar")
parser.add_argument("--params_name", required=False)
parser.add_argument("--ps_out", default = "peptideShaker_out")
parser.add_argument("--settings_dir", default = "settings")
parser.add_argument("--moff_path", default = "/z/home/aoj/opt/moFF/", help="Path to moFF repository")
parser.add_argument("--search_engines", default = "comet", help="list of search_engines separated by comma. Name must be identical to tag used in SearchGUI. Ex comet,msgf,xtandem")
parser.add_argument("--root_dir", default = "/z/home/aoj/thesis/genedata/")

scripts = np.array(["check_flags.sh", "create_decoy_database.sh", "create_settings_file.sh", "search_mgf.sh", "call_peptide_shaker.sh", "fetch_mgf_metadata.sh", "call_moFF.sh"])
help_message = "\n".join(["{}: {}".format(i, script_name) for i, script_name in enumerate(scripts)])
help_message = "String of integers separated by space. Ex \"0 1 2\"\n" + help_message
parser.add_argument("--steps", default = "0", help = help_message)

args = parser.parse_args()
arguments = vars(args)
arguments["steps"] = [int(x) for x in arguments["steps"].split(" ")]
if arguments["params_name"] is None:
    arguments["params_name"] = arguments["exp_name"]

arguments["searchgui_engines"]["-{} 1" for e in arguments["search_engines"]]

scripts = scripts[arguments["steps"]]
arguments["steps"] = " ".join([str(e) for e in arguments["steps"]])

handle = open("{}/pipeline_settings.txt".format(arguments["root_dir"]), "w")
for key, value in arguments.items():
    handle.write("{}:{}\n".format(key.upper(), value.replace(" ", ",")))
handle.close()

print(scripts)
print("Starting pipeline")
for scr in scripts:

    if scr == "call_peptide_shaker.sh":
        continue
    else:
        cmd = r"nohup {}/scripts/bash/{} > {}/{}/{}_{}.out".format(arguments["root_dir"], scr, arguments["root_dir"], arguments["exp_name"], arguments["exp_name"], scr)
        print(cmd)
        subprocess.check_call(cmd, shell=True)
print("Pipeline ended")

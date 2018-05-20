import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--database_names", default = "Homo_sapiens NZ_Contaminants")
parser.add_argument("--spectra", default = "data/mgf_input")
parser.add_argument("--exp_name", required=True)
parser.add_argument("--searchgui_path", default = "/z/home/aoj/opt/SearchGUI-3.3.1/SearchGUI-3.3.1.jar")
parser.add_argument("--peptideshaker_path", default = "/z/home/aoj/opt/PeptideShaker-1.16.23/PeptideShaker-1.16.23.jar")
parser.add_argument("--params_name", required=False)
parser.add_argument("--ps_out", default = "peptideShaker_out")
parser.add_argument("--settings_dir", default = "settings")
parser.add_argument("--root_dir", default = "/z/home/aoj/thesis/genedata/")
parser.add_argument("--steps", default = "0")

args = parser.parse_args()
arguments = vars(args)
arguments["steps"] = [int(x) for x in arguments["steps"].split(" ")]
if arguments["params_name"] is None:
    arguments["params_name"] = arguments["exp_name"]

import subprocess
import numpy as np

scripts = np.array(["check_flags.sh", "create_decoy_database.sh", "create_settings_file.sh", "search_mgf.sh", "call_peptide_shaker.sh"])
scripts = scripts[arguments["steps"]]
arguments["steps"] = " ".join([str(e) for e in arguments["steps"]])

handle = open("{}/pipeline_settings.txt".format(arguments["root_dir"]), "w")
for key, value in arguments.items():
    handle.write("{}:{}\n".format(key.upper(), value.replace(" ", ",")))
handle.close()
#flags  = r"'{}' {} {} {} {} {} {} {} {} &".format(arguments["database_names"], arguments["spectra"], arguments["params_name"],
#                                                  arguments["searchgui_path"], arguments["peptideshaker_path"],
#                                                  arguments["exp_name"], arguments["ps_out"], arguments["settings_dir"], arguments["root_dir"])

print(scripts)
print("start")
for scr in scripts:
    cmd = r"nohup {}/scripts/bash/{} > {}/{}/{}.out".format(arguments["root_dir"], scr, arguments["root_dir"], arguments["exp_name"], arguments["exp_name"])
    print(cmd)
    subprocess.check_call(cmd, shell=True)
print("end")

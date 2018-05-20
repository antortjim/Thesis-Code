import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--database_names", default = "Homo_sapiens NZ_Contaminants")
parser.add_argument("--spectra", default = "data/mgf_input")
parser.add_argument("--params_name", default = "thp1")
parser.add_argument("--searchgui_path", default = "/z/home/aoj/opt/SearchGUI-3.3.1/SearchGUI-3.3.1.jar")
parser.add_argument("--peptideshaker_path", default = "/z/home/aoj/opt/PeptideShaker-1.16.23/PeptideShaker-1.16.23.jar")
parser.add_argument("--exp_name", default = "thp1")
parser.add_argument("--ps_out", default = "peptideShaker_out")
parser.add_argument("--settings_dir", default = "settings")
parser.add_argument("--root_dir", default = "/z/home/aoj/thesis/genedata/thp1")

args = parser.parse_args()
arguments = vars(args)

import subprocess
print("start")
cmd = r"nohup {}/proteomics_pipeline.sh '{}' {} {} {} {} {} {} {} {} &".format(arguments["root_dir"], arguments["database_names"], arguments["spectra"], arguments["params_name"],
                                                        arguments["searchgui_path"], arguments["peptideshaker_path"],
                                                        arguments["exp_name"], arguments["ps_out"], arguments["settings_dir"], arguments["root_dir"])
print(cmd)
subprocess.check_call(cmd, shell=True)
print("end")

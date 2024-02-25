#!/usr/bin/env python
import json
import sys
import os
import re
from os.path import basename
import SharedConsts
from script.guidance_scoring_and_visualization import remove_low_sp_sites_no_bioperl
from script.guidance_sequence_functions import codes2name_fasta_from1


def print_message_to_output(msg):
    with open(output_page, "a") as output_file:
        output_file.write(f"\n<ul><li>{msg}</li></ul>\n")

VARS = {}
if len(sys.argv) < 2:
    sys.exit("USAGE: perl {} --MSA <Base MSA> --Scores <Scores_File> --FilterdMSA <Out_MSA_File> --Cutoff <Cutoff> --RemovedPos <Out_File_With_Removed_Pos>\n".format(sys.argv[0]))

# Command line mode
if sys.argv[1].startswith("--"):
    if len(sys.argv) < 11:
        sys.exit("USAGE: python3 {} --MSA <Base MSA> --Scores <Scores_File> --FilterdMSA <Out_MSA_File> --Cutoff <Cutoff> --RemovedPos <Out_File_With_Removed_Pos>\n".format(sys.argv[0]))
    # Extract options
    options = {sys.argv[i]: sys.argv[i+1] for i in range(1, len(sys.argv), 2)}
    is_server = "NO"
    cutoff = options["--Cutoff"]
    VARS['Alignment_File_without_low_SP_Col'] = options["--FilterdMSA"]
    VARS['Col_Scores_File'] = options["--Scores"]
    VARS['Alignment_File'] = options["--MSA"]
    VARS['removed_low_SP_SITE'] = options['--RemovedPos']

# Server mode
else:
    stored_data_file = sys.argv[1]
    cutoff = float(sys.argv[2])
    # with open(stored_data_file, 'rb') as f:
    #     VARS = pickle.load(f)
    with open(stored_data_file, 'r') as vars_file:
        json_string = vars_file.read()
        VARS = json.loads(json_string)

    # Set file paths
    VARS['Alignment_File_without_low_SP_Col'] = VARS["WorkingDir"] + "/" + VARS["Alignment_File_without_low_SP_Col"] + f".{cutoff}"
    VARS['Col_Scores_File'] = VARS["WorkingDir"] + "/" + VARS["Output_Prefix"] + "_res_pair_col.scr"
    VARS['Alignment_File'] = VARS['WorkingDir'] + VARS['Alignment_File']
    VARS['removed_low_SP_SITE'] = VARS['WorkingDir']+VARS['removed_low_SP_SITE'] + f".{cutoff}"
    is_server = "YES"


# Remove sites with SP-score < Cutoff
if is_server == "YES":
    log_file = VARS["OutLogFile"]
    try:
        with open(log_file, 'a') as log_file:
            log_file.write(f"Guidance::removeLowSPsites_NoBioPerl ('{VARS['Alignment_File']}', '{VARS['Col_Scores_File']}', '{VARS['Alignment_File_without_low_SP_Col']}', {cutoff}, '{VARS['removed_low_SP_SITE']}');\n")
    except Exception as e:
        print(f"Can't open Log File: {log_file} {e} \n")
        sys.exit()

ans = remove_low_sp_sites_no_bioperl(f"{VARS['Alignment_File']}",
                                     f"{VARS['Col_Scores_File']}",
                                     f"{VARS['Alignment_File_without_low_SP_Col']}",
                                     cutoff,
                                     f"{VARS['removed_low_SP_SITE']}")

if ans[0] == "OK":
    VARS['REMOVED_SITES'] = ans[1]
    VARS['MSA_LENGTH'] = ans[2]

if is_server == "NO":
    print(f"REMOVED_SITES:{VARS['REMOVED_SITES']}\n")
    print(f"MSA_LENGTH:{VARS['MSA_LENGTH']}\n")

# server
else:
    log_file = VARS["OutLogFile"]
    with open(log_file, "a") as log_file:
        log_file.write(f"REMOVED_SITES:{VARS['REMOVED_SITES']}\n")
        log_file.write(f"MSA_LENGTH:{VARS['MSA_LENGTH']}\n")

    VARS['Alignment_File_without_low_SP_Col_with_Names'] = f"{VARS['Alignment_File_without_low_SP_Col']}.With_Names"
    if os.path.getsize(f"{VARS['Alignment_File_without_low_SP_Col']}") > 0:  # Not EMPTY
        ans = codes2name_fasta_from1(f"{VARS['Alignment_File_without_low_SP_Col']}",
                                     f"{VARS['WorkingDir']+VARS['code_fileName']}",
                                     f"{VARS['Alignment_File_without_low_SP_Col_with_Names']}")
        if ans[0] != "OK":
            log_file.write(f"Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(\"{VARS['Alignment_File_without_low_SP_Col']}\",\"{VARS['WorkingDir']+ VARS['code_fileName']}\",\"{VARS['Alignment_File_without_low_SP_Col_with_Names']}\") failed: {''.join(ans)}\n",)

    # Update the output page
    #######################################
    output_page = VARS["WorkingDir"] + "/" + VARS["output_page"]

    with open(output_page, "r") as output:
        out_lines = output.readlines()

    with open(output_page, "w") as output:
        remove_pos_section = 0
        for line in out_lines:
            if "Remove unreliable columns below confidence score" in line:
                remove_pos_section = 1
                output.write(line)
            elif "form" in line and remove_pos_section == 1:
                output.write(line)
                Alignment_File_without_low_SP_Col_with_Names_NO_PATH = os.path.basename(
                    VARS["Alignment_File_without_low_SP_Col_with_Names"])
                removed_low_SP_SITE_NO_PATH = os.path.basename(VARS["removed_low_SP_SITE"])

                print_message_to_output(f"<A HREF='{Alignment_File_without_low_SP_Col_with_Names_NO_PATH}' TARGET=_blank>The MSA after the removal of unreliable columns (below {cutoff})</A><font size=-1> (see list of removed columns <A HREF='{removed_low_SP_SITE_NO_PATH}' TARGET=_blank>here</A>)</font><br>")
                remove_pos_section = 0
            else:
                output.write(line)


#!/usr/bin/env python
import json
import sys
import os
from os.path import basename

from script.guidance_scoring_and_visualization import remove_low_sp_seq
from script.guidance_sequence_functions import codes2name_fasta_from1


def print_message_to_output(msg):
    with open(output_page, "a") as output_file:
        output_file.write(f"\n<ul><li>{msg}</li></ul>\n")

VARS = {}

if len(sys.argv) < 3:
    sys.exit("USAGE: python3 {} --MSA <Base MSA> --Scores <Scores_File> --FilterdSeq <Out_Seq_File> --Cutoff <Cutoff> --RemovedSeq <Out_File_With_Removed_Seq> --Type <BySeqName|ByRowNum [optional]>\n".format(sys.argv[0]))

# command line mode
if sys.argv[1].startswith("--"):
    if len(sys.argv) < 11:
        sys.exit("USAGE: python3 {} --MSA <Base MSA> --Scores <Scores_File> --FilterdSeq <Out_Seq_File> --Cutoff <Cutoff> --RemovedSeq <Out_File_With_Removed_Seq> --Type <BySeqName|ByRowNum [optional]>\n".format(sys.argv[0]))
    options = {sys.argv[i]: sys.argv[i+1] for i in range(1, len(sys.argv), 2)}

    VARS["Alignment_File"] = options["--MSA"]
    VARS["Seq_Scores_File"] = options["--Scores"]
    VARS["Seq_File_without_low_SP_SEQ"] = options["--FilterdSeq"]
    cutoff = options["--Cutoff"]
    VARS["removed_low_SP_SEQ"] = options["--RemovedSeq"]
    if "--Type" not in options:
        type = "BySeqName"
    else:
        type = options["--Type"]
    is_server = "NO"

# server mode
else:
    stored_data_file = sys.argv[1]
    stored_form_file = sys.argv[2]
    cutoff = float(sys.argv[3])

    with open(stored_data_file, 'r') as vars_file:
        json_string = vars_file.read()
        VARS = json.loads(json_string)

    with open(stored_form_file, 'r') as form_file:
        json_string = form_file.read()
        FORM = json.loads(json_string)

    VARS["Seq_Scores_File"] = f"{VARS['WorkingDir']}{VARS['Output_Prefix']}_res_pair_seq.scr"
    VARS["Seq_File_without_low_SP_SEQ"] = f"{VARS['WorkingDir']}{VARS['Seq_File_without_low_SP_SEQ']}.{cutoff}"
    VARS["removed_low_SP_SEQ"] = f"{VARS['WorkingDir']}/{VARS['removed_low_SP_SEQ']}.{cutoff}"
    is_server = "YES"
    type = "ByRowNum"


#remove sites with SP-score < Col sp_cutoff
############################################

if is_server == "YES":
    log_file = VARS["OutLogFile"]
    try:
        with open(log_file, "a") as log_file:
            log_file.write(f"Guidance::removeLowSPseq({VARS['WorkingDir'] + VARS['Alignment_File']}, {VARS['Seq_Scores_File']}, {VARS['Seq_File_without_low_SP_SEQ']}, {cutoff}, {VARS['removed_low_SP_SEQ']}, {type});\n")
    except Exception as e:
        print(f"Can't open Log File: {log_file} {e} \n")
        sys.exit()

ans = remove_low_sp_seq(VARS['WorkingDir'] + VARS['Alignment_File'], VARS['Seq_Scores_File'], VARS['Seq_File_without_low_SP_SEQ'], cutoff, VARS['removed_low_SP_SEQ'], type)
print("ANS: ", "".join(ans), "\n")

if is_server == "YES":
    VARS["Seq_File_without_low_SP_SEQ_with_Names"] = f"{VARS['Seq_File_without_low_SP_SEQ']}.With_Names"
    VARS["removed_low_SP_SEQ_With_Names"] = f"{VARS['removed_low_SP_SEQ']}.With_Names"

    if os.path.getsize(VARS["Seq_File_without_low_SP_SEQ"]) > 0:
        ans = codes2name_fasta_from1(VARS["Seq_File_without_low_SP_SEQ"], f"{VARS['WorkingDir']}{VARS['code_fileName']}", VARS["Seq_File_without_low_SP_SEQ_with_Names"])
        if ans[0] != "OK":
            log_file.write(f"Guidance::codes2name_fasta_from1({VARS['Seq_File_without_low_SP_SEQ']}, \"{VARS['WorkingDir']}{VARS['code_fileName']}\", {VARS['Seq_File_without_low_SP_SEQ_with_Names']}) failed: {' '.join(ans)}\n")
            sys.exit()

    if os.path.getsize(VARS["removed_low_SP_SEQ"]) > 0:
        ans = codes2name_fasta_from1(VARS["removed_low_SP_SEQ"], f"{VARS['WorkingDir']}{VARS['code_fileName']}", VARS["removed_low_SP_SEQ_With_Names"])
        if ans[0] != "OK":
            log_file.write(f"Guidance::codes2nameFastaFrom1({VARS['removed_low_SP_SEQ']}, \"{VARS['WorkingDir']}{VARS['code_fileName']}\", {VARS['removed_low_SP_SEQ_With_Names']}) failed: {' '.join(ans)}\n")
            sys.exit()

    # Update the output page
    #######################################
    output_page = VARS["WorkingDir"] + "/" + VARS["output_page"]

    with open(output_page, "r") as output:
        out_lines = output.readlines()

    with open(output_page, "w") as output:
        remove_seq_section = 0
        Seq_File_without_low_SP_SEQ_with_Names_NO_PATH = basename(
            VARS["Seq_File_without_low_SP_SEQ_with_Names"])
        removed_low_SP_SEQ_With_Names_NO_PATH = basename(VARS["removed_low_SP_SEQ_With_Names"])

        for line in out_lines:
            if "Remove unreliable sequences below confidence score" in line:
                remove_seq_section = 1
                output.write(line)
            elif "form" in line and remove_seq_section == 1 and "form.data" not in line:
                output.write(line)
                if FORM['Redirect_From_MAFFT'] == 1:
                    print_message_to_output(f"<A HREF='{Seq_File_without_low_SP_SEQ_with_Names_NO_PATH}' TARGET=_blank>The input sequences after the removal of unreliable sequences (with confidence score below {Cutoff})</A><font size=-1> (see list of removed sequences <A HREF='{removed_low_SP_SEQ_With_Names_NO_PATH}' TARGET=_blank>here</A></font>)&nbsp;&nbsp;&nbsp;<INPUT TYPE=\"BUTTON\" VALUE=\"run GUIDANCE on the confidently-aligned sequences only\" ONCLICK=\"var answer = confirm('ATTENTION: Running GUIDANCE on the confidently-aligned sequences only, ignores the parameters used for the original run on MAFFT server. It is therefore recommended to adjust these parameters or aligning the confidently-aligned sequences on MAFFT server and run GUIDANCE again from there');if (answer){{window.open('http://guidance.tau.ac.il/index_rerun.php?run={VARS['run_number']}&file={Seq_File_without_low_SP_SEQ_with_Names_NO_PATH}')}}\"><br>")
                else:
                    print_message_to_output(f"<A HREF='{Seq_File_without_low_SP_SEQ_with_Names_NO_PATH}' TARGET=_blank>The input sequences after the removal of unreliable sequences (with confidence score below {Cutoff})</A><font size=-1> (see list of removed sequences <A HREF='{removed_low_SP_SEQ_With_Names_NO_PATH}' TARGET=_blank>here</A></font>)&nbsp;&nbsp;&nbsp;<INPUT TYPE=\"BUTTON\" VALUE=\"run GUIDANCE on the confidently-aligned sequences only\" ONCLICK=\"window.open('http://guidance.tau.ac.il/index_rerun.php?run={VARS['run_number']}&file={Seq_File_without_low_SP_SEQ_with_Names_NO_PATH}')\"><br>")
            else:
                output.write(line)


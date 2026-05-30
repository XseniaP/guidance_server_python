"""Sequence filtering utilities: mask low-score residues, remove low-SP columns/sequences."""
import json
import math
import os
import sys
from os.path import basename
from warnings import warn

import Bio.SeqIO as SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from guidance3.pipeline.scoring import remove_low_sp_sites_no_bioperl, remove_low_sp_seq
from guidance3.sequences.fasta import codes2name_fasta_from1


def _print_to_output(output_page, msg):
    with open(output_page, "a") as f:
        f.write(f"\n<ul><li>{msg}</li></ul>\n")


def mask_residues_main(argv=None):
    """Mask low-score residues in an MSA — called by the Flask mask route."""
    if argv is None:
        argv = sys.argv

    if len(argv) < 4:
        sys.exit(f"USAGE: {argv[0]} STORED_DATA_FILE ALPHABET CUTOFF\n"
                 "ALPHABET can be either aa or nuc\n")

    stored_data_file, alphabet, cutoff = argv[1], argv[2], argv[3]

    with open(stored_data_file) as f:
        VARS = json.loads(f.read())

    msa_file = os.path.join(VARS["WorkingDir"], VARS["Alignment_File"])
    score_file = os.path.join(VARS["WorkingDir"], f"{VARS['Output_Prefix']}_res_pair_res.scr")
    out_website = f"Mask_Residues_Res_{cutoff}.aln"
    out_file = os.path.join(VARS["WorkingDir"], out_website)
    seq_names_index = os.path.join(VARS["WorkingDir"], "Seqs.Codes")
    output_page = os.path.join(VARS["WorkingDir"], VARS["output_page"])

    if alphabet == "aa":
        missing_data_char = "X"
    elif alphabet == "nuc":
        missing_data_char = "N"
    else:
        sys.exit("ALPHABET must be either 'aa' or 'nuc'\n")

    seqs, ids = [], []
    with open(msa_file) as fasta_file:
        for rec in SeqIO.parse(fasta_file, "fasta"):
            seqs.append(list(str(rec.seq)))
            ids.append(rec.id)

    with open(score_file) as in_file:
        for line in in_file:
            line = line.strip()
            if line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) == 3:
                col, row, score = int(cols[0]), int(cols[1]), float(cols[2])
                if not math.isnan(score) and score < float(cutoff):
                    seqs[row - 1][col - 1] = missing_data_char
            else:
                warn(f"WARNING: failed to parse line: '{line}'\n")

    id_names = {}
    with open(seq_names_index) as seq_index:
        for line in seq_index:
            name, seq_id = line.strip().split("\t")
            id_names[seq_id] = name

    with open(out_file, "w") as out:
        for i, seq_chars in enumerate(seqs):
            seq = "".join(seq_chars)
            original_name = id_names.get(ids[i], ids[i])
            SeqIO.write(SeqRecord(Seq(seq), id=original_name, description=""), out, "fasta")

    os.chmod(out_file, 0o664)

    with open(output_page) as f:
        out_lines = f.readlines()
    with open(output_page, "w") as f:
        for line in out_lines:
            f.write(line)
            if "Mask specific residues below a certain cutoff:" in line:
                os.chmod(out_file, 0o664)
                _print_to_output(output_page,
                                 f"<A HREF={out_website} TARGET=_blank>The MSA after masking unreliable residues (below {cutoff})</A>")


def remove_pos_main(argv=None):
    """Remove low-SP-score columns — called by the Flask remove_pos route."""
    if argv is None:
        argv = sys.argv

    VARS = {}
    if len(argv) < 2:
        sys.exit("USAGE: python3 {} VARS_JSON CUTOFF\n".format(argv[0]))

    if argv[1].startswith("--"):
        if len(argv) < 11:
            sys.exit("USAGE: python3 {} --MSA <MSA> --Scores <Scores> --FilterdMSA <Out> --Cutoff <C> --RemovedPos <Out2>\n".format(argv[0]))
        options = {argv[i]: argv[i + 1] for i in range(1, len(argv), 2)}
        is_server = "NO"
        cutoff = options["--Cutoff"]
        VARS['Alignment_File_without_low_SP_Col'] = options["--FilterdMSA"]
        VARS['Col_Scores_File'] = options["--Scores"]
        VARS['Alignment_File'] = options["--MSA"]
        VARS['removed_low_SP_SITE'] = options['--RemovedPos']
    else:
        stored_data_file = argv[1]
        cutoff = float(argv[2])
        with open(stored_data_file) as f:
            VARS = json.loads(f.read())
        VARS['Alignment_File_without_low_SP_Col'] = (
            VARS["WorkingDir"] + "/" + VARS["Alignment_File_without_low_SP_Col"] + f".{cutoff}")
        VARS['Col_Scores_File'] = VARS["WorkingDir"] + "/" + VARS["Output_Prefix"] + "_res_pair_col.scr"
        VARS['Alignment_File'] = VARS['WorkingDir'] + VARS['Alignment_File']
        VARS['removed_low_SP_SITE'] = VARS['WorkingDir'] + VARS['removed_low_SP_SITE'] + f".{cutoff}"
        VARS['code_fileName'] = 'Seqs.Codes'
        is_server = "YES"

    if is_server == "YES":
        try:
            with open(VARS["OutLogFile"], 'a') as log:
                log.write(f"remove_pos: ({VARS['Alignment_File']}, {VARS['Col_Scores_File']}, "
                          f"{VARS['Alignment_File_without_low_SP_Col']}, {cutoff}, {VARS['removed_low_SP_SITE']})\n")
        except Exception as e:
            print(f"Can't open log: {e}")
            sys.exit()

    ans = remove_low_sp_sites_no_bioperl(
        VARS['Alignment_File'], VARS['Col_Scores_File'],
        VARS['Alignment_File_without_low_SP_Col'], cutoff, VARS['removed_low_SP_SITE'])

    if ans[0] == "OK":
        VARS['REMOVED_SITES'] = ans[1]
        VARS['MSA_LENGTH'] = ans[2]

    if is_server == "NO":
        print(f"REMOVED_SITES:{VARS['REMOVED_SITES']}\nMSA_LENGTH:{VARS['MSA_LENGTH']}\n")
        return

    with open(VARS["OutLogFile"], "a") as log:
        log.write(f"REMOVED_SITES:{VARS['REMOVED_SITES']}\nMSA_LENGTH:{VARS['MSA_LENGTH']}\n")

    VARS['Alignment_File_without_low_SP_Col_with_Names'] = f"{VARS['Alignment_File_without_low_SP_Col']}.With_Names"
    if os.path.getsize(VARS['Alignment_File_without_low_SP_Col']) > 0:
        ans = codes2name_fasta_from1(
            VARS['Alignment_File_without_low_SP_Col'],
            VARS['WorkingDir'] + VARS['code_fileName'],
            VARS['Alignment_File_without_low_SP_Col_with_Names'])

    output_page = VARS["WorkingDir"] + "/" + VARS["output_page"]
    with open(output_page) as f:
        out_lines = f.readlines()
    with open(output_page, "w") as f:
        remove_pos_section = 0
        for line in out_lines:
            if "Remove unreliable columns below confidence score" in line:
                remove_pos_section = 1
                f.write(line)
            elif "form" in line and remove_pos_section == 1:
                f.write(line)
                col_no_path = basename(VARS["Alignment_File_without_low_SP_Col_with_Names"])
                rem_no_path = basename(VARS["removed_low_SP_SITE"])
                _print_to_output(output_page,
                                 f"<A HREF='{col_no_path}' TARGET=_blank>The MSA after removing unreliable columns "
                                 f"(below {cutoff})</A><font size=-1> (see removed columns "
                                 f"<A HREF='{rem_no_path}' TARGET=_blank>here</A>)</font><br>")
                remove_pos_section = 0
            else:
                f.write(line)


def remove_seq_main(argv=None):
    """Remove low-SP-score sequences — called by the Flask remove_seq route."""
    if argv is None:
        argv = sys.argv

    VARS = {}
    if len(argv) < 3:
        sys.exit("USAGE: python3 {} VARS_JSON FORM_JSON CUTOFF\n".format(argv[0]))

    if argv[1].startswith("--"):
        options = {argv[i]: argv[i + 1] for i in range(1, len(argv), 2)}
        VARS["Alignment_File"] = options["--MSA"]
        VARS["Seq_Scores_File"] = options["--Scores"]
        VARS["Seq_File_without_low_SP_SEQ"] = options["--FilterdSeq"]
        cutoff = options["--Cutoff"]
        VARS["removed_low_SP_SEQ"] = options["--RemovedSeq"]
        seq_type = options.get("--Type", "BySeqName")
        is_server = "NO"
    else:
        stored_data_file, stored_form_file = argv[1], argv[2]
        cutoff = float(argv[3])
        with open(stored_data_file) as f:
            VARS = json.loads(f.read())
        with open(stored_form_file) as f:
            FORM = json.loads(f.read())
        VARS["Seq_Scores_File"] = f"{VARS['WorkingDir']}{VARS['Output_Prefix']}_res_pair_seq.scr"
        VARS["Seq_File_without_low_SP_SEQ"] = f"{VARS['WorkingDir']}{VARS['Seq_File_without_low_SP_SEQ']}.{cutoff}"
        VARS["removed_low_SP_SEQ"] = f"{VARS['WorkingDir']}/{VARS['removed_low_SP_SEQ']}.{cutoff}"
        VARS["code_fileName"] = "Seqs.Codes"
        is_server = "YES"
        seq_type = "ByRowNum"

    if is_server == "YES":
        try:
            with open(VARS["OutLogFile"], "a") as log:
                log.write(f"remove_seq: ({VARS['WorkingDir'] + VARS['Alignment_File']}, "
                          f"{VARS['Seq_Scores_File']}, {VARS['Seq_File_without_low_SP_SEQ']}, "
                          f"{cutoff}, {VARS['removed_low_SP_SEQ']}, {seq_type})\n")
        except Exception as e:
            print(f"Can't open log: {e}")
            sys.exit()

    ans = remove_low_sp_seq(
        VARS['WorkingDir'] + VARS['Alignment_File'], VARS['Seq_Scores_File'],
        VARS['Seq_File_without_low_SP_SEQ'], cutoff, VARS['removed_low_SP_SEQ'], seq_type)
    print("ANS:", "".join(ans))

    if is_server != "YES":
        return

    VARS["Seq_File_without_low_SP_SEQ_with_Names"] = f"{VARS['Seq_File_without_low_SP_SEQ']}.With_Names"
    VARS["removed_low_SP_SEQ_With_Names"] = f"{VARS['removed_low_SP_SEQ']}.With_Names"

    if os.path.getsize(VARS["Seq_File_without_low_SP_SEQ"]) > 0:
        codes2name_fasta_from1(VARS["Seq_File_without_low_SP_SEQ"],
                               f"{VARS['WorkingDir']}{VARS['code_fileName']}",
                               VARS["Seq_File_without_low_SP_SEQ_with_Names"])

    if os.path.getsize(VARS["removed_low_SP_SEQ"]) > 0:
        codes2name_fasta_from1(VARS["removed_low_SP_SEQ"],
                               f"{VARS['WorkingDir']}{VARS['code_fileName']}",
                               VARS["removed_low_SP_SEQ_With_Names"])

    output_page = VARS["WorkingDir"] + "/" + VARS["output_page"]
    with open(output_page) as f:
        out_lines = f.readlines()
    with open(output_page, "w") as f:
        remove_seq_section = 0
        seq_no_path = basename(VARS["Seq_File_without_low_SP_SEQ_with_Names"])
        rem_no_path = basename(VARS["removed_low_SP_SEQ_With_Names"])
        for line in out_lines:
            if "Remove unreliable sequences below confidence score" in line:
                remove_seq_section = 1
                f.write(line)
            elif "form" in line and remove_seq_section == 1 and "form.data" not in line:
                f.write(line)
                _print_to_output(output_page,
                                 f"<A HREF='{seq_no_path}' TARGET=_blank>The input sequences after removing "
                                 f"unreliable sequences (below {cutoff})</A>"
                                 f"<font size=-1> (see removed sequences "
                                 f"<A HREF='{rem_no_path}' TARGET=_blank>here</A></font>)<br>")
                remove_seq_section = 0
            else:
                f.write(line)


def mask_residues_simple(msa_file, score_file, out_file, cutoff, alphabet):
    """Mask low-score residues given direct file paths (no JSON config, no HTML output).

    Replacement characters: 'X' for amino acids, 'N' for nucleotides.
    """
    if alphabet == "aa":
        missing_data_char = "X"
    elif alphabet == "nuc":
        missing_data_char = "N"
    else:
        sys.exit("ALPHABET must be either 'aa' or 'nuc'\n")

    seqs, ids = [], []
    with open(msa_file) as f:
        for rec in SeqIO.parse(f, "fasta"):
            seqs.append(list(str(rec.seq)))
            ids.append(rec.id)

    with open(score_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) == 3:
                col, row, score = int(cols[0]), int(cols[1]), float(cols[2])
                if not math.isnan(score) and score < float(cutoff):
                    seqs[row - 1][col - 1] = missing_data_char

    with open(out_file, "w") as f:
        for i, seq_chars in enumerate(seqs):
            SeqIO.write(SeqRecord(Seq("".join(seq_chars)), id=ids[i], description=""), f, "fasta")


def mask_residues_simple_main(argv=None):
    """Entry point for guidance3-mask-simple console script."""
    if argv is None:
        argv = sys.argv
    if len(argv) < 6:
        sys.exit(f"USAGE: {argv[0]} CODON_MSA_FILE GUIDANCE_RESIDUE_SCORES_FILE OUT_FILE CUTOFF ALPHABET\n"
                 "ALPHABET can be either aa or nuc\n")
    mask_residues_simple(argv[1], argv[2], argv[3], argv[4], argv[5])


if __name__ == "__main__":
    # Allow running as: python -m guidance3.sequences.filters mask|remove_pos|remove_seq|mask_simple [args]
    if len(sys.argv) < 2:
        sys.exit("USAGE: python -m guidance3.sequences.filters <mask|remove_pos|remove_seq|mask_simple> [args...]\n")
    mode = sys.argv.pop(1)
    if mode == "mask":
        mask_residues_main()
    elif mode == "remove_pos":
        remove_pos_main()
    elif mode == "remove_seq":
        remove_seq_main()
    elif mode == "mask_simple":
        mask_residues_simple_main()
    else:
        sys.exit(f"Unknown mode: {mode}\n")

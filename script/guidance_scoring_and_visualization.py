import math
import os

from matplotlib import pyplot as plt

from guidance_sequence_functions import *
import numpy as np

from time_decorator import timeit


def remove_pos_no_bioperl(pos_to_remove, msa_hash_ref):
    for seq_id, seq_list in msa_hash_ref.items():
        seq_list[pos_to_remove] = ""
    return msa_hash_ref


def print_msa_hash(ref_to_hash):
    for key, seq_list in ref_to_hash.items():
        print(f">{key}")
        print("".join(seq_list))


def read_msa_to_hash(msa_file):
    msa_hash = {}
    msa_length = 0
    msa_order = []

    with open(msa_file, 'r') as in_msa:
        fasta_line = in_msa.readline()
        while fasta_line:
            header = fasta_line[1:].strip()
            fasta_line = in_msa.readline()

            seq = ""
            while fasta_line and fasta_line[0] != ">":
                seq += fasta_line.strip()
                fasta_line = in_msa.readline()

            seq_arr = list(seq)
            if header in msa_hash:
                print(f"[WARNING] The sequence name '{header}' appears more than once in the MSA '{msa_file}'. "
                      "Consider only the first instance....\n")
            else:
                msa_hash[header] = seq_arr
                msa_order.append(header)

            if msa_length == 0:
                msa_length = len(seq)

    return "OK", msa_hash, msa_length, msa_order

@timeit
def remove_low_sp_sites_no_bioperl(msa_file, sp_file, out_file, cutoff, pos_removed_file):
    num_pos_removed = 0

    ans = read_msa_to_hash(msa_file)
    if ans[0] != "OK":
        joined_answer = " ".join(ans)
        return f"removeLowSPsites error: {joined_answer}\n"
    msa_hash_ref, msa_length, msa_order_array_ref = ans[1], ans[2], ans[3]
    try:
        with open(sp_file, 'r') as in_file, open(out_file, 'w') as out_file_handle:
            if pos_removed_file != "":
                with open(pos_removed_file, 'w') as removed_pos_file:
                    header = in_file.readline()
                    for line in in_file:
                        line = line.strip()
                        match = re.match(r'^\s*(\d+)\s+(\d+(\.\d+)?)', line)
                        if match:
                            pos, score = int(match.group(1)), float(match.group(2))
                            if score < float(cutoff):
                                msa_hash_ref = remove_pos_no_bioperl(pos - 1,
                                                                     msa_hash_ref)  # Pos-1 because array starts from 0
                                removed_pos_file.write(f'Remove Pos: {pos}\tScore: {score}\n')
                                num_pos_removed += 1

                    for seq_id in msa_order_array_ref:
                        seq = ''.join(msa_hash_ref[seq_id])
                        if seq.startswith('-' * len(seq)):
                            print(
                                "WARNING: After removing positions scored below {}, the sequence {} comprises only gap characters...\n".format(
                                    cutoff, seq_id))
                        out_file_handle.write(f'>{seq_id}\n{seq}\n')

        return 'OK', num_pos_removed, msa_length

    except Exception as e:
        print(f"remove_low_sp_sites_no_bioperl: exception occurred: {e}\n")
        sys.exit()


def calculate_mean_and_std(data_file, column):
    # get a delimited file and calculate the std and average

    column -= 1  # Adjust for zero-based indexing
    data = []

    try:
        with open(data_file, 'r') as data_file:
            header = next(data_file)  # Skip the header
            for line in data_file:
                line = line.strip()
                values = line.split()
                if len(values) > column and values[column].isdigit():
                    data.append(float(values[column]))
    except Exception as e:
        sys.exit("Error: exception occurred while opening file in calculate_mean_and_std: " + str(e) + "\n")

    if data:
        avg = np.mean(data)
        std = np.std(data)
    else:
        avg = np.nan
        std = np.nan

    return avg, std

@timeit
def remove_low_sp_sites_consider_z(msa_file, sp_file, out_file, cutoff, z_cutoff, pos_removed_file):
    num_removed_pos = 0
    result, msa_hash_ref, msa_length, msa_order_array_ref = read_msa_to_hash(msa_file)
    if result != "OK":
        return f"removeLowSPsites_Consider_Z:{result}\n"

    mean, std = calculate_mean_and_std(sp_file, 2)

    with open(sp_file, 'r') as in_file, open(out_file, 'w') as out_file:
        if pos_removed_file != "" or pos_removed_file != None:
            with open(pos_removed_file, 'w') as removed_pos:
                num_removed_pos = 0
                for line in in_file:
                    match = re.match(r'\s*(\d+)\s+(\d+(\.\d+)?)', line)
                    if match:
                        site_num, site_score = map(float, match.groups())
                        z_score = "NaN" if std == 0 else (site_score - mean) / std

                        if z_score != "NaN" and z_score < z_cutoff:
                            if site_score < cutoff:
                                msa_hash_ref = remove_pos_no_bioperl(int(site_num) - 1, msa_hash_ref)
                                removed_pos.write(f"Remove Pos: {site_num}\tScore: {site_score}\tZ_Score: {z_score}\n")
                                num_removed_pos += 1

                for seq_id in msa_order_array_ref:
                    seq = ''.join(msa_hash_ref[seq_id])
                    if seq.startswith('-'):
                        print(f"WARNING: After removing positions scored below {cutoff} and Z_Score: {z_cutoff}, "
                              f"the sequence {seq_id} comprised of only gap characters...\n")

                    out_file.write(f">{seq_id}\n{seq}\n")

    return "OK", num_removed_pos, msa_length

@timeit
def convert_to_csv(input_file, output_file):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            # header = infile.readline()
            for line in infile:
                line = line.strip()
                line = trim(line)
                elements = line.split()
                new_line = ",".join(elements)
                line = line.replace("#", "")
                if "END" not in line:
                    outfile.write(new_line + "\n")
        return ["OK"]
    except Exception as e:
        return f"Error: {str(e)}"

@timeit
def print_colored_alignment_with_css(in_msa_file, out_html_file, scores_file, codes_file, col_scores_csv, x_label,
                                     seq_scores):
    with open(seq_scores, 'r') as seq_scores_file:
        seq_scores_data = {}
        for line in seq_scores_file:
            if not line.startswith('#'):
                seq, score = line.split()
                seq_scores_data[seq] = score

    with open(scores_file, 'r') as scores_file:
        scores_data = {}
        code_names = {}
        for line in scores_file:
            if not line.startswith('#'):
                line.strip()
                line = trim(line)
                col, seq, score = line.split()
                if seq not in scores_data:
                    scores_data[seq] = {}
                scores_data[seq][col] = score

        if codes_file != "":
            try:
                with open(codes_file, 'r') as codes_file:
                    for line in codes_file:
                        seq_name, code = line.strip().split('\t')
                        code_names[code] = seq_name
            except Exception as e:
                return f"Guidance::printColoredAlignment Can't open the Codes file: '{codes_file}' {e}"

    try:
        alignment = AlignIO.read(in_msa_file, 'fasta')
    except Exception as e:
        sys.exit(f"Can't open {in_msa_file}: {e}\n")

    alignment.verbose = True  # HAIM COMMNET
    # Otherwise, bioperl adds sequence start/stop values
    set_displayname_flat(alignment)  # HAIM COMMENT
    ans = check_msa_licit_and_size(in_msa_file, "fasta", "no")
    if ans[0] == "OK":
        msa_depth = ans[1]
    else:
        answer_joined = " ".join(ans)
        return f"printColoredAlignment:  {answer_joined}\n"

    # Print HTML start
    # Code from Conseq colored MSA: ~/pupkoSVN/trunk/www/conseq/runCalc_Conseq.pl line 985

    # msa_depth = len(alignment)
    sequence_length_for_display = 400000
    fontsize = 2
    msaRightOrderCounter = 0
    tdWidth = 5
    msaRightOrder = []
    colorstep = ["Score1", "Score2", "Score3", "Score4", "Score5", "Score6", "Score7", "Score8", "Score9", "Score9"]
    colorstep_code = ["#10C8D1", "#8CFFFF", "#D7FFFF", "#EAFFFF", "#FFFFFF", "#FCEDF4", "#FAC9DE", "#F07DAB", "#A02560"]

    try:
        with open(out_html_file, "w") as msa_colored_html:
            msa_colored_html.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\"\n")
            msa_colored_html.write("\"http://www.w3.org/TR/html4/strict.dtd\">\n")
            msa_colored_html.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\"\n")
            msa_colored_html.write("\"http://www.w3.org/TR/html4/loose.dtd\">\n")
            msa_colored_html.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Frameset//EN\"\n")
            msa_colored_html.write("\"http://www.w3.org/TR/html4/frameset.dtd\">\n")
            msa_colored_html.write("<head>\n")
            msa_colored_html.write("<meta http-equiv=\"X-UA-Compatible\" content=\"IE=EmulateIE7\"/>\n")
            msa_colored_html.write(f"<link rel=\"stylesheet\" type=\"text/css\" href=\"{MSA_Score_CSS}\"/>\n")
            msa_colored_html.write("</head>\n")
            msa_colored_html.write("<H1 align=center><u>MSA color-coded by GUIDANCE scores</u></H1>\n\n")
            msa_colored_html.write("<table>\n")

            # get Align max_seq_length
            seq = alignment.alignment.sequences[1]
            align_width = alignment.get_alignment_length()
            # Align_width = len(seq.seq)

            # Print upper Scale
            msa_colored_html.write(
                "<tr>\n<td class=\"Score5\">&nbsp</td><td class=\"Seq_Name\">&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp</td><td>1</td>\n")
            i = 2
            while i < align_width:
                if (i % 10) == 0:
                    digits = [char for char in str(i)]
                    msa_colored_html.write("\n")
                    for digit in digits:
                        msa_colored_html.write(f"<td>{digit}</td>")
                        i += 1
                else:
                    msa_colored_html.write("<td></td>")
                    i += 1

            msa_colored_html.write("\n</tr>\n")

            for block_start in range(0, align_width, sequence_length_for_display):
                block_end = min(block_start + sequence_length_for_display, align_width)

                # Iterate over sequences and print up to sequenceLengthForDisplay residues
                for i in range(1, int(msa_depth) + 1):
                    seq = alignment[i - 1]

                    msa_colored_html.write("<tr>\n")
                    if seq_scores_data[str(i)] == "nan":
                        color_class = "ScoreNaN"
                    else:
                        color_class = colorstep[int(9 * float(seq_scores_data[str(i)]))]
                    # color_class = get_color_class(seq_scores_data.get(str(i), "nan"))
                    msa_colored_html.write(f"<td class=\"{color_class}\">&nbsp</td>\n")

                    if codes_file == "" or codes_file == None:
                        msa_colored_html.write(f"<td class=\"Seq_Name\">{seq.id}</td>\n")
                    else:
                        msa_colored_html.write(f"<td class=\"Seq_Name\">{code_names[seq.id][:25]}</td>\n")

                    # seq_block = seq.seq[block_start:block_end]
                    seq_block = seq[block_start:block_end]

                    for pos in range(len(seq_block)):
                        res = seq_block[pos]
                        if res == '-':
                            msa_colored_html.write("<td>{}</td>\n".format(res))
                        else:
                            color_class = ""
                            # if scores_data[i][block_start + pos] == "nan":
                            check = scores_data[str(i)][str(pos + 1)]
                            if scores_data[str(i)][str(pos + 1)] == 'nan':
                                color_class = "ScoreNaN"
                            else:
                                color_class = colorstep[int(9 * float(scores_data[str(i)][str(pos + 1)]))]
                            # color_class = get_color_class(scores_data[i][block_start + pos])
                            msa_colored_html.write("<td class=\"{}\">{}</td>\n".format(color_class, res))

                    msa_colored_html.write("</tr>\n\n")

                msa_colored_html.write("<tr></tr>\n\n")

            # Print lower Scale
            print_scale(msa_colored_html, align_width)

            # msa_colored_html.write("</body>\n</table>\n")

            # create_html_graph(col_scores_csv, out_html_file, x_label)
            create_html_graph(col_scores_csv, msa_colored_html, x_label)

            msa_colored_html.write("\n<br><b><u>Legend:</u><br><br>\nThe alignment confidence scale:</b><br>\n")
            print_color_scale(msa_colored_html, fontsize, colorstep_code)

            msa_colored_html.write("</body>\n</table>\n")

    except Exception as e:
        sys.exit(f"Can't open colored html file {out_html_file}: {e}\n")
    return ["OK"]


def print_scale(file, align_width):
    file.write("<tr>\n<td class=\"Score5\">&nbsp</td><td class=\"Seq_Name\"></td><td>1</td>\n")
    i = 2
    while i < align_width:
        if i % 10 == 0:
            digits = list(str(i))
            for digit in digits:
                file.write("\n<td>{}</td>".format(digit))
                i += 1
        else:
            file.write("<td></td>")
            i += 1
    file.write("\n</tr>\n")

@timeit
def create_html_graph(CSV_File, OUT, X_LABLE):
    # Create an HTML BARs graph
    # GET: 1. CSV FILE (The X var is the first Col and the Y is the second, X must be sorted)
    #      2. HTML OUTPUT
    #      3. X Lable (OPTIONAL)
    ############################################################################################
    try:
        # with open(Out, "a") as OUT, open(CSV_File) as DATA:
        with open(CSV_File) as DATA:
            # GRAPH PROPERTIES
            graphHeight = 200  # target height of graph
            BarWidth = 15  # width of bars
            maxResult = 1
            scale = 1
            last_x = 0

            # Set the scale
            scale = graphHeight / maxResult

            # SPACER BEFORE THE GRAPH
            OUT.write(
                "<tr><td class=\"Score5\">&nbsp</td><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr>\n")
            OUT.write(
                "<tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr>\n")
            OUT.write(
                "<tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr><tr><td class=\"Seq_Name\"></td></tr>\n")

            OUT.write("<tr>\n")
            OUT.write(
                "<td class=\"Score5\">&nbsp</td><td class=\"Seq_Name\" style = 'text-align: right'>{0}<br>SCORE</td>\n".format(
                    X_LABLE))

            header = DATA.readline()

            for line in DATA:
                line = line.strip()
                data = line.split(",")

                if data[1] != "NaN" and (int(data[0]) - last_x) == 1:
                    OUT.write("<td valign = bottom style = 'border-bottom: 1px solid black")
                    if data[0] == 1:  # The First Point, plot also the Y bar
                        OUT.write(";border-left: 1px solid black")
                    OUT.write(
                        ";'><img title=\"{}:{}\" src=\"https://taux.evolseq.net/guidance/static/images/blue.gif\" width=\"{}\" height=\"{}\" border=\"1\"></td>\n".format(
                            data[0], data[1], BarWidth, float(data[1]) * scale))
                    last_x = int(data[0])
                elif data[1] != "NaN":
                    while (int(data[0]) - last_x) != 1:
                        OUT.write(
                            "<td valign = bottom style = 'border-bottom: 1px solid black;'><img src=\"https://taux.evolseq.net/guidance/static/images/blue.gif\" width=\"{}\" height=\"0\" border=\"0\"></td>\n".format(
                                BarWidth))
                        last_x += 1
                    OUT.write(
                        "<td valign = bottom style = 'border-bottom: 1px solid black;'><img title=\"{}:{}\" src=\"https://taux.evolseq.net/guidance/static/images/blue.gif\" width=\"{}\" height=\"{}\" border=\"1\"></td>\n".format(
                            data[0], data[1], BarWidth, float(data[1]) * scale))
                    last_x = int(data[0])
                elif data[1] == "NaN":
                    OUT.write(
                        "<td valign = bottom style = 'border-bottom: 1px solid black;'><img title=\"{}:{}\" src=\"https://taux.evolseq.net/guidance/static/images/blue.gif\" width=\"{}\" height=\"{}\" border=\"1\"></td>\n".format(
                            data[0], data[1], BarWidth, float(data[1]) * scale))
                    last_x = int(data[0])

            OUT.write("</tr>\n")
            OUT.write("</table>\n")
            OUT.write("<b><P align=\"center\">Column</p></b>\n")

    except Exception as e:
        return f"Guidance::CreateHTML_Graph: Can't open file: {e}"


def print_color_scale(file, fontsize, colorstep_code):
    file.write(
        "<table style='table-layout: auto;margin-left: 0em;margin-right: 0em;padding:1px 1px 1px 1px; margin:1px 1px 1px 1px; border-collapse: collapse;' border=0 cols=1 width=310>\n<tr><td align=center>\n<font face='Courier New' color='black' size=+1><center>\n")
    for i in range(8, -1, -1):
        if i == 0:
            file.write(
                f"<font face='Courier New' color='white' size={fontsize}><span style='background: {colorstep_code[i]};'>&nbsp;{i + 1}&nbsp;</span></font>")
        else:
            file.write(
                f"<font face='Courier New' color='black' size={fontsize}><span style='background: {colorstep_code[i]};'>&nbsp;{i + 1}&nbsp;</span></font>")
    file.write(
        "</font></center>\n<center><table style = 'table-layout: auto;margin-left:0em;margin-right: 0em;padding:1px 1px 1px 1px; margin:1px 1px 1px 1px; border-collapse: collapse;' border=0 cols=3 width=310>\n<tr>\n<td align=left><b>Confident</b></td>\n<td align=center><b><---></b></td>\n<td align=right><b>Uncertain</b></td>\n</tr>\n</table></center>\n</td>\n</tr>\n</table>\n")
    file.write("<left><table style = 'table-layout: auto;margin-left: 0em;margin-right:0em;padding:1px 1px 1px 1px; margin:1px 1px 1px 1px; border-collapse:collapse;' border=0 cols=2 width=100>\n<tr>\n<td align=center class=\"ScoreNaN\">&nbsp;</td><td align=left>Insufficient Data</b></td>")

@timeit
def remove_low_sp_seq(msa_file, seq_sp_file, out_file, cutoff, removed_seq_file, seq_type="ByRowNum"):
    # Alow removal of sequences by their row number ($type=ByRowNum)in the MSA (when using the MSA set score raw file the scores are for MSA row)
    # or by sequence name ($type=BySeqName)
    # when using ByRowNum - we use BioPerl object [default]
    # when using BySeqName - we use hash to represent MSA (currently without extensive MSA validation)
    # print(msa_file,"\n", seq_sp_file,"\n", out_file,"\n", cutoff,"\n", removed_seq_file,"\n", seq_type)
    msa_hash_ref = None
    msa_length = None
    msa_order_array_ref = None
    in_handle = None
    aln = None

    if seq_type.upper() == "BYSEQNAME":
        ans = read_msa_to_hash(msa_file)
        if ans[0] != "OK":
            return f"remove_low_sp_seq: {ans[0]}\n"
        msa_hash_ref, msa_length, msa_order_array_ref = ans[1:]

    elif seq_type.upper() == "BYROWNUM":
        aln = AlignIO.read(msa_file, "fasta")
        # aln = in_handle[0]
        aln.verbose = True  # HAIM COMMENT
        set_displayname_flat(aln)  # HAIM COMMENT

    with open(seq_sp_file, "r") as seq_sp_scores, open(out_file, "w") as out, open(removed_seq_file,
                                                                                   "w") as removed_seq:
        seq_sp_scores.readline()  # Skip header
        for line in seq_sp_scores:
            # line = line.strip()
            match = re.match(r"^\s*(\d+)\s+(\d+\.\d+)", line)
            if match:
                seq_id, seq_sp_score = match.groups()
                seq, seq_name = None, None

                if seq_type.upper() == "BYSEQNAME":
                    seq = "".join(msa_hash_ref[seq_id-1])
                    seq_name = seq_id

                elif seq_type.upper() == "BYROWNUM":
                    seq_obj = aln[int(seq_id)-1]
                    seq = str(seq_obj.seq)
                    seq_name = seq_obj.id

                seq = seq.replace("-", "")
                seq = re.sub(r"(.{60,60})", r"\1\n", seq)
                seq += "\n" if not seq.endswith("\n") else ""

                if float(seq_sp_score) >= float(cutoff):
                    out.write(f">{seq_name}\n{seq}\n")
                else:
                    removed_seq.write(f">{seq_name}\n{seq}\n")

    return ["OK"]

@timeit
def remove_low_sp_seq_consider_z_score(msa_file, seq_sp_file, out_file, sp_cutoff, z_cutoff, removed_seq_file, seq_type="ByRowNum"):
    # Will remove all sequences in which their Z score is below cutoff and their SP score is also below cutoff.
    # Allow removal of sequences by their row number ($type=ByRowNum) in the MSA
    # or by sequence name ($type=BySeqName).

    mean, std = calculate_mean_and_std(seq_sp_file, 2)
    msa_hash_ref, msa_length, msa_order_array_ref = read_msa_to_hash(msa_file) if seq_type.upper() == "BYSEQNAME" else (None, None, None)

    alignment = AlignIO.read(msa_file, "fasta") if seq_type.upper() == "BYROWNUM" else None

    with open(seq_sp_file, "r") as seq_sp_scores, open(out_file, "w") as out, open(removed_seq_file, "w") as removed_seq:
        seq_sp_scores.readline()  # Skip header
        for line in seq_sp_scores:
            match = re.match(r"^\s*(\d+)\s+(\d+\.\d+)", line)
            # match = re.match(r"^\s*(\S+)\s+(\d+(\.\d+)?)", line)
            if match:
                seq_id, seq_sp_score = match.groups()
                if seq_type.upper() == "BYSEQNAME":
                    seq = "".join(msa_hash_ref[int(seq_id)-1])
                    seq_name = seq_id
                elif seq_type.upper() == "BYROWNUM":
                    seq_obj = alignment[int(seq_id) - 1]
                    seq = str(seq_obj.seq)
                    # seq = str(alignment[int(seq_id)-1].seq)
                    seq_name = seq_obj.id

                seq = seq.replace("-", "")
                seq = re.sub(r"(.{60,60})", r"\1\n", seq)
                seq += "\n" if not seq.endswith("\n") else ""
                z_score = "NaN"
                if std > 0:
                    z_score = (float(seq_sp_score) - mean) / std
                if seq_sp_score.lower() == "nan":
                    out.write(f">{seq_name}_SP_{seq_sp_score}_Z_{z_score}\n{seq}\n")
                elif z_score != "NaN" and z_score < z_cutoff:  # a negative outlier
                    if float(seq_sp_score) < float(sp_cutoff):  # To avoid filtering highly scored positions
                        removed_seq.write(f">{seq_name}_SP_{seq_sp_score}_Z_{z_score}\n{seq}\n")
                    else:
                        out.write(f">{seq_name}_SP_{seq_sp_score}_Z_{z_score}\n{seq}\n")
                else:
                    out.write(f">{seq_name}_SP_{seq_sp_score}_Z_{z_score}\n{seq}\n")
    return ["OK"]

@timeit
def make_JalView_output(JalView_Applet_Page, WorkingDir, http, inMsa, inMsa_With_names, scores, outJalviewFeaturesFile, NamesCodeFile, Jalview_AnnotFile, Data_File, Y_label):
    try:
        with open(JalView_Applet_Page, "w") as jalview_file:
            make_Jalview_Color_MSA(os.path.join(WorkingDir, inMsa), scores, os.path.join(WorkingDir, outJalviewFeaturesFile), NamesCodeFile)
            make_Jalview_AnnotationGraph(os.path.join(WorkingDir, Jalview_AnnotFile), Data_File, Y_label)

            jalview_file.write("<HTML>\n")
            jalview_file.write("<applet CODEBASE=\"http://guidance.tau.ac.il/\" CODE=\"jalview.bin.JalviewLite\" width=100% height=100% ARCHIVE=\"jalviewApplet.jar\">\n")
            jalview_file.write(f"<param name=\"file\" value=\"{http}{inMsa_With_names}\">\n")
            jalview_file.write(f"<param name=\"features\" value=\"{http}{outJalviewFeaturesFile}\">\n")
            jalview_file.write(f"<param name=\"annotations\" value=\"{http}{Jalview_AnnotFile}\">\n")
            jalview_file.write("<param name=\"application_url\" value=\"http://www.jalview.org/services/launchApp\">\n")
            jalview_file.write("<param name=\"showbutton\" VALUE=\"false\">\n")
            jalview_file.write("<param name=\"showConservation\" VALUE=\"false\">\n")
            jalview_file.write("<param name=\"showQuality\" VALUE=\"false\">\n")
            jalview_file.write("<param name=\"showConsensus\" VALUE=\"false\">\n")
            jalview_file.write("</APPLET>\n")
            jalview_file.write("</HTML>\n")

        return ["OK"]

    except Exception as e:
        return f"Error: {str(e)}"

@timeit
def make_Jalview_Color_MSA(inMsaFile, scoresFile, outJalviewFeaturesFile, codesFile=""):
    sequenceLengthForDisplay = 400000
    # Print HTML start
    with (open(outJalviewFeaturesFile, "w") as jalview_features):
        # Color steps
        jalview_features.write("Score1\t10C8D1\n")
        jalview_features.write("Score2\t8CFFFF\n")
        jalview_features.write("Score3\tD7FFFF\n")
        jalview_features.write("Score4\tEAFFFF\n")
        jalview_features.write("Score5\tFFFFFF\n")
        jalview_features.write("Score6\tFCEDF4\n")
        jalview_features.write("Score7\tFAC9DE\n")
        jalview_features.write("Score8\tF07DAB\n")
        jalview_features.write("Score9\tA02560\n")
        jalview_features.write("ScoreNaN\tC0C0C0\n")

        jalview_features.write("STARTGROUP\tGUIDANCE\n")

        with open(scoresFile, 'r') as scores_file:
            scores_data = {}
            code_names = {}
            for line in scores_file:
                if not line.startswith('#'):
                    line.strip()
                    line = trim(line)
                    col, seq, score = line.split()
                    if seq not in scores_data:
                        scores_data[seq] = {}
                    scores_data[seq][col] = score

            # Read Codes
            if codesFile != "" or codesFile != None:
                with open(codesFile, "r") as codes_file:
                    for line in codes_file:
                        # line = line.strip()
                        # seq_name, code = line.split('\t')
                        seq_name, code = line.strip().split('\t')
                        code_names[code] = seq_name

        # Read MSA
        try:
            alignment = AlignIO.read(inMsaFile, "fasta")
        except Exception as e:
            sys.exit(f"Can't open {inMsaFile}: {e}\n")
        alignment.verbose = True
        set_displayname_flat(alignment)
        ans = check_msa_licit_and_size(inMsaFile, "fasta", "no")
        if ans[0] == "OK":
            msa_depth = ans[1]
        else:
            answer_joined = " ".join(ans)
            return f"make_Jalview_Color_MSA:  {answer_joined}\n"
        # msa_depth = len(alignment)
        # the second Score 9 is the most confident (the score is exactly 1)
        color_step = ["Score1", "Score2", "Score3", "Score4", "Score5", "Score6", "Score7", "Score8", "Score9", "Score9"]
        colorstep_code = ["#10C8D1", "#8CFFFF", "#D7FFFF", "#EAFFFF", "#FFFFFF", "#FCEDF4", "#FAC9DE", "#F07DAB",
                          "#A02560"]
        align_width = alignment.get_alignment_length()

        # Print HTML start
        # msaColors = {}
        # msaPrintColors = {}
        # lineCounter = 0

        for blockStart in range(0, align_width, sequenceLengthForDisplay):
            blockEnd = min(blockStart + sequenceLengthForDisplay, align_width)

            for i in range(1, int(msa_depth) + 1):
                seq = alignment[i - 1]
                seq_block = seq[blockStart : blockEnd]

                # NEW
                color_class = ""
                # Print seq
                # seq_str = str(seq.seq[blockStart - 1:blockEnd])
                gaps = 0

                # for pos, res in enumerate(seq_block):
                for pos in range(len(seq_block)):
                    res = seq_block[int(pos)]
                    # color_class = ""
                    prob = ""
                    check = scores_data[str(i)][str(pos + 1)]

                    if res == '-':
                        gaps += 1
                        continue
                    else:
                        color_class = ""
                    # elif scores_data[str(i)][pos] != "nan":
                        if scores_data[str(i)][str(pos + 1)] != "nan":
                            color_class = color_step[int(9 * float(scores_data[str(i)][str(pos + 1)]))]
                    # get_color_class(float(scores_data[str(i)][pos]))
                            prob = scores_data[str(i)][str(pos + 1)]

                            if color_class != "Score5":
                                jalview_features.write(f"{prob}\tID_NOT_SPECIFIED\t{i - 1}\t{pos + 1 - gaps}\t{pos + 1 - gaps}\t{color_class}\t{prob}\n")

                        elif scores_data[str(i)][str(pos + 1)] == "nan":
                            color_class = "ScoreNaN"
                            jalview_features.write(f"NA\tID_NOT_SPECIFIED\t{i - 1}\t{pos + 1 - gaps}\t{pos + 1 - gaps}\t{color_class}\t{prob}\n")

        jalview_features.write("ENDGROUP\tGUIDANCE\n")

    # return ["OK"]


# def get_color_class(score):
#     color_step = ["Score1", "Score2", "Score3", "Score4", "Score5", "Score6", "Score7", "Score8", "Score9"]
#     return color_step[int(9 * score)]

@timeit
def make_Jalview_AnnotationGraph(Jalview_AnnotFile, Data_File, Y_label, Y_data_Col=1):
    last_x = 0
    with open(Jalview_AnnotFile, "w") as out_file:
        out_file.write("JALVIEW_ANNOTATION\n")
        out_file.write(f"BAR_GRAPH\t{Y_label}\t")

        with open(Data_File) as data_file:
            header = data_file.readline()  # skip header
            for line in data_file:
                line = line.strip()
                data = line.split(',')
                if data[Y_data_Col] != "nan" and (int(data[0]) - int(last_x)) == 1:
                    out_file.write(f"{data[Y_data_Col]},{data[Y_data_Col]},{data[Y_data_Col]}|")
                    last_x = int(data[0])
                elif data[1] != "nan":
                    while int(data[0]) - int(last_x) != 1:
                        out_file.write("0,0,0|")
                        last_x += 1
                    out_file.write(f"{data[Y_data_Col]},{data[Y_data_Col]},{data[Y_data_Col]}|")
                    last_x = int(data[0])
                elif data[Y_data_Col] == "nan":
                    out_file.write(f"{data[Y_data_Col]},{data[Y_data_Col]},{data[Y_data_Col]}|")
                    last_x = int(data[0])

    return ["OK"]

@timeit
def calculate_sp_scores(args_library):
    # if args_library.isServer == 1:
    #     if args_library.PROGRAM == "GUIDANCE":
    #         # print_message_to_output("Calculating GUIDANCE scores", args_library)
    #         message = "Calculating GUIDANCE scores"
    #     elif args_library.PROGRAM == "HoT":
    #         # print_message_to_output("Calculating HoT scores", args_library)
    #         message = "Calculating HoT scores"
    #     elif args_library.PROGRAM == "GUIDANCE2":
    #         # print_message_to_output("Calculating GUIDANCE2 scores", args_library)
    #         message = "Calculating GUIDANCE2 scores"
    #     elif args_library.PROGRAM == "GUIDANCE3":
    #         # print_message_to_output("Calculating GUIDANCE3 scores", args_library)
    #         message = "Calculating GUIDANCE3 scores"
    #
    #     print_message_to_output(message, args_library)

    if args_library.PROGRAM in ["GUIDANCE", "HoT"]:
        args_library.Output_Prefix = f"{args_library.dataset}.{args_library.MSA_Program}.Guidance"
    elif args_library.PROGRAM == "GUIDANCE2":
        args_library.Output_Prefix = f"{args_library.dataset}.{args_library.MSA_Program}.Guidance2"

    cmd = ""
    if args_library.userMSA_File != "" and args_library.Seq_Type == "Codons":       # user gave codon alignment
        args_library.Alignment_File_translated_from_user_codon_alignmet = f"{args_library.Alignment_File}.TranslatedProt"
        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(
                f"Codon_Aln_to_AA_Aln({args_library.WorkingDir}{args_library.Alignment_File},{args_library.WorkingDir}{args_library.Alignment_File_translated_from_user_codon_alignmet},{args_library.CodonTable},XCodonsFromALN.html)\n")
        ans = codon_alignment_to_aminoacids_alignment(
            f"{args_library.WorkingDir}{args_library.Alignment_File}",
            f"{args_library.WorkingDir}{args_library.Alignment_File_translated_from_user_codon_alignmet}",
            args_library.CodonTable,
            "XCodonsFromALN.html", args_library
        )
        if ans[0] != "OK":
            exit_on_error("user_error", ans, args_library)  # error
        elif ans[1] != "":  # warning
            if args_library.is_server == 1:  # server
                with open(args_library.output_page, "a") as output_file:
                    output_file.write(
                        f"<br><b><font color=\"red\" size4=>Warning:</b></font><font size=\"4\"> {ans[1]}</font>\n")
                log_file.write(f"Warning: {ans[1]}\n")
            else:
                print(f"Warning: {ans[1]}\n")
                with open(args_library.OutLogFile, "a") as log_file:
                    log_file.write(f"Warning: {ans[1]}\n")
        if os.path.getsize(
                os.path.join(f"{args_library.WorkingDir}",
                             f"{args_library.Alignment_File_translated_from_user_codon_alignmet}")) == 0 or not os.path.exists(
            os.path.join(
                f"{args_library.WorkingDir}", f"{args_library.Alignment_File_translated_from_user_codon_alignmet}")):
            exit_on_error("sys_error",
                          f"{args_library.WorkingDir}{args_library.Alignment_File_translated_from_user_codon_alignmet} does not exist/empty\n", args_library)
        # cmd = f"{args_library.msa_set_score_prog}   {os.path.join(args_library.WorkingDir, args_library.Alignment_File_translated_from_user_codon_alignmet)}   {os.path.join(args_library.WorkingDir, args_library.Output_Prefix)}   -d {args_library.Scoring_Alignments_Dir}  >  {args_library.WorkingDir}{args_library.dataset}.{args_library.MSA_Program}.msa_set_score.std"
        cmd = f"{args_library.msa_set_score_prog}   {os.path.join(args_library.WorkingDir, args_library.Alignment_File_translated_from_user_codon_alignmet)}   {os.path.join(args_library.WorkingDir, args_library.Output_Prefix)}   -d {args_library.Scoring_Alignments_Dir}  >  {args_library.WorkingDir}{args_library.dataset}.{args_library.MSA_Program}.msa_set_score.std"
    else:
        # cmd = f"{args_library.msa_set_score_prog}   {os.path.join(args_library.WorkingDir, args_library.Alignment_File)}   {os.path.join(args_library.WorkingDir, args_library.Output_Prefix)}   -d {args_library.Scoring_Alignments_Dir}  >  {args_library.WorkingDir}{args_library.dataset}.{args_library.MSA_Program}.msa_set_score.std"
        cmd = f"{args_library.msa_set_score_prog}   {os.path.join(args_library.WorkingDir, args_library.Alignment_File)}   {os.path.join(args_library.WorkingDir, args_library.Output_Prefix)}   -d {args_library.Scoring_Alignments_Dir}"
        # cmd = f"{args_library.msa_set_score_prog}   {os.path.join(args_library.WorkingDir, args_library.Alignment_File)}   {os.path.join(args_library.WorkingDir, args_library.Output_Prefix + f'_tree_{countTrees}')}   -d {args_library.Scoring_Alignments_Dir}"

    with open(args_library.OutLogFile, "a") as log_file:
        log_file.write(f"calculating SP scores: {cmd}\n")
        print(f"calculating SP scores: {cmd}\n")

    if os.path.exists(f"{args_library.Scoring_Alignments_Dir}/.DS_Store"):
        os.remove(f"{args_library.Scoring_Alignments_Dir}/.DS_Store")

    subprocess.call(cmd, shell=True)

    if not os.path.exists(f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_res.scr") or os.path.getsize(
            f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_res.scr") == 0:
        for i in range(3):
            try:
                subprocess.call(cmd, shell=True)
            except Exception as e:
                with open(args_library.OutLogFile, "a") as log_file:
                    log_file.write(f"Failed to calculate final scores {e}\n")
                    print(f"Failed to calculate final scores {e}\n")
                continue
            break

    if not os.path.exists(
                    f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_res.scr") or os.path.getsize(
                    f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_res.scr") == 0:
            exit_on_error("sys_error",
                      f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_res.scr does not exist/empty\n",
                      args_library)


    if args_library.PROGRAM == "HoT":
        with open(f"{os.path.join(args_library.WorkingDir, args_library.Alignment_File)}", "r") as orig_align, open(
                f"{os.path.join(args_library.WorkingDir, args_library.Alignment_File)}.NEW", "w") as new_align:
            for line in orig_align:
                if ">seq" in line:
                    if re.search(r">seq[0]+([1-9]+[0-9]*)$", line):
                        new_align.write(f">{int(re.search(r'>seq[0]+([1-9]+[0-9]*)$', line).group(1)) + 1}\n")
                    else:
                        new_align.write(">1\n")
                else:
                    new_align.write(line)
        os.rename(f"{os.path.join(args_library.WorkingDir, args_library.Alignment_File)}",
                  f"{args_library.WorkingDir}{args_library.Alignment_File}.ORIG")
        os.rename(f"{os.path.join(args_library.WorkingDir, args_library.Alignment_File)}.NEW",
                  f"{args_library.WorkingDir}{args_library.Alignment_File}")

    if args_library.isServer == 1:
        update_progress(f"{args_library.WorkingDir}{args_library.progress_report}", f"Finished Calculating {args_library.PROGRAM} scores")

def round_scores_file(score_file):
    with open(score_file, 'r') as file:
        lines = file.readlines()

    with open(score_file, 'w') as file:
        for line in lines:
            line = trim(line)
            parts = line.split()
            if len(parts) > 1 and parts[1].replace('.', '').isdigit():
                file.write(f"{parts[0]}\t{float(parts[1]):.3f}\n")
            elif len(parts) == 1:
                file.write(f"{parts[0]}\n")
            else:
                file.write(f"{parts[0]}\t{parts[1]}\n")

@timeit
def modify_score_files_for_codons_and_server(args_library):
    if args_library.Seq_Type == 'Codons':
        # We should modify the Scores files to be for CODONS - i.e each col score is repeated 3 times for the col and col+1,col+2

        # the alignment file should be back to CODONS alignment
        if args_library.userMSA_File == '':
            shutil.copy(os.path.join(args_library.WorkingDir, args_library.Alignment_File),
                        os.path.join(args_library.WorkingDir, args_library.Alignment_File_PROT))
            with open(args_library.OutLogFile, "a") as log_file:
                log_file.write("cp {} {}\n".format(os.path.join(args_library.WorkingDir, args_library.Alignment_File),
                                                   os.path.join(args_library.WorkingDir, args_library.Alignment_File_PROT)))
            # copy the hash - AA_to_DNA_aligned change the hash, so we copy it for later usage
            DNA_AA_cp = {key: value[:] for key, value in args_library.DNA_AA.items()}
            with open(args_library.OutLogFile, "a") as log_file:
                log_file.write("AA_to_DNA_aligned({}, {}, {})\n".format(
                    os.path.join(args_library.WorkingDir, args_library.Alignment_File_PROT),
                    os.path.join(args_library.WorkingDir, args_library.Alignment_File),
                    DNA_AA_cp))
            ans = AA_to_DNA_aligned(
                os.path.join(args_library.WorkingDir, args_library.Alignment_File_PROT),
                os.path.join(args_library.WorkingDir, args_library.Alignment_File),
                DNA_AA_cp)
            if ans[0] != 'ok':
                if ans[1] == 'user':
                    exit_on_error('user_error', ' '.join(ans), args_library)
                elif ans[1] == 'sys':
                    exit_on_error('sys_error', ' '.join(ans), args_library)

        # CP the scores files calculated for PROT
        #############################################
        for file_suffix in ["_col_col", "_res_pair_col", "_res_pair_seq_pair", "_res_pair_res", "_res_pair"]:
            prot_file = os.path.join(args_library.WorkingDir, "{}.PROT.scr".format(args_library.Output_Prefix + file_suffix))
            new_file = os.path.join(args_library.WorkingDir, "{}.scr".format(args_library.Output_Prefix + file_suffix))
            # print("mv {} {}".format(prot_file, new_file))
            os.rename(prot_file, new_file)

        # Convert the scores files calculated for PRO to codons
        #########################################################
        for file_suffix in ["_res_pair_col", "_res_pair_res", "_res_pair"]:
            prot_file = os.path.join(args_library.WorkingDir, "{}.PROT.scr".format(args_library.Output_Prefix + file_suffix))
            new_file = os.path.join(args_library.WorkingDir, "{}.scr".format(args_library.Output_Prefix + file_suffix))
            with open(args_library.OutLogFile, "a") as log_file:
                log_file.write("Convert_to_Codons_Numbering({}, {})\n".format(prot_file, new_file))
            ans = convert_to_codons_numbering(prot_file, new_file)
            if ans[0] != 'OK':
                exit_on_error('sys_error', ' '.join(ans), args_library)

    if args_library.isServer == 1:
        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(f"Round_scores_file {os.path.join(args_library.WorkingDir,args_library.Output_Prefix)}_res_pair_seq.scr \n")
        round_scores_file(os.path.join(args_library.WorkingDir, "{}_res_pair_seq.scr".format(args_library.Output_Prefix)))
        round_scores_file(os.path.join(args_library.WorkingDir, "{}_res_pair_col.scr".format(args_library.Output_Prefix)))

@timeit
def remove_sites(args_library):
    # Return names to the Seq score file
    ############################################
    args_library.Seq_Scores = args_library.Output_Prefix + "_res_pair_seq.scr" + "_with_Names"

    log_message = f"Guidance::codes2name_scoresFile_NEW(\"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_seq.scr\",\"{args_library.WorkingDir}{args_library.code_fileName}\",\"{args_library.WorkingDir}{args_library.Alignment_File}\",\"{args_library.WorkingDir}{args_library.Seq_Scores}\");\n"
    with open(args_library.OutLogFile, "a") as log_file:
        log_file.write(f"{log_message}")
    ans = codes2name_scoresFile_NEW(f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_seq.scr",
                                    f"{args_library.WorkingDir}{args_library.code_fileName}",
                                    f"{args_library.WorkingDir}{args_library.Alignment_File}",
                                    f"{args_library.WorkingDir}{args_library.Seq_Scores}")
    if ans[0] != "OK":
        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write("Guidance::codes2name_scoresFile: " + ''.join(ans))

    # Remove sites with SP-score < Col sp_cutoff
    ############################################
    args_library.Alignment_File_without_low_SP_Col = f"{args_library.dataset}.{args_library.MSA_Program}.Without_low_SP_Col"
    args_library.removed_low_SP_SITE = f"{args_library.SeqsFile}.{args_library.dataset}.{args_library.MSA_Program}.Removed_Col"

    log_message = f"Guidance::removeLowSPsites_NoBioPerl (\"{args_library.WorkingDir}{args_library.Alignment_File}\",\"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr\",\"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Col}\",{args_library.SP_COL_CUTOFF},\"{args_library.WorkingDir}{args_library.removed_low_SP_SITE}\");\n"
    with open(args_library.OutLogFile, "a") as log_file:
        log_file.write(f"{log_message}\n")
    ans = remove_low_sp_sites_no_bioperl(f"{args_library.WorkingDir}{args_library.Alignment_File}",
                                         f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr",
                                         f"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Col}",
                                         args_library.SP_COL_CUTOFF,
                                         f"{args_library.WorkingDir}{args_library.removed_low_SP_SITE}")
    if ans[0] == "OK":
        args_library.REMOVED_SITES = ans[1]
        args_library.MSA_LENGTH = ans[2]
    with open(args_library.OutLogFile, "a") as log_file:
        log_file.write(f"REMOVED_SITES:{args_library.REMOVED_SITES}\n")
        log_file.write(f"MSA_LENGTH:{args_library.MSA_LENGTH}\n")

    # JS - print for flask version
    msa_length_file = os.path.join(args_library.WorkingDir, 'MSA_LENGTH')
    with open(msa_length_file, "w") as length_file:
        length_file.write(f"{args_library.MSA_LENGTH}\n")

    args_library.Alignment_File_without_low_SP_Col_with_Names = f"{args_library.Alignment_File_without_low_SP_Col}.With_Names"
    if os.path.getsize(f"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Col}") > 0:  # Not EMPTY
        ans = codes2name_fasta_from1(f"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Col}",
                                     f"{args_library.WorkingDir}{args_library.code_fileName}",
                                     f"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Col_with_Names}")
        if ans[0] != "OK":
            exit_on_error("sys_error",
                          f"Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(\"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Col}\",\"{args_library.WorkingDir}{args_library.code_fileName}\",\"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Col_with_Names}\") failed: {''.join(ans)}\n",
                          args_library)

    if args_library.Z_Col_Cutoff != "NA":
        # Remove sites with SP-score < Col sp_cutoff if Z is below cutoff
        ##################################################################
        args_library.Alignment_File_without_low_SP_Z_Col = f"{args_library.dataset}.{args_library.MSA_Program}.Without_low_SP_and_low_Z_Col"
        args_library.removed_low_SP_Z_SITE = f"{args_library.SeqsFile}.{args_library.dataset}.{args_library.MSA_Program}.Removed_SP_and_Z_Col"

        log_message = f"Guidance::removeLowSPsites_NoBioPerl_Consider_Z (\"{args_library.WorkingDir}{args_library.Alignment_File}\",\"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr\",\"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Z_Col}\",{args_library.SP_COL_CUTOFF},{args_library.Z_Col_Cutoff},\"{args_library.WorkingDir}{args_library.removed_low_SP_Z_SITE}\");\n"
        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(f"{log_message}\n")
        ans = remove_low_sp_sites_consider_z(f"{args_library.WorkingDir}{args_library.Alignment_File}",
                                             f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr",
                                             f"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Z_Col}",
                                             args_library.SP_COL_CUTOFF, args_library.Z_Col_Cutoff,
                                             f"{args_library.WorkingDir}{args_library.removed_low_SP_Z_SITE}")
        if ans[0] == "OK":
            args_library.REMOVED_Z_SITES = ans[1]
            args_library.MSA_LENGTH = ans[2]
        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(f"REMOVED_SITES_SP_AND_Z: {args_library.REMOVED_Z_SITES}\n")
            log_file.write(f"MSA_LENGTH: {args_library.MSA_LENGTH}\n")

        # JS - print for flask version
        msa_length_file = os.path.join(args_library.WorkingDir, 'MSA_LENGTH')
        with open(msa_length_file, "w") as length_file:
            length_file.write(f"{args_library.MSA_LENGTH}\n")

        args_library.Alignment_File_without_low_SP_Z_Col_with_Names = f"{args_library.Alignment_File_without_low_SP_Z_Col}.With_Names"
        alignment_file_path = f"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Z_Col}"
        if os.path.getsize(alignment_file_path) > 0:
            ans = codes2name_fasta_from1(f"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Z_Col}",
                                         f"{args_library.WorkingDir}{args_library.code_fileName}",
                                         f"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Z_Col_with_Names}")
            if ans[0] != "OK":
                joined_answer = "".join(ans)
                exit_on_error("sys_error",
                              f"Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(\"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Z_Col}\",\"{args_library.WorkingDir}{args_library.code_fileName}\",\"{args_library.WorkingDir}{args_library.Alignment_File_without_low_SP_Z_Col_with_Names}\") failed: {joined_answer}\n",
                              args_library)

@timeit
def prepare_plots(args_library):
    try:
        with open(f"{args_library.OutLogFile}", "a") as log_file:
            log_file.write(
                f"ans=Convert_to_CSV(\"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr\",\"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr.csv\");\n")
            ans = convert_to_csv(f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr",
                                 f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr.csv")
            joined_ans = " ".join(ans)
            log_file.write(f"ANS: {joined_ans}\n")
            args_library.res_pair_res_html_file = f"{args_library.dataset}.{args_library.MSA_Program}.Guidance_res_pair_res.html"
            log_file.write(
                f"Guidance::printColoredAlignment_With_CSS(\"{args_library.WorkingDir}{args_library.Alignment_File}\",\"{args_library.WorkingDir}{args_library.res_pair_res_html_file}\",\"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_res.scr\",\"{args_library.WorkingDir}{args_library.code_fileName}\",\"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr.csv\",{args_library.PROGRAM},\"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_seq.scr\");\n")
            ans = print_colored_alignment_with_css(f"{args_library.WorkingDir}{args_library.Alignment_File}",
                                                   f"{args_library.WorkingDir}{args_library.res_pair_res_html_file}",
                                                   f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_res.scr",
                                                   f"{args_library.WorkingDir}{args_library.code_fileName}",
                                                   f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr.csv",
                                                   args_library.PROGRAM,
                                                   f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_seq.scr")
            joined_ans = " ".join(ans)
            log_file.write(f"ANS: {joined_ans}\n")
    except Exception as e:
        sys.exit(f"ERROR: Could not open log file in prepare_plots(): {e}\n")


@timeit
def remove_sequences_sp_score(args_library):
    # remove seq with SP-score < Seq sp_cutoff
    ############################################
    args_library.Seq_File_without_low_SP_SEQ = args_library.SeqsFile + ".Without_low_SP_Seq"
    args_library.removed_low_SP_SEQ = args_library.SeqsFile + ".Removed_Seq"
    args_library.Seq_File_without_low_SP_SEQ_with_Names = args_library.Seq_File_without_low_SP_SEQ + ".With_Names"
    args_library.removed_low_SP_SEQ_With_Names = args_library.removed_low_SP_SEQ + ".With_Names"

    try:
        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(
                f"Guidance::removeLowSPseq(\"{args_library.WorkingDir}{args_library.Alignment_File}\",\"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_seq.scr\",\"{args_library.WorkingDir}{args_library.Seq_File_without_low_SP_SEQ}\",{args_library.SP_SEQ_CUTOFF},\"{args_library.WorkingDir}{args_library.removed_low_SP_SEQ}\")\n")
    except Exception as e:
        sys.exit(f"Error occurred while running add_original_seq_names_to_the_MSA(): {e}\n")

    remove_low_sp_seq(f"{args_library.WorkingDir}{args_library.Alignment_File}",
                      f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_seq.scr",
                      f"{args_library.WorkingDir}{args_library.Seq_File_without_low_SP_SEQ}", args_library.SP_SEQ_CUTOFF,
                      f"{args_library.WorkingDir}{args_library.removed_low_SP_SEQ}")
    if os.path.getsize(f"{args_library.WorkingDir}{args_library.Seq_File_without_low_SP_SEQ}") > 0:  # NOT EMPTY
        ans = codes2name_fasta_from1(f"{args_library.WorkingDir}{args_library.Seq_File_without_low_SP_SEQ}",
                                     f"{args_library.WorkingDir}{args_library.code_fileName}",
                                     f"{args_library.WorkingDir}{args_library.Seq_File_without_low_SP_SEQ_with_Names}")
        if ans[0] != "OK":
            exit_on_error("sys_error",
                          f"Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(\"{args_library.WorkingDir}{args_library.Seq_File_without_low_SP_SEQ}\",\"{args_library.WorkingDir}{args_library.code_fileName}\",\"{args_library.WorkingDir}{args_library.Seq_File_without_low_SP_SEQ_with_Names}\") failed:" +
                          ''.join(ans) + "\n", args_library)

    if os.path.getsize(f"{args_library.WorkingDir}{args_library.removed_low_SP_SEQ}") > 0:  # NOT EMPTY Seq were removed
        ans = codes2name_fasta_from1(f"{args_library.WorkingDir}{args_library.removed_low_SP_SEQ}",
                                     f"{args_library.WorkingDir}{args_library.code_fileName}",
                                     f"{args_library.WorkingDir}{args_library.removed_low_SP_SEQ_With_Names}")
        if ans[0] != "OK":
            exit_on_error("sys_error",
                          f"Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(Guidance::codes2nameFastaFrom1(\"{args_library.WorkingDir}{args_library.removed_low_SP_SEQ}\",\"{args_library.WorkingDir}{args_library.code_fileName}\",\"{args_library.WorkingDir}{args_library.removed_low_SP_SEQ_With_Names}\") failed:" + ''.join(
                              ans) + "\n", args_library)

@timeit
def remove_sequences_sp_and_z_score(args_library):
    # remove seq with SP-score < Seq sp_cutoff if Z<(-Z_Cutoff)
    ##############################################################
    if args_library.Z_Seq_Cutoff != "NA":
        seq_file_without_low_sp_z_seq = args_library.SeqsFile + ".Without_low_Z_and_low_SP_Seq"
        removed_low_sp_z_seq = args_library.SeqsFile + ".Removed_Z_SP_Seq"
        seq_file_without_low_sp_z_seq_with_names = seq_file_without_low_sp_z_seq + ".With_Names"
        removed_low_sp_z_seq_with_names = removed_low_sp_z_seq + ".With_Names"

        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(f"remove_low_sp_seq_consider_z_score(\"{args_library.WorkingDir}{args_library.Alignment_File}\","
                           f"\"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_seq.scr\","
                           f"\"{args_library.WorkingDir}{seq_file_without_low_sp_z_seq}\","
                           f"{args_library.SP_SEQ_CUTOFF},{args_library.Z_Seq_Cutoff},\"{args_library.WorkingDir}{removed_low_sp_z_seq}\");\n")
            remove_low_sp_seq_consider_z_score(
                f"{args_library.WorkingDir}{args_library.Alignment_File}",
                f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_seq.scr",
                f"{args_library.WorkingDir}{seq_file_without_low_sp_z_seq}",
                args_library.SP_SEQ_CUTOFF,
                float(args_library.Z_Seq_Cutoff),
                f"{args_library.WorkingDir}{removed_low_sp_z_seq}"
            )
        if os.path.getsize(f"{args_library.WorkingDir}{seq_file_without_low_sp_z_seq}") > 0:
            ans = codes2name_fasta_from1(f"{args_library.WorkingDir}{args_library.seq_file_without_low_sp_z_seq}",
                                         f"{args_library.WorkingDir}{args_library.code_fileName}",
                                         f"{args_library.WorkingDir}{args_library.seq_file_without_low_sp_z_seq_with_names}")
            if ans[0] != "OK":
                exit_on_error("sys_error",
                              f"Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(\"{args_library.WorkingDir}{args_library.seq_file_without_low_sp_z_seq}\",\"{args_library.WorkingDir}{args_library.code_fileName}\",\"{args_library.WorkingDir}{args_library.seq_file_without_low_sp_z_seq_with_names}\") failed:" + " ".join(
                                  ans) + "\n", args_library)

        if os.path.getsize(f"{args_library.WorkingDir}{args_library.removed_low_sp_z_seq}") > 0:
            ans = codes2name_fasta_from1(f"{args_library.WorkingDir}{args_library.removed_low_sp_z_seq}",
                                         f"{args_library.WorkingDir}{args_library.code_fileName}",
                                         f"{args_library.WorkingDir}{args_library.removed_low_sp_z_seq_with_names}")
            if ans[0] != "OK":
                exit_on_error("sys_error",
                              f"Guidance::codes2nameFastaFrom1: Guidance::codes2nameFastaFrom1(\"{args_library.WorkingDir}{args_library.removed_low_sp_z_seq}\",\"{args_library.WorkingDir}{args_library.code_fileName}\",\"{args_library.WorkingDir}{args_library.removed_low_sp_z_seq_with_names}\") failed:" + " ".join(
                                  ans) + "\n", args_library)


def make_jalview(args_library):
    args_library.JalView_page = args_library.Output_Prefix + "_JalView.html"
    args_library.JalView_Features = args_library.Output_Prefix + "_JalView_Features"
    args_library.JalView_Annotations = args_library.Output_Prefix + "_JalView_Annot"
    if args_library.isServer == 1:  # JalView Applet output (for web-server)
        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(
                f"Prepare JalView outputs: Guidance::make_JalView_output({args_library.WorkingDir}{args_library.JalView_page},{args_library.WorkingDir},{args_library.run_url},{args_library.Alignment_File},{args_library.Alignment_File_With_Names},{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_res.scr,{args_library.JalView_Features},{args_library.WorkingDir}{args_library.code_fileName},{args_library.JalView_Annotations},{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr.csv,{args_library.PROGRAM} scores);\n")
        make_JalView_output(f"{args_library.WorkingDir}{args_library.JalView_page}", f"{args_library.WorkingDir}",
                            f"{args_library.run_url}", f"{args_library.Alignment_File}",
                            f"{args_library.Alignment_File_With_Names}",
                            f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_res.scr",
                            f"{args_library.JalView_Features}", f"{args_library.WorkingDir}{args_library.code_fileName}",
                            f"{args_library.JalView_Annotations}",
                            f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr.csv",
                            f"{args_library.PROGRAM} scores")
    else:
        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(f"Guidance::make_Jalview_Color_MSA(\"{args_library.WorkingDir}{args_library.Alignment_File}\",\"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_res.scr\",\"{args_library.WorkingDir}{args_library.JalView_Features}\",\"{args_library.WorkingDir}{args_library.code_fileName}\");\n")
            make_Jalview_Color_MSA(f"{args_library.WorkingDir}{args_library.Alignment_File}", f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_res.scr", f"{args_library.WorkingDir}{args_library.JalView_Features}", f"{args_library.WorkingDir}{args_library.code_fileName}")
            log_file.write(f"make_Jalview_AnnotationGraph(\"{args_library.WorkingDir}{args_library.JalView_Annotations}\",\"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr.csv\",{args_library.PROGRAM}.\" scores\");\n")
            make_Jalview_AnnotationGraph(f"{args_library.WorkingDir}{args_library.JalView_Annotations}", f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_col.scr.csv", f"{args_library.PROGRAM} scores")

def print_output_to_the_server(args_library):
    pass

# def create_png_for_seqscores(prefix):
@timeit
def create_png_for_seqscores(args_library):
    data_file = f"{args_library.WorkingDir}{args_library.Output_Prefix}_res_pair_seq.scr_with_Names"
    # data_file = f"{prefix}_res_pair_seq.scr_with_Names"
    scores = []
    names = []
    with open(data_file, "r") as in_file:
        for line in in_file:
            line = line.strip()
            if line.startswith("SEQUENCE_NAME"):
                continue
            cols = line.split()
            if len(cols) == 2:
                name, score = cols[0], float(cols[1])
                if score != float('nan') and math.isnan(score) != True:
                    scores.append(score)
                    names.append(name)
    fig1 = plt.gcf()
    n, bins, patches = plt.hist(scores)
    plt.xlabel('Sequence score')
    plt.ylabel('Number of sequences')
    plt.title('Histogram with Sequence scores distribution')
    # plt.show()
    fig1.savefig(os.path.join(args_library.WorkingDir,'histogram_seq_scores_distribution.png'))
    plt.close(fig1)
    # fig1.savefig(f'{prefix}_histogram_seq_scores_distribution.png')

    # Create a boxplot
    plt.boxplot(scores)
    # Calculate the first and third quartiles
    q1 = np.percentile(scores, 25)
    q3 = np.percentile(scores, 75)
    # Calculate the interquartile range (IQR)
    iqr = q3 - q1
    # Define the outlier threshold (1.5 times the IQR)
    outlier_threshold = 1.5 * iqr

    # Identify the outliers
    outliers = []
    for i, value in enumerate(scores):
        if value > q3 + outlier_threshold or value < q1 - outlier_threshold:
            outliers.append((i, value))

    sorted_outliers = sorted(outliers, key=lambda x: x[1])
    index = 1
    # Get y-axis tick marks
    yticks = plt.yticks()[0]
    # Find ticks distance
    dist = yticks[1] - yticks[0]
    # Annotate outliers with their names
    for i, value in sorted_outliers:
        # if value < 0.94:
        #     plt.annotate(f'{names[i]}', xy=(float(f"1.0{index}"), value), xytext=(float(f"1.11{index}"), value + float(f"0.{index}")* dist),
        #                  arrowprops=dict(facecolor='red', shrink=0.05), fontsize = 'small')
        # else:
        #     plt.annotate(f'{names[i]}', xy=(float(f"1.0{index}"), value),
        #                  xytext=(float(f"1.11{index}"), value + float(f"0.{index}") * dist),
        #                  arrowprops=dict(facecolor='red', shrink=0.05), fontsize='small')
        plt.annotate(f'{names[i]}', xy=(float(f"1.0{index}"), value),
                     xytext=(float(f"1.11{index}"), value + float(f"0.{index}") * 1.7 * dist),
                     arrowprops=dict(facecolor='red', shrink=0.05), fontsize='small')
        index += 1

    fig2 = plt.gcf()
    # plt.xlabel('Data')
    plt.ylabel('Sequence score')
    plt.title('Boxplot with Sequence scores and outlier names')
    fig = plt.figure(figsize=(10, 7))
    plt.tight_layout()
    # plt.show()
    fig2.savefig(os.path.join(args_library.WorkingDir,'boxplot_seq_scores_and_outliers.png'))
    plt.close(fig2)
    # fig2.savefig(f'{prefix}_boxplot_seq_scores_and_outliers.png')

@timeit
def calculate_sp_scores_convergence(args_library, countTrees):
    # if args_library.isServer == 1:
    #     if args_library.PROGRAM == "GUIDANCE":
    #         # print_message_to_output(f"Calculating GUIDANCE scores for tree # {countTrees}", args_library)
    #         update_progress( f"{args_library.WorkingDir}{args_library.progress_report}",f"Calculating GUIDANCE scores for tree # {countTrees}")
    #     elif args_library.PROGRAM == "HoT":
    #         # print_message_to_output(f"Calculating HoT scores for tree # {countTrees}", args_library)
    #         update_progress(f"{args_library.WorkingDir}{args_library.progress_report}",
    #                         f"Calculating HoT scores for tree # {countTrees}")
    #     elif args_library.PROGRAM == "GUIDANCE2":
    #         # print_message_to_output(f"Calculating GUIDANCE2 scores for tree # {countTrees}", args_library)
    #         update_progress(f"{args_library.WorkingDir}{args_library.progress_report}",
    #                         f"Calculating GUIDANCE2 scores for tree # {countTrees}")
    #     elif args_library.PROGRAM == "GUIDANCE3":
    #         # print_message_to_output(f"Calculating GUIDANCE3 scores for tree # {countTrees}", args_library)
    #         update_progress(f"{args_library.WorkingDir}{args_library.progress_report}",
    #                         f"Calculating GUIDANCE3 scores for tree # {countTrees}")

    if args_library.PROGRAM in ["GUIDANCE", "HoT"]:
        args_library.Output_Prefix = f"{args_library.dataset}.{args_library.MSA_Program}.Guidance"
    elif args_library.PROGRAM == "GUIDANCE2":
        args_library.Output_Prefix = f"{args_library.dataset}.{args_library.MSA_Program}.Guidance2"

    cmd = ""
    if args_library.userMSA_File != "" and args_library.Seq_Type == "Codons":
        args_library.Alignment_File_translated_from_user_codon_alignmet = f"{args_library.Alignment_File}.TranslatedProt"
        with open(args_library.OutLogFile, "a") as log_file:
            log_file.write(
                f"Codon_Aln_to_AA_Aln({args_library.WorkingDir}{args_library.Alignment_File},{args_library.WorkingDir}{args_library.Alignment_File_translated_from_user_codon_alignmet},{args_library.CodonTable},XCodonsFromALN.html)\n")
        ans = codon_alignment_to_aminoacids_alignment(
            f"{args_library.WorkingDir}{args_library.Alignment_File}",
            f"{args_library.WorkingDir}{args_library.Alignment_File_translated_from_user_codon_alignmet}",
            args_library.CodonTable,
            "XCodonsFromALN.html", args_library
        )
        if ans[0] != "OK":
            exit_on_error("user_error", ans, args_library)  # error
        elif ans[1] != "":  # warning
            if args_library.is_server == 1:  # server
                with open(args_library.output_page, "a") as output_file:
                    output_file.write(
                        f"<br><b><font color=\"red\" size4=>Warning:</b></font><font size=\"4\"> {ans[1]}</font>\n")
                log_file.write(f"Warning: {ans[1]}\n")
            else:
                print(f"Warning: {ans[1]}\n")
                with open(args_library.OutLogFile, "a") as log_file:
                    log_file.write(f"Warning: {ans[1]}\n")
        if os.path.getsize(
                os.path.join(f"{args_library.WorkingDir}",
                             f"{args_library.Alignment_File_translated_from_user_codon_alignmet}")) == 0 or not os.path.exists(
            os.path.join(
                f"{args_library.WorkingDir}", f"{args_library.Alignment_File_translated_from_user_codon_alignmet}")):
            exit_on_error("sys_error",
                          f"{args_library.WorkingDir}{args_library.Alignment_File_translated_from_user_codon_alignmet} does not exist/empty\n", args_library)
        # cmd = f"{args_library.msa_set_score_prog}   {os.path.join(args_library.WorkingDir, args_library.Alignment_File_translated_from_user_codon_alignmet)}   {os.path.join(args_library.WorkingDir, args_library.Output_Prefix) + f'_tree_{countTrees}'}   -d {args_library.Scoring_Alignments_Dir}  >  {args_library.WorkingDir}{args_library.dataset}.{args_library.MSA_Program}.msa_set_score.std"
        cmd = f"{args_library.msa_set_score_prog}   {os.path.join(args_library.WorkingDir, args_library.Alignment_File_translated_from_user_codon_alignmet)}   {os.path.join(args_library.WorkingDir, args_library.Output_Prefix) + f'_tree_{countTrees}'}   -d {args_library.Scoring_Alignments_Dir}"
    else:
        cmd = f"{args_library.msa_set_score_prog}   {os.path.join(args_library.WorkingDir, args_library.Alignment_File)}   {os.path.join(args_library.WorkingDir, args_library.Output_Prefix + f'_tree_{countTrees}')}   -d {args_library.Scoring_Alignments_Dir}"
    with open(args_library.OutLogFile, "a") as log_file:
        log_file.write(f"calculating SP scores for tree # {countTrees}: {cmd}\n")
        print(f"calculating SP scores for tree # {countTrees}: {cmd}\n")
    if os.path.exists(f"{args_library.Scoring_Alignments_Dir}/.DS_Store"):
        os.remove(f"{args_library.Scoring_Alignments_Dir}/.DS_Store")
    subprocess.call(cmd, shell=True)
    alt_msas = len(os.listdir(args_library.Scoring_Alignments_Dir))
    if not os.path.exists(f"{args_library.WorkingDir}{args_library.Output_Prefix + f'_tree_{countTrees}'}_res_pair_res.scr") or os.path.getsize(
            f"{args_library.WorkingDir}{args_library.Output_Prefix + f'_tree_{countTrees}'}_res_pair_res.scr") == 0:
        exit_on_error("sys_error",
                      f"{args_library.WorkingDir}{args_library.Output_Prefix + f'_tree_{countTrees}'}_res_pair_res.scr does not exist/empty\n",
                      args_library)
    if args_library.PROGRAM == "HoT":
        with open(f"{os.path.join(args_library.WorkingDir, args_library.Alignment_File)}", "r") as orig_align, open(
                f"{os.path.join(args_library.WorkingDir, args_library.Alignment_File)}.NEW", "w") as new_align:
            for line in orig_align:
                if ">seq" in line:
                    if re.search(r">seq[0]+([1-9]+[0-9]*)$", line):
                        new_align.write(f">{int(re.search(r'>seq[0]+([1-9]+[0-9]*)$', line).group(1)) + 1}\n")
                    else:
                        new_align.write(">1\n")
                else:
                    new_align.write(line)
        os.rename(f"{os.path.join(args_library.WorkingDir, args_library.Alignment_File)}",
                  f"{args_library.WorkingDir}{args_library.Alignment_File}.ORIG")
        os.rename(f"{os.path.join(args_library.WorkingDir, args_library.Alignment_File)}.NEW",
                  f"{args_library.WorkingDir}{args_library.Alignment_File}")

    return alt_msas

@timeit
def add_scores_to_dict(args_library, epsilon, countTrees, lock):
    # score = 10*epsilon  #just random score not satisfying the condition of convergence
    MSA_score_file = os.path.join(args_library.WorkingDir, f"{args_library.Output_Prefix + f'_tree_{countTrees}'}_msa.scr")

    with open(MSA_score_file, 'r') as f:
        for line in f:
            if "#MEAN_RES_PAIR_SCORE" in line:
                with lock:
                    args_library.mean_res_pair_score.append(float(line.strip().split()[1]))
                    args_library.mean_col_score.append(float(line.strip().split()[3]))
    f.close()
    # print(args_library.mean_res_pair_score)
    # print(args_library.mean_col_score)
@timeit
def check_convergence(args_library, epsilon):
    score1, score2 = 10 * epsilon, 10 * epsilon
    if len(args_library.mean_col_score)>1 and len(args_library.mean_res_pair_score)>1:
        score1 = abs(args_library.mean_col_score[-1] - args_library.mean_col_score[-2])
        score2 = abs(args_library.mean_res_pair_score[-1] - args_library.mean_res_pair_score[-2])


    if score1 <= epsilon and score2 <= epsilon:
        return 1
    else:
        return 0
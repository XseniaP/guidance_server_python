from guidance_args_library import *
from guidance_scoring_and_visualization import *
from guidance_program_call_functions import *

Bin = os.path.dirname(sys.argv[0])
BIN_DIR = os.path.dirname(Bin)
RESULTS = os.path.join(BIN_DIR, "results/Guidance")
# --seqFile  /Users/kpolonsky/PycharmProjects/HoT_Py/Seqs.Orig.fas --msaProgram MAFFT --seqType aa --outDir /Users/kpolonsky/PycharmProjects/HoT_Py/ENSG00000017260_1/ --msaFile /Users/kpolonsky/PycharmProjects/HoT_Py/MSA.MAFFT.aln --program HoT


if __name__ == "__main__":

    args_library = Library()
    args_library.check_and_set_input_and_output_variables(sys.argv)

    # NOW WE ALWAYS WITH AA SEQ
    if args_library.isServer == 1:
        with open(f"{args_library.WorkingDir}{args_library.output_page}", "a") as output_file:
            output_file.write("<h4><font face=Verdana><u>Running Messages:</u></h4></font>\n")

    # if is_server == 1:
    #     change_qued_to_running(VARS['WorkingDir'] + VARS['output_page'])
    #     with open(VARS['WorkingDir'] + VARS['output_page'], "a") as OUTPUT:
    #         OUTPUT.write("<h4><font face=Verdana><u>Running Progress:</u></h4></font>\n")
            print_initial_running_progress(args_library)
            output_file.write("<div id='includedContent'>{{ progress_report|safe }}</div>\n")


    if args_library.PROGRAM in ["GUIDANCE", "GUIDANCE3"]:
        run_guidance(args_library)

    elif args_library.PROGRAM == "HoT":
        run_hot(args_library)

    elif args_library.PROGRAM in ["GUIDANCE2", "GUIDANCE3_HOT"]:
        run_guidance2(args_library)

    calculate_sp_scores(args_library)
    modify_score_files_for_codons_and_server(args_library)
    remove_sites(args_library)
    prepare_plots(args_library)
    add_original_seq_names_to_the_MSA(args_library)
    remove_sequences_sp_score(args_library)
    remove_sequences_sp_and_z_score(args_library)
    make_jalview(args_library)
    flag_that_finished_ok(args_library)
    create_tar_archives(args_library)
    create_png_for_seqscores(args_library)

    if os.path.exists(os.path.join(RESULTS, "timestamps.txt")):
        shutil.move(os.path.join(RESULTS, "timestamps.txt"), os.path.join(args_library.WorkingDir, "timestamps.txt"))


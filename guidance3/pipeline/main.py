from guidance3.config import RunConfig
from guidance3.pipeline.scoring import (
    calculate_sp_scores,
    modify_score_files_for_codons_and_server,
    remove_sites,
    prepare_plots,
    remove_sequences_sp_score,
    remove_sequences_sp_and_z_score,
    make_jalview,
)
from guidance3.pipeline.alignment import (
    run_guidance,
    run_hot,
    run_guidance3,
    prepare_rerun_parameters,
)
from guidance3.sequences.fasta import (
    add_original_seq_names_to_the_MSA,
    create_tar_archives,
    select_best_msa,
)
from guidance3.utils.common import (
    flag_that_finished_ok,
    print_initial_running_progress,
    send_finish_email_to_user,
)
import os
import sys
import shutil

Bin = os.path.dirname(sys.argv[0])
BIN_DIR = os.path.dirname(Bin)
RESULTS = os.path.join(BIN_DIR, "results/Guidance")


def main(argv=None):
    if argv is None:
        argv = sys.argv

    config = RunConfig()
    config.check_and_set_input_and_output_variables(argv)

    try:
        _run(config)
    except Exception as e:
        import traceback
        msg = traceback.format_exc()
        sys.stderr.write(msg)
        if getattr(config, 'WorkingDir', None):
            errors_path = os.path.join(config.WorkingDir, 'errors.txt')
            with open(errors_path, 'w') as ef:
                ef.write(f"Unhandled error:\n{msg}")
        sys.exit(1)


def _run(config):
    if config.isServer == 1:
        with open(f"{config.WorkingDir}{config.output_page}", "a") as output_file:
            output_file.write("<h4><font face=Verdana><u>Running Messages:</u></h4></font>\n")
            print_initial_running_progress(config)
            output_file.write("<div id='includedContent'>{{ progress_report|safe }}</div>\n")

    if config.PROGRAM == "GUIDANCE":
        run_guidance(config)

    elif config.PROGRAM == "HoT":
        run_hot(config)

    elif config.PROGRAM in ["GUIDANCE3", "GUIDANCE3_HOT"]:
        run_guidance3(config)

    calculate_sp_scores(config)
    modify_score_files_for_codons_and_server(config)
    remove_sites(config)
    prepare_plots(config)
    add_original_seq_names_to_the_MSA(config)
    remove_sequences_sp_score(config)
    remove_sequences_sp_and_z_score(config)
    make_jalview(config)
    create_tar_archives(config)
    if config.PROGRAM == "GUIDANCE3":
        select_best_msa(config)
    flag_that_finished_ok(config)
    print(config)
    if config.isServer == 1:
        prepare_rerun_parameters(config)
        send_finish_email_to_user(config)

    if os.path.exists(os.path.join(RESULTS, "timestamps.txt")):
        shutil.move(os.path.join(RESULTS, "timestamps.txt"), os.path.join(config.WorkingDir, "timestamps.txt"))


if __name__ == "__main__":
    main()

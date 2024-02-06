import sys
import os
import time
import shutil
import re
import subprocess
import inspect
import tarfile
from time import time, localtime, strftime
import signal
from datetime import datetime
from os import times
from hot_cos_file_handler import *

debug = int(os.environ.get('DEBUG_LEVEL', 0))

# def my_sigtrap(sig):
#     log_print(0,0,f"{os.path.basename(sys.argv[0])} {os.getpid()} caught a SIG{sig} -- {time.strftime('%Y-%m-%d %H:%M:%S')} \n", file_handler)
#     sys.exit(0)


def log_print(level, to, log_entry, file_handler):
    if level > int(debug):
        return

    log_entry = log_entry.rstrip() + "\n"
    package, filename, line, *_ = inspect.stack()[1]

    if not log_entry.startswith("-"):
        log_entry = f"@ line {line} of {filename}\n{log_entry}@---\n"

    # if log_file:
    with open(file_handler.log_file, "a") as log_file_handler:
        log_file_handler.write(log_entry)

    if to > 0:
        sys.stdout.write(log_entry)

    if to > 1:
        sys.stderr.write(log_entry)

def write_to_file(filename, content):
    with open(filename, 'w') as file:
        file.write(content)

def run_command_line(command_line, file_handler):
    package, filename, line = inspect.stack()[1][1:4]
    log_print(0, debug, f"---- sh: line {line} of {filename}\n{command_line}\n", file_handler)
    rc = subprocess.getoutput(command_line)
    log_print(1, debug - 1, f"---- output:\n{rc}\n----\n", file_handler)
    return rc

def cleanup(state, file_handler):
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    with open(file_handler.time_file, 'a') as time_file_handler:
        time_file_handler.write(f"Total: {','.join(map(str, times()))},{time() - file_handler.t0}\n{localtime()}\n")

    if state == 0:
        if debug < 4: run_command_line(f"rm -f prof* tr* in* *.atsaf temp* pre", file_handler)
        os.chdir("..")
        run_command_line(f"mkdir -p {file_handler.output_base_dir}", file_handler)

        if file_handler.should_create_archive:
            run_command_line(f"tar -czf {file_handler.done_tar_file} {file_handler.output_dir}; rm -rf {file_handler.output_dir}", file_handler)
            log_print(0, 1,
                      f"---\n {file_handler.output_dir} {file_handler.current_script_file} done : dir saved to {file_handler.done_tar_file}\n{strftime('%Y-%m-%d %H:%M:%S', localtime())}\n", file_handler)
        else:
            if file_handler.basedir not in file_handler.output_base_dir:
                run_command_line(f"rm -rf {file_handler.output_base_dir}/{file_handler.output_dir}; mv -f {file_handler.output_dir} {file_handler.output_base_dir}", file_handler)
            log_print(0, 1,
                      f"---\n {file_handler.output_dir} {file_handler.current_script_file} done : dir saved to {file_handler.output_base_dir}/{file_handler.output_dir}\n{strftime('%Y-%m-%d %H:%M:%S', localtime())}\n", file_handler)
        sys.exit()
    else:
        os.chdir("..")
        run_command_line(f"mkdir -p {file_handler.output_base_dir}_err", file_handler)

        if file_handler.should_create_archive:
            run_command_line(f"tar -czf {file_handler.err_tar_file} {file_handler.output_dir}", file_handler)
            log_print(0, 1,
                      f"{file_handler.output_dir} {file_handler.current_script_file} error : tmp dir saved to {file_handler.err_tar_file}\n{strftime('%Y-%m-%d %H:%M:%S', localtime())}\n", file_handler)
        else:
            run_command_line(f"cp -rf {file_handler.output_dir} {file_handler.output_base_dir}_err", file_handler)
            log_print(0, 1,
                      f"{file_handler.output_dir} {file_handler.current_script_file} error : tmp dir saved to {file_handler.output_base_dir}_err/{file_handler.output_dir}\n{strftime('%Y-%m-%d %H:%M:%S', localtime())}\n", file_handler)

        # if config.debug < 4: run_command_line(f"rm -rf {file_handler.output_dir}")
        sys.exit()

def print_usage():
    usage = f"""\n\nCOS v2.05\n\nUsage:\n
              \n   {os.path.basename(sys.argv[0])} case_id msa_method seq_type input_fasta_file output_dir [MSA_Status File] ['tgz'] [msa_program_path] [user guide_tree] [split number]--- [Additional MSA program parameters]\n\n
              \tmsa_method= \t 'CLW' : ClustalW2
                            \t 'MFT' : mafft
                            \t 'PRK' : prank
                            \t          Append 'h' to do just HoT, e.g.: 'CLWh'\n\n
              \tseq_type=   \t 'aa' | 'nt'\n\n
              \tsplit_number=\t'all' | specific split number\n
              \tSuccessful output at: (output_dir)/(case_id)_cos_(msa_method)/ or (output_dir)/(case_id)_cos_(msa_method).tgz\n
              \tError output at: (output_dir)_err/(case_id)_cos_(msa_method)_err/ or (output_dir)_err/(case_id)_cos_(msa_method)_err.tgz\n
            """
    print(usage + "\n")

def handle_termination_signal(file_handler):
    def my_sigtrap(sig):
        log_print(0, 0,
                  f"{file_handler.current_script_file} {os.getpid()} caught a SIG{sig} -- {time.strftime('%Y-%m-%d %H:%M:%S')} \n",
                  file_handler)
        sys.exit(0)

    signal.signal(signal.SIGTERM, my_sigtrap)

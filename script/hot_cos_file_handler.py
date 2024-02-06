import os
import sys
from time import time, localtime, strftime

class FileHandler:
    def __init__(self, input_file, output_base_dir , status_file, output_dir, tgz):
        self.input_file = os.path.abspath(input_file)
        self.output_dir = os.path.join(os.getcwd(),output_dir)
        self.output_base_dir = os.path.abspath(output_base_dir)
        self.status_file = status_file
        self.log_file = os.path.join(os.getcwd(),output_dir, f"{output_dir}.log")
        self.time_file = 'times.txt'
        self.done_tar_file = os.path.join(output_base_dir, f"{output_dir}.tgz")
        self.err_tar_file = os.path.join(output_base_dir + '_err', f"{output_dir}_err.tgz")
        self.basedir = os.path.abspath('./')
        self.t0 = time()
        self.should_create_archive = 1 if tgz and tgz.lower().startswith('t') else 0
        self.current_script_file = os.path.basename(sys.argv[0])


        if not os.path.exists(input_file):
            print(f"\nERROR: File not found: {input_file}\n")
            # print_usage()
            sys.exit()

    # Instance method
    # def func(self):
    #     pass
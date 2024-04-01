import re
import sys

from hot_cos_logger import print_usage
from hot_cos_mafft_functions import align_sequences_MAF, align_profiles_MAF, met_init_MAF, make_guide_tree_MAF
from hot_cos_prank_functions import align_sequences_PRK, align_profiles_PRK, met_init_PRK, make_guide_tree_PRK


class SequencingMethod:
    # MSA methods vocabulary (from method to command?)
    method_to_cmd_mapping = {
        'CLW': 'CW2',
        'CW2': 'CW2',
        'CW3': 'CW2',
        'MFT': 'MAF',
        'MFM': 'MAF',
        'PRK': 'PRK'
    }
    def __init__(self, sequencing_name, msa_program_path, parameters):
        # if sequencing_name not in self.method_to_cmd_mapping:
        #     print(f"\nERROR: Unknown msa method {sequencing_name}\n{config.usage}")
        #     sys.exit()
        if sequencing_name not in self.method_to_cmd_mapping:
            print(f"\nERROR: Unknown MSA method {sequencing_name}\n")
            print_usage()
            sys.exit()
        # self.name = re.match(r'^(.{3})(.?)', sequencing_name).group(1)
        self.name = sequencing_name
        # self.run_only_hot = re.match(r'^(.{3})(.?)', sequencing_name).group(2)
        self.command = self.method_to_cmd_mapping[self.name]
        self.path = msa_program_path
        self.parameters = ' '.join(parameters + [' ']).replace('---', '')
        self.version = None
        self.prkver = None


    # Instance method
    def met_init(self, seqtype, file_handler):
        if self.name == 'MFT' or self.name == 'MFM':
            return met_init_MAF(seqtype, self, file_handler)
        elif self.name == 'PRK':
            met_init_PRK(seqtype, self, file_handler)
        # elif self.name == 'CLW' or self.name == 'CW2' or self.name == 'CW3':
        #     met_init_CW2(seqtype, msa_program_path)

    def make_guide_tree(self, infile, treefile, sequence, file_handler):
        if self.name == 'MFT' or self.name == 'MFM':
            return make_guide_tree_MAF(infile, treefile, self.version, sequence, file_handler)
        if self.name == 'PRK':
            return make_guide_tree_PRK(infile, treefile, self, sequence, file_handler)
    def align_sequences(self, infile, treefile, outfile, sequence, file_handler):
        if self.version is None:
            print(f"\nERROR: MSA method version is EMPTY\n")
            sys.exit()
        if self.name == 'MFT' or self.name == 'MFM':
            return align_sequences_MAF(infile, treefile, outfile, self.version, sequence, file_handler)
        if self.name == 'PRK':
            return align_sequences_PRK(infile, treefile, outfile, self, sequence, file_handler)
    def align_profiles(self,pfile1, pfile2,tfile1, tfile2, tfile3, ofile, sequence, file_handler):
        if self.path is None:
            print(f"\nERROR: MSA program path is EMPTY\n")
            sys.exit()
        if self.name == 'MFT' or self.name == 'MFM':
            return align_profiles_MAF(pfile1, pfile2, ofile, self.path, sequence, file_handler)
        if self.name == 'PRK':
            return align_profiles_PRK(pfile1, pfile2, tfile3, ofile, self, sequence, file_handler)

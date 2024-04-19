import os
import sys

from hot_cos_logger import log_print, cleanup
import re

class Sequence:
    def __init__(self, sequence_type, sequence_file_name):
        self.sequence_type = 1 if sequence_type.lower().startswith('n') else 0
        self.sequence_heads = None
        self.heads_sequence_file_name = sequence_file_name
        self.sequence_tails = None
        self.tails_sequence_file_name = None
        self.bidirectional_sequences_manager = {'assigned_ids': [], 'original_ids': [], 'sequences': [[], []], 'fasta': [[], []]}
        self.file_extensions = ['.fasta', '.atsaf']
        self.hot_direction = ['h', 't']


    # Instance method
    def construct_bidirectional_sequence_manager(self, file_handler):
        """Read fasta file and Fill in the bidirectional sequence manager with the original sequences and reverse sequences fasta format strings,
          as well as arrays with the relevant sequences and their corresponding IDs"""
        # Read fasta file, write sequence names file, fill sequence structure hot_seqs

        if not os.path.exists(self.heads_sequence_file_name):
            log_print(0, 2, f"ERROR: {file_handler.current_script_file} $$ : File not found: {self.heads_sequence_file_name}", file_handler)
            cleanup(1, file_handler)

        # catch names containing '>'
        with open(self.heads_sequence_file_name, 'r') as heads_sequence_file_handler:
            # delimiter = ">"
            record = heads_sequence_file_handler.read().split('>')
            if record[0] != '':
                heads_sequence_file_handler.close()
                log_print(0, 2, f"ERROR: {file_handler.current_script_file} $$ : File not in fasta format: {self.heads_sequence_file_name}", file_handler)
                cleanup(1, file_handler)
            names_txt = ''
            i = 0
            for element in record[1:]:
                sn = f'seq{i:04d}'
                self.bidirectional_sequences_manager['assigned_ids'].append(sn)
                name, seq = element.split('\n', 1)
                self.bidirectional_sequences_manager['original_ids'].append(name)
                names_txt += f"{sn} {name}\n"
                seq = ''.join(seq.split())
                seq = seq.upper()
                self.bidirectional_sequences_manager['sequences'][0].append(seq)
                # every 60 characters in the seq1 string with the same substring followed by a newline character
                seq1 = '\n'.join([seq[k:k + 60] for k in range(0, len(seq), 60)]) + '\n'
                self.bidirectional_sequences_manager['fasta'][0].append(f">{sn}\n{seq1}")
                seq = seq[::-1]
                self.bidirectional_sequences_manager['sequences'][1].append(seq)
                seq = '\n'.join([seq[k:k + 60] for k in range(0, len(seq), 60)]) + '\n'
                self.bidirectional_sequences_manager['fasta'][1].append(f">{sn}\n{seq}")
                i += 1
        output_sequence_names_file = 'seq_names.txt'
        self.bidirectional_sequences_manager['number_of_sequences'] = i
        with open(output_sequence_names_file, 'w') as output_file_handler:
            output_file_handler.write(names_txt)

        return self.bidirectional_sequences_manager
    def reverse_sequence(self, file, mode, file_handler):
        """Check if it is a valid sequence and return string with Tails (reversed) version of the sequence or if mode=1 return 0 in case of a failure"""
        # Read fasta MSA file and check that sequences match global %hot_seqs
        # Reverse residue order of sequences and return in string
        # reverse_sequence(self, file, mode, file_handler)
        # mode=[0:return tails or shutdown on error, |
        #       1:just check, dont shutdown and return boolean, this is for resuming aborted runs]

        pid = os.getpid()
        if not os.path.exists(file):
            if mode:
                return 0
            log_print(0, 2, f"ERROR: {file_handler.current_script_file} {pid} : File not found: {file}", file_handler)
            cleanup(1, file_handler)

        with open(file, "r") as infile:
            records = infile.read().split('>')
            if records[0] != '':
                if mode:
                    return 0
                log_print(0, 2, f"ERROR: {file_handler.current_script_file} {pid} : File not in fasta format: {file}", file_handler)
                cleanup(1, file_handler)
            sequence_direction = 0 if file.endswith(self.file_extensions[0]) else 1
            fastar = ''
            for record in records[1:]:
                lines = record.split('\n', 1)
                name, seq = lines[0], lines[1].replace('\n', '').replace('\r', '')
                # reverse
                seqr = seq[::-1]
                # split into groups of 60 characters
                seqr = '\n'.join([seqr[i:i + 60] for i in range(0, len(seqr), 60)])
                fastar += f">{name}\n{seqr}\n"

                sid = re.search(r'\d{4}$', name).group(0)
                seq = seq.replace('-', '')
                log_print(6, 0,
                          f"{name} {sequence_direction} {sid}:\n{seq}\n ref:\n{self.bidirectional_sequences_manager['sequences'][sequence_direction][int(sid)]}\n-----\n", file_handler)

                if seq != self.bidirectional_sequences_manager['sequences'][sequence_direction][int(sid)]:
                    # infile.close()
                    if mode:
                        return 0
                    log_print(0, 2,
                              f"ERROR: {file_handler.current_script_file} {pid} : Sequence mismatch in file: {file} {name} :\n---\n{seq}\nshould be:\n{self.bidirectional_sequences_manager['sequences'][sequence_direction][int(sid)]}\n---\n", file_handler)
                    cleanup(1, file_handler)
        return fastar



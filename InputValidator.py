import os
import sys
import logging
import re
import Bio
from Bio.Data import CodonTable

if os.path.exists('/home/josefspr/bioseq'):  # remote run
    sys.path.insert(0, '/home/josefspr/bioseq/guidance/guidance.v2.02/www/Guidance')

# this line is used and important here
import SharedConsts as CONSTS
from utils import *

_VALID_NUC = frozenset('ACGTURYSWKMBDHVNX')

class InputValidator:
    def __init__(self):
        return
    def validate_Seqs(self, working_dir, seqFile, seqType, isMSA, codonTable=None):

        seqFilePath = os.path.join(working_dir, seqFile)
        seqFile_fixed = f'{seqFile}.FIXED'
        seqFilePath_fixed = os.path.join(working_dir, seqFile_fixed)

        # ── Pass 1: read & validate (no writing yet) ──────────────────────────
        try:
            with open(seqFilePath, "r", encoding="ISO-8859-1") as f_in:
                raw_lines = f_in.readlines()
        except IOError:
            return 'sys_error', f"Validate_Seqs:Can't open {seqFilePath}"

        seq = ''
        seq_name = ''
        seq_length = 0
        counter = 0
        warning = ''
        errors = ""
        valid_seqs = []  # (name, seq) pairs that passed all checks

        def _validate_one(seq_name, seq):
            nonlocal errors, warning, counter, seq_length
            if isMSA:
                if seq_length == 0:
                    seq_length = len(seq)
                if len(seq) != seq_length:
                    errors += (f"The sequences of the provided MSA are not properly aligned, "
                               f"For example the seq: '{seq_name}' is not aligned with all others. "
                               f"Please fix the alignment and run GUIDANCE again or provide GUIDANCE sequences only<br>")
                if seqType == "Codons":
                    ans = self.validate_seq_in_CodonAlign(seq, seq_name, codonTable)
                    if ans != 'OK':
                        errors += ans
            if not isMSA:
                m = re.search(r'([-]+)$', seq)
                if m:
                    seq = re.sub(m.group(1), '', seq)
                    warning = "Gap characters (-) were removed from the end of the sequences"
                if re.search('[-]', seq):
                    errors += (f"Seq: named '{seq_name}' contain a gap character '-' which is illegal "
                               f"when sequences are submited to GUIDANCE. If you intended to submit an "
                               f"alignment, please upload the file using the 'Upload MSA file for evaluation' option<br>")
                if seqType == "Codons" and codonTable:
                    ans_codon = self.validate_seq_in_CodonAlign(seq, seq_name, codonTable)
                    if ans_codon != 'OK':
                        errors += ans_codon
            if re.search(r'\*+', seq):
                seq = re.sub(r'\*+', '', seq)
                warning = "Star character (*) were removed from the end of the sequences"
            ans = self.validate_single_seq(seq_name, seq, seqType)
            if ans == 'OK':
                valid_seqs.append((seq_name, seq))
                counter += 1
            else:
                errors += ans
            return seq

        for line in raw_lines:
            line = line.rstrip()
            line = re.sub(r'^\s+|\s+$', '', line)
            if line == '':
                continue
            if line[0] != '>':
                seq += line
            else:
                m = re.search(r'^>(.*)', line)
                if m:
                    if seq == '' and seq_name != '':
                        errors += f"The sequence named '{seq_name}' is missing<br>"
                    if seq != '' and seq_name != '':
                        seq = _validate_one(seq_name, seq)
                    seq_name = m.group(1)
                    seq_name = re.sub(r'^\s+|\s+$', '', seq_name)
                    if seq_name == '':
                        seqNum = counter + 1
                        errors += f"Seq number {seqNum} has no sequence name; Please fix and resubmit<br>"
                    else:
                        seq = ''

        # validate last sequence
        if seq == '' and seq_name != '':
            errors += f"The sequence named '{seq_name}' is missing<br>"
        elif seq_name != '':
            _validate_one(seq_name, seq)

        # ── Return validation errors without touching the FIXED file ──────────
        if errors:
            try:
                with open(f'{working_dir}/errors.txt', "w") as ef:
                    ef.write(errors.replace("<br>", "\n"))
            except IOError:
                pass
            return errors

        # ── Pass 2: write FIXED only after validation passed ─────────────────
        try:
            with open(seqFilePath_fixed, "w", encoding="ISO-8859-1") as f_out:
                for sname, sseq in valid_seqs:
                    f_out.write(f">{sname}\n{sseq}\n")
        except IOError:
            return 'sys_error', f"Validate_Seqs:Can't open {seqFilePath_fixed} for writing"

        return 'OK', warning, seqFile_fixed, counter

    def validate_single_seq(self, seqName, seq, seqType):

        if seqType == 'AminoAcids' and not re.search('[ABRNDCQEGHILKMFPSTWYVXZabrndcqeghilkmfpstwyvxz]+', seq): 
            return f"Seq: '{seqName}' is empty<br>"
        elif seqType != 'AminoAcids' and  not re.search('[ACTGUNactgun]+', seq): 
            return f"Seq: '{seqName}' is empty<br>"
            
        m = re.search('([^ABRNDCQEGHILKMFPSTWYVXZabrndcqeghilkmfpstwyvxz-])', seq)
        if seqType == 'AminoAcids' and m: #Maybe allow: _*-?
            return f"Seq: '{seqName}' contained the character '{m.group(1)}', which is not a standard Amino Acid<br>"
            
        m = re.search('([^ACGTRYWSMKHBVDNUXacgtrywsmkhbvdnux-])', seq)
        if seqType != 'AminoAcids' and m: 
            wrong_char = m.group(1)
            if re.search('[Uu]', seq) and seqType == 'Nucleotides': 
                return f"Currently GUIDANCE does not accept 'U's in nucleotide sequences, you may consider replacing the 'U's by 'T's and re-submit. <br> In addition, seq: '{seqName}' contained the character '{wrong_char}', which is not a standard Nucleotide <br>"
            return f"Seq: '{seqName}' contained the character '{wrong_char}', which is not a standard Nucleotide<br>"
            
        m = re.search('[Uu]', seq)
        if seqType == "Nucleotides" and m: #Maybe allow: _*-?
            return f"Currently GUIDANCE does not accept 'U's in nucleotide sequences, you may consider replacing the 'U's by 'T's and re-submit.<br>"
        return 'OK'

    def validate_seq_in_CodonAlign(self, DNA_seq, seqName, codonTableIndex):
    
        AASeq=""
        codonTable_obj = CodonTable.generic_by_id[int(codonTableIndex)]
        
        DNA_seq = DNA_seq.rstrip()
        seq_length = len(DNA_seq)
        if seq_length % 3 > 0:
            return f"Sequence '{seqName}' is not a valid coding sequence: the sequence is of length {seq_length} which is not divided by 3\n"
        
        i = 0
        while i < seq_length - 2:

            codon = DNA_seq[i: i+3]
            if codon == '---':
                AA = '-'
            else:
                AA = self.translate_codon(codonTable_obj, codon)
            AASeq = AASeq + AA
            i = i + 3

        if '?' in AASeq:
            return f"Sequence: '{seqName}' contains non-DNA characters - did you submit amino acid sequences as Codons? Please select the correct sequence type and resubmit\n"
        if '*' in AASeq:
            return f"Sequence: '{seqName}' contains a stop codon, please remove all stop codons (from all sequences) and submit to GUIDANCE again\n"
        return 'OK'

    def translate_codon(self, table, codon):

        codon = codon.upper()
        if not all(c in _VALID_NUC for c in codon):
            return '?'  # truly non-DNA character (e.g. amino acid-specific: E, F, L, Q, P)
        if codon in table.stop_codons:
            return '*'
        elif codon in table.forward_table:
            return table.forward_table[codon]
        else:
            return 'X'  # valid IUPAC ambiguous codon - not an error

    def countSeq(self, seqFile):
        numSeqs = 0
        try:
            with open(seqFile, "r", encoding="ISO-8859-1") as f:
                for line in f:
                    if line[0] == '>':
                        numSeqs += 1
            f.close()
            return numSeqs
        except:
            raise Exception( f"InputValidator.countSeq: Error in sequence file {seqFile}", "user")

    def get_max_seq_length(self, seqFile):

        longest_seq = 0
        try: 
            with open(seqFile, "r", encoding="ISO-8859-1") as f:
                try: 
                    line = f.readline()
                    while line: 
                        # next line is sequence
                        line = f.readline()
                        seq = ''
                        while str_defined_and_not_equal_char(line, '>'):
                            line = line.rstrip()
                            seq = seq + line
                            line = f.readline()
                        if longest_seq == 0 or len(seq) > longest_seq: 
                            longest_seq = len(seq)
                except:
                    raise Exception (f"InputValidator.get_max_seq_length: error in sequence file", "system")
            f.close()
            return longest_seq
            
        except:
           raise Exception (f"InputValidator.get_max_seq_length: can't open FASTA: '{seqFile}'", "system")

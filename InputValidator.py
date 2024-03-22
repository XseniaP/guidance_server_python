import os
import sys
import logging
import re
import Bio

if os.path.exists('/home/josefspr/bioseq'):  # remote run
    sys.path.insert(0, '/home/josefspr/bioseq/guidance/guidance.v2.02/www/Guidance')

import SharedConsts as CONSTS  
from utils import *

class InputValidator: 

    def validate_Seqs (working_dir, seqFile, seqType, isMSA, codonTable=None):

        seqFilePath = os.path.join( working_dir, seqFile)
        seq = ''
        seq_name = ''
        seq_length = 0
        counter = 0
        warning = ''
        errors = ""
        
        try: 
            with open (seqFilePath, "r", encoding="ISO-8859-1") as f_in: 
                try:
                    seqFile_fixed = f'{seqFile}.FIXED'
                    seqFilePath_fixed = os.path.join( working_dir, seqFile_fixed)
                    #print (f'validate_Seqs: opening {seqFilePath_fixed}')
                    with open (seqFilePath_fixed, "w", encoding="ISO-8859-1") as f_out:
                        for line in f_in:
                            # remove newline, leading and trailing blanl spaces
                            line = line.rstrip()
                            line = re.sub( r'^\s+|\s+$', '', line)
                            if line == '':
                                continue
                            # append sequences separated by newline
                            if line[0] != '>':
                               seq += line
                            else:
                                m = re.search( r'^>(.*)', line)
                                if m: 
                                    # validate previous seq
                                    if seq == '' and seq_name != '':
                                        # return f"The sequence named '{seq_name}' is missing<br>"
                                        errors += f"The sequence named '{seq_name}' is missing<br>"
                                        
                                    if seq != '' and seq_name != '':
                                        # validate seq according to if is MSA or not
                                        if isMSA: 
                                            # Make sure alignment length equal
                                            if seq_length == 0: 
                                                seq_length = len(seq)# initialize the first one
                                            if len(seq) != seq_length:
                                                # return f"The sequences of the provided MSA are not properly aligned, For example the seq: '{seq_name}' does not aligned to all others. Please fix the alignment and run GUIDANCE again or provide GUIDANCE sequences only<br>"
                                                errors += f"The sequences of the provided MSA are not properly aligned, For example the seq: '{seq_name}' does not aligned to all others. Please fix the alignment and run GUIDANCE again or provide GUIDANCE sequences only<br>"
                                            if seqType == "Codons": 
                                                # Make sure that in Codon Alignment there are no stop Codons and all seq are divided by 3
                                                ans = InputValidator.validate_seq_in_CodonAlign( seq, seq_name, codonTable)
                                                if ans != 'OK': 
                                                    # return ans
                                                    errors += ans
                                        if not isMSA:
                                            # check gap characters
                                            m = re.search (r'([-]+)$', seq)
                                            if m:
                                                seq = re.sub (m.group(1), '', seq)
                                                warning = "Gap characters (-) were removed from the end of the sequences"
                                            if re.search('[-]', seq):
                                                # return f"Seq: named '{seq_name}' contain a gap character '-' which is illegal when sequences are submited to GUIDANCE. If you intended to submit an alignment, please upload the file using the 'Upload MSA file for evaluation' option<br>"
                                                errors += f"Seq: named '{seq_name}' contain a gap character '-' which is illegal when sequences are submited to GUIDANCE. If you intended to submit an alignment, please upload the file using the 'Upload MSA file for evaluation' option<br>"
                                        if re.search('\*+', seq):
                                            seq = re.sub ('\*+', '', seq)
                                            warning = "Star character (*) were removed from the end of the sequences"
                                        ans = InputValidator.validate_single_seq( seq_name, seq, seqType)
                                        if ans == 'OK':
                                            f_out.write(f">{seq_name}\n")  # prev seq
                                            f_out.write(f"{seq}\n")        # prev seq
                                            counter += 1
                                        else:
                                            # return ans
                                            errors += ans
                                            
                                    # start new seq
                                    m = re.search( r'^>(.*)', line)
                                    if m:
                                        seq_name = m.group(1)
                                        seq_name = re.sub('^\s+|\s+$', '', seq_name) # remove leading/trailing blanks
                                        if seq_name == '':
                                            seqNum = counter + 1
                                            # return "Seq number {seqNum} has no sequence name; Please fix and resubmit<br>"
                                            errors += "Seq number {seqNum} has no sequence name; Please fix and resubmit<br>"
                                        else:
                                            seq = ''
                                            
                        # end of loop: validate last seq
                        #print (f'validate_Seqs: end of loop')
                        if seq == '' and seq_name != '':
                            # return f"The sequence named '{seq_name}' is missing<br>"
                            errors += f"The sequence named '{seq_name}' is missing<br>"
                        else:
                            # validate seq according to if is MSA or not
                            if isMSA: 
                                # Make sure alignment length equal
                                if seq_length == 0: 
                                    seq_length = len(seq)# initialize the first one
                                if len(seq) != seq_length:
                                    # return f"The sequences of the provided MSA are not properly aligned, For example the seq: '{seq_name}' does not aligned to all others. Please fix the alignment and run GUIDANCE again or provide GUIDANCE sequences only<br>"
                                    errors += f"The sequences of the provided MSA are not properly aligned, For example the seq: '{seq_name}' does not aligned to all others. Please fix the alignment and run GUIDANCE again or provide GUIDANCE sequences only<br>"
                                if seqType == "Codons": 
                                    # Make sure that in Codon Alignment there are no stop Codons and all seq are divided by 3
                                    ans = InputValidator.validate_seq_in_CodonAlign( seq, seq_name, codonTable)
                                    if ans != 'OK': 
                                        # return ans
                                        errors += ans
                            if not isMSA:
                                # check gap characters
                                m = re.search (r'([-]+)$', seq)
                                if m:
                                    seq = re.sub (m.group(1), '', seq)
                                    warning = "Gap characters (-) were removed from the end of the sequences"
                            if re.search('\*+', seq):
                                seq = re.sub ('\*+', '', seq)
                                warning = "Star character (*) were removed from the end of the sequences"
                            ans = InputValidator.validate_single_seq( seq_name, seq, seqType)
                            if ans == 'OK':
                                f_out.write(f">{seq_name}\n")  # prev seq
                                f_out.write(f"{seq}\n")        # prev seq
                                counter += 1
                            else:
                                # return ans
                                errors += ans
                    f_out.close()
                except:
                    error = f"Validate_Seqs:Can't open {seqFilePath_fixed} for writing"
                    return 'sys_error', error
                    
            f_in.close()
            
        except: 
            error = f"Validate_Seqs:Can't open {seqFilePath}"
            return 'sys_error', error
        if errors != "":
            try:
                f = open(f'{working_dir}/errors.txt', "w")
                f.write(errors.replace("<br>","\n"))
                f.close()
            except:
                error = f"Validate_Seqs:Can't open {f} for writing"
                return 'sys_error', error
            return errors
        return 'OK', warning, seqFile_fixed, counter

    def validate_single_seq (seqName, seq, seqType):

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

    def validate_seq_in_CodonAlign ( DNA_seq, seqName, codonTableIndex):
    
        AASeq=""
        codonTable_obj = Bio.Data.CodonTable.generic_by_id[codonTableIndex]
        
        DNA_seq = DNA_seq.rstrip()
        seq_length = len(DNA_seq)
        if seq_length % 3 > 0:
            return f"Sequence '{seqName}' is not a valid coding sequence: the sequence is of length {seq_length} which is not divided by 3"
        
        i = 0
        while i < seq_length - 2:

            codon = DNASequence[i, i+3]
            if codon == '---':
                AA = '-'
            else:
                AA = translate_codon(codonTable_obj, codon)
            AASeq = AASeq + AA
            i = i + 3

        if '*' in AASeq:
            return f"Sequence: '{seqName}' contains a stop codon, please remove all stop codons (from all sequences) and submit to GUIDANCE again"
        return 'OK'

    def translate_codon (table, codon):

        codon = codon.upper()
        if codon in table.stop_codons:
            return '*'
        else:
            return forward_table[codon]
        
    def countSeq(seqFile): 
        numSeqs = 0
        try:
            with open( seqFile, "r", encoding="ISO-8859-1") as f:
                for line in f:
                    if line[0] == '>':
                        numSeqs += 1
            f.close()
            return numSeqs
        except:
            raise Exception( f"InputValidator.countSeq: Error in sequence file {seqFile}", "user")

    def get_max_seq_length (seqFile):

        longest_seq = 0
        try: 
            with open (seqFile, "r", encoding="ISO-8859-1") as f: 
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

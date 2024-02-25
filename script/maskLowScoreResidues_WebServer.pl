#!/usr/bin/perl -w

#  MSA.MAFFT.Guidance_res_pair.scr  MASK_OUT_FILE 0.3 aa
#UserMSA.FIXED.ORIG MSA.MAFFT.Guidance_res_pair_res.scr  MASK_OUT_FILE 0.1 aa
#UserMSA.FIXED.With_Names  MSA.MAFFT.Guidance_res_pair_res.scr  Mask_Out_File 0.3 aa
#MSA.MAFFT.Guidance_res_pair_res.scr

use strict;
use Storable;
use FindBin qw($Bin);

use lib $Bin; # == www/Guidance
use Guidance;

use Bio::SeqIO;
use Bio::Align::AlignI;
use Bio::AlignIO;
use warnings;



my $stored_data_file=shift;
my $alphabet=shift;
my $cutoff=shift;

my $vars_ref = retrieve($stored_data_file);
my %VARS = %$vars_ref;


my $msaFile="$VARS{WorkingDir}$VARS{Alignment_File}"; 
my $scoreFile="$VARS{WorkingDir}$VARS{Output_Prefix}_res_pair_res.scr";
my $out_website="Mask_Residues_Res_".$cutoff.".aln";
my $outFile="$VARS{WorkingDir}$out_website";
my $seq_names_index="$VARS{WorkingDir}$VARS{code_fileName}";


my $missingDataChar;
if ($alphabet eq "aa") {
    $missingDataChar = "X";
} elsif ($alphabet eq "nuc") {
    $missingDataChar = "N";
} else { die "ALPHABET must be either 'aa' or 'nuc'\n" }

my $in_fasta = Bio::SeqIO->new(-file => $msaFile, '-format' => 'fasta');
my @seqs;
my @ids;
while (my $seqObj = $in_fasta->next_seq()) {
    my @seq_chars = split(//,$seqObj->seq());
    push(@seqs,\@seq_chars);
    push(@ids,$seqObj->id());
}

#my $in_fasta  = Bio::AlignIO->new(-file => $msaFile , '-format' => 'fasta');
#my $aln = $in_fasta->next_aln();


open IN, "<$scoreFile" or die "can't open file $scoreFile";
#print "cutoff: $cutoff\n";
while (my $line = <IN>) #COL_NUMBER	#ROW_NUMBER	#RES_PAIR_RESIDUE_SCORE
{
    chomp $line;
    # print $line;
    if ($line =~ m/^#/)
	{
	    next;
	}
	if ($line =~ m/^\s*(\d+)\s+(\d+)\s+(\S+)$/) 
	{
	    if ($3 ne 'nan' and $3 < $cutoff) 
	    {
			my $col=$1-1;
			my $row=$2-1;
			$seqs[$row][$col] = $missingDataChar;
			#warn "DEBUG: masking $row,$col\n";
	    }
	} 
	else 
	{ 
	    warn "WARNING: failed to parse line: '$line'\n";
	}
}
close IN;

my %id_names=();
open (SEQ_INDEX,$seq_names_index) || die "Can't open seq_index: '$seq_names_index' $!";
while (my $line=<SEQ_INDEX>)
{
	chomp $line;
	my ($name,$id)=split(/\t/,$line);
	$id_names{$id}=$name;
}
close SEQ_INDEX;
open OUT, ">$outFile";
for (my $i=0; $i<@seqs; ++$i) 
{
	my $id = $ids[$i];
	my $id_name=$id_names{$id};
	my $seqRef=$seqs[$i];
	my @seq_arr = @$seqRef;
	my $seq = join('',@seq_arr);
	print OUT ">$id_name\n";
	print OUT "$seq\n";
}
close OUT;
chmod (0664,$outFile);
# Update the output page
#######################################
open (OUTPUT,"$VARS{WorkingDir}$VARS{output_page}");
my @out=<OUTPUT>;
close (OUTPUT);
open (OUTPUT,">$VARS{WorkingDir}$VARS{output_page}");
my $Remove_Pos_Section=0;
foreach my $line (@out)
{
    
    if ($line=~/Mask specific residues below a certain cutoff:/)
    {
        print OUTPUT $line ;
		system ("chmod +r $outFile");
		print_message_to_output("<A HREF=$out_website TARGET=_blank>The MSA after masking unreliable residues (below $cutoff)</A>"); 
	}
    else
    {
		print OUTPUT $line ;
    }
}

close (OUTPUT);


#---------------------------------------------
sub print_message_to_output{
#---------------------------------------------
    my $msg = shift;
    print OUTPUT "\n<ul><li>$msg</li></ul>\n";
}


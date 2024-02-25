# !/usr/bin/perl

use strict;
use warnings;

# use JSON::XS;
use Storable;
use File::Slurp;
use Data::Dumper;
use JSON::XS;

my $file = shift;
my $outFile = shift;
my $file_content = read_file($file);

my $arrayref = decode_json $file_content;
#print Data::Dumper->Dump([$arrayref], [qw(arrayref)]);

store $arrayref, $outFile;



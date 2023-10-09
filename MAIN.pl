#!/usr/bin/perl

use warnings;
use strict;
use English;
use File::Temp qw/ tempdir /;
use Getopt::Long 'HelpMessage';
use File::Basename;
use Cwd;
use POSIX;

#my $FASTA = $ARGV[0];
my $FASTA="/home/mwyler/GenotapeClassifier/testNEWseq.fa";
# getwd
my $dir = getcwd . "/";
# script folder
my $SCRIPTS = dirname $0;

### INPUT ----------------------------

# make temp folder for working
my $TEMPfolder = tempdir( DIR => $dir, CLEANUP => 1 );

# read in fasta file, check header and put into temp folder
open(IN, $FASTA) or die "can't open $FASTA";
my %seqs = ();
my $header = '';

#open(IN, $FASTA) or die "can't open $FASTA";
while (my $line = <IN>){
    chomp $line;
	$line =~ s/^\s*$//;
    if ($line =~ m/^>(.*)$/){
        $header = $line;
       # Check if segments are annotated (as '|1|)
        if ($header !~ m/\|[[:digit:]]\|/){
                print "How are each segment defined? Not line '|number|'.\n";
                exit;
        }
    } else {
        $seqs{"$header"} .= $line;
    }
}
close (IN);


 ### ALIGNMENTS ----------------------

# copy out alignments from script folder
system "cp $SCRIPTS/RefSegment_*.aln $TEMPfolder/";

# loop over segments (pull out from hash and find alignment
my @NRsegments = (1..8);
foreach (@NRsegments){
        my $NUMERO = $_;
	# name alignment file
        my $ALN = "$TEMPfolder/RefSegment_" . $NUMERO . ".aln";
        # get right sequence out of hash
        my $SEQSTRING = "|" . $NUMERO . "|";
	my @ALLsequences = keys %seqs;
	my $SEQNAME ;
	foreach (@ALLsequences) {if (index($_, $SEQSTRING) != -1) { $SEQNAME = $_}};
	# print fasta segment
        my $SEQfile = $TEMPfolder . "/segment_" . $NUMERO . ".fa";
	my $sequenza = $seqs{"$SEQNAME"};
	open(FH, '>', $SEQfile) or die $!;
        print FH "$SEQNAME\n$sequenza\n";
	close(FH);
	# attach sequence to alignment
	system "muscle -quiet -profile -in1 $ALN -in2 $SEQfile -out $TEMPfolder/combined${NUMERO}_comb.aln";
}


# make a new tree
system "ls $TEMPfolder/combined*_comb.aln | parallel -j8 'fasttree -quiet -nt {} > {.}.tree'  > /dev/null 2>&1 ";

# calculate distance matrix
system "ls $TEMPfolder/combined*.tree | parallel -j8 '~/Software/gotree matrix -i {} -o {.}.distMatrix' ";

### CLASSIFICATION -------------------

# run R script to get closest
system "Rscript $SCRIPTS/subScript.R $SCRIPTS $TEMPfolder";

system "rm -r uscita";
system "cp -r $TEMPfolder uscita";

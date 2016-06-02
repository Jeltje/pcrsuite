#! /usr/bin/perl -wT

# Copyright Notice and Disclaimer for PCR Suite
#
# Copyright (c) 2003 Erasmus MC Rotterdam. All rights reserved.
#
# The PCR Suite is based on the Primer3 program of the Whitehead Institute (Copyright notice and Disclaimer: see below). The Primer3 core program is identical to the original software distributed by the Whitehead Institute. The HTML forms used in PCR Suite are based on the Whitehead's Primer3 HTML form, but are not identical. The PCR Suite CGI script is new and property of the Erasmus MC.
#
# Redistribution and use in source and binary forms of the PCR Suite CGI script, with or without modification, are permitted provided that the following conditions are met:
#
#   1. Redistributions must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. Redistributions of source code must also reproduce this information in the source code itself.
#   2. If the program is modified, redistributions must include a notice (in the same places as above) indicating that the redistributed program is not identical to the version distributed by Whitehead Institute.
#   3. All advertising materials mentioning features or use of this software must display the following acknowledgment:
#      This product includes software developed by the Erasmus Medical Centre.
#   4. The name of the Erasmus MC may not be used to endorse or promote products derived from this software without specific prior written permission. 
#
# The code of the script is available at http://www.eur.nl/fgg/primer3/source.txt.
#
# THE PCR SUITE SCRIPT IS PROVIDED BY THE ERASMUS MC ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ERASMUS MC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# Copyright Notice and Disclaimer for Primer3
# Copyright (c) 1996,1997,1998 Whitehead Institute for Biomedical Research. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
#   1. Redistributions must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. Redistributions of source code must also reproduce this information in the source code itself.
#   2. If the program is modified, redistributions must include a notice (in the same places as above) indicating that the redistributed program is not identical to the version distributed by Whitehead Institute.
#   3. All advertising materials mentioning features or use of this software must display the following acknowledgment:
#      This product includes software developed by the Whitehead Institute for Biomedical Research.
#   4. The name of the Whitehead Institute may not be used to endorse or promote products derived from this software without specific prior written permission. 

# We also request that use of this software be cited in publications as
# Steve Rozen and Helen J. Skaletsky (2000) Primer3 on the WWW for general users and for biologist programmers. In: Krawetz S, Misener S (eds) Bioinformatics Methods and Protocols: Methods in Molecular Biology. Humana Press, Totowa, NJ, pp 365-386
# (Code available at http://www-genome.wi.mit.edu/genome_software/other/primer3.html.) THIS SOFTWARE IS PROVIDED BY THE WHITEHEAD INSTITUTE ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE WHITEHEAD INSTITUTE BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI::Pretty qw(-unique_headers);

$| = 1;

our @allout = ();	#array for output, will be filled by subroutines
my $params_general = "";
my @results = "";


# Main program
# tests from which form the input comes: genomic, overlap, cDNAs, SNPs
# formats the input from the forms and calls the corresponding 
# subroutines 
# Also creates the HTML output page

my $q = new CGI;

my $ident = $q-> param("INFROM") 
  || error ($q, "It seems you are not using the proper input form");

create_html($q);

# create shared parameters list

if ($q ->param('Pick Primers')){
  $params_general = shared_parameters($q);
} else {
  error ($q, "Did not see the 'Pick Primers' query parameter")
}


# get the user id and add it to a list
my$userAddr=$ENV{REMOTE_ADDR};
my$file="/usr/local/apache/share/eur/htdocs/fgg/kgen/primer/stats.txt";

open (STATS, ">>$file");
print STATS scalar localtime;
print STATS "\t$userAddr\t$ident\n";
close STATS;

# extract specific info and run program for each instance of INPUT

if ($ident eq "genomic"){

# get specific parameters from form

  my$PRIMER_PRODUCT_SIZE_RANGE=  $q-> param("PRIMER_PRODUCT_SIZE_RANGE") 
    ||"250-450";
  my$FLANKING_SIZE = $q-> param("FLANKING_SIZE") || "40";
  if($FLANKING_SIZE > 100){
    error( $q, "Flanking size must be less than 100.");
  }  

  my $gene = $q -> param("gene") || error( $q, "which gene?");

# get data from GenBank file with a separate subroutine

  my $filedata = get_filedata($q);

# correct for capitalization
  (my$geneMatch) = $filedata =~ /gene=\"($gene)\"/i;
  
  error( $q, "Gene $gene is not in input file.") unless $geneMatch;
  $gene = $geneMatch;

  push (@allout,  "<b>primers for $gene<\/b>");

# get the exon sequences with a subroutine

  my@exons = select_exons($filedata, $gene,
    $PRIMER_PRODUCT_SIZE_RANGE);

# send the input to the genomic_primers subroutine, which will
# do more formatting and runs primer3

  @allout = find_genomic_primers($params_general, $FLANKING_SIZE, 
    $PRIMER_PRODUCT_SIZE_RANGE, @exons);
}

elsif ($ident eq "SNP"){

# get specific parameters from form

  my$PRIMER_PRODUCT_SIZE_RANGE=  $q-> param("PRIMER_PRODUCT_SIZE_RANGE") 
    ||"250-450";

  my $filedata = get_filedata($q);
  unless ($filedata =~ /variation/){
    error( $q, "No SNPs in input file.");
  }

# get the SNPs and their flanking sequence

  my@snps = select_range($filedata, $PRIMER_PRODUCT_SIZE_RANGE);

# send the input to the snp_primers subroutine, which will
# do more formatting and runs primer3

  @allout = find_snp_primers($params_general, 
       $PRIMER_PRODUCT_SIZE_RANGE, @snps);
}


elsif ($ident eq "overlap"){

# get specific parameters from form

  my$PRIMER_PRODUCT_SIZE_RANGE=  $q-> param("PRIMER_PRODUCT_SIZE_RANGE") 
    ||"250-450"; 
  my$OVERLAP_SIZE_RANGE = $q-> param("OVERLAP_SIZE_RANGE")|| "30-90";
  my$TARGET = $q-> param("TARGET") ||
    error($q, "Please fill in Target");

# remove heading, numbers and spaces from the inputsequence

  my$seq = clean_sequence($q); 

# send the input to the overlap_primers subroutine, which will
# do more formatting and runs primer3

  @allout = find_overlap_primers($params_general, 
     $PRIMER_PRODUCT_SIZE_RANGE, $OVERLAP_SIZE_RANGE,
     $TARGET, $seq);
}

elsif ($ident eq "cDNA"){

# get specific parameters from form

  my$PRIMER_PRODUCT_SIZE_RANGE=  $q-> param("PRIMER_PRODUCT_SIZE_RANGE")
    ||"250-800"; #check format!
  my $filedata = get_filedata($q);

# separate the cDNAs

  my @filedata = split /\/\/\n/, $filedata;
  pop @filedata;
  my $locus;
  my $file;

  foreach $file(@filedata){
    print $q->p(" ");
    $file =~ s/\s*LOCUS/LOCUS/;
    $file .= "\/\/\n";

# get sequence, name, and the ORF position for each cDNA

    my($dna, $gene, $from, $to) = get_cDNA_ID_and_sequence($file); 
    push (@allout, "<hr>","<b>primers for $gene<\/b>");

# send the input to the cDNA_primers subroutine, which will
# do more formatting and runs primer3

    my@prims = find_cDNA_primers($params_general, 
      $PRIMER_PRODUCT_SIZE_RANGE, $from, $to, $dna);
    @allout = @prims;
  }
}

else {error ($q, "It seems you are not using the proper input form")};

print '<FONT face="courier new" size="2">';
foreach (@allout){
print $q->p($_);
}
print '</FONT>';
print "<hr>";
print $q->end_html;

exit;


##################### subroutines ##########################

############################################################
#                                                          #
#    subroutines for HTML formatting                       #
#                                                          #
############################################################


sub create_html{

# creates the start of the HTML output

  my ($q) = @_;

  print $q-> header( "text/html" );
  print $q-> start_html (-title => "your primers",
                         -bgcolor => "#ffcc66");
  print $q->h2("Your Primers");
  print $q->hr;
  print $q->p("If you want to design other primers use the");
  print $q->a({-href => 
    "//www-genome.wi.mit.edu/cgi-bin/primer/primer3_www.cgi"}, 
    "Primer3 input form");
  print $q->p("Note: these primers may contain repeats.");
  print $q->p("wait while the program is running...");
}

sub error {

# used by main program, clean_sequence and get_filedata
# creates HTML output for error messages

 my( $q, $reason ) = @_;

 print $q->header( "text/html" ),
 $q->start_html( "Error" ),
 $q->h1( "Error" ),
 $q->p( "Your upload was not processed because the following error",
 "occured: "),
 $q->p( $q->i( $reason )),
 $q->end_html;
 exit;
}


sub format_sequence{

# format sequence takes in the primer3 output, a line length, 
# a color for the primer and a color for the product, 
# and a trim flag (A means do not trim)
# it outputs lines of the required length with HTML tags
# for font coloring
# when asked, the sequence is trimmed.
# To run the subroutine on a cDNA (in genomic_primers)
# the results hash should be empty, and the sequence should
# be given as the last argument

  my($src, $results, $linesize, $primerColor, $productColor, $trim, $dna) = @_;

# Create a list of action items for inserting linebreaks and coloring
  my%action = ();

  my$seq = $dna || $$results{'SEQUENCE'};
  if($results){

# lowercase inputsequence
    $seq = lc($seq);

    my$target = $$results{'TARGET'};

    my($start_left, $length_left) = 
       split (/\,/, $$results{'PRIMER_LEFT'});
    my($end_right, $length_right) = 
       split (/\,/, $$results{'PRIMER_RIGHT'});
    $end_right++;
    my($from, $size) = split /,/, $target;
    my$to = $size+$from -1;
    my$start_right = $end_right - $length_right;

    if($src eq "genomic" || $src eq "cdna"){
# uppercase target
      substr($seq, $from-1) =~ tr/[a-z]/[A-Z]/;
      substr($seq, $to) =~ tr/[A-Z]/[a-z]/;
# uppercase primers
      substr($seq, $start_left, $length_left) =~ tr/actg/ACTG/;
      substr($seq, $end_right-$length_right, $length_right) =~ tr/actg/ACTG/; 
    }elsif($src eq "snp"){
 
# uppercase primers
      substr($seq, $start_left, $length_left) =~ tr/actg/ACTG/;
      substr($seq, $end_right-$length_right, $length_right) =~ tr/actg/ACTG/; 
# uppercase all SNPs
      $seq =~ tr/swmkry/SWMKRY/;
    }elsif($src eq "overlap"){

# uppercase product including primers
      substr($seq, $start_left, $$results{'PRIMER_PRODUCT_SIZE'}) =~ tr/[a-z]/[A-Z]/;

    }

# left primer
    if (exists($action{$start_left})){
      $action{$start_left} .= "<font color=$primerColor>";
    }else{
      $action{$start_left} = "<font color=$primerColor>";
    }
    my$end_left = $start_left + $length_left;
    if (exists($action{$end_left})){
      $action{$end_left} .= "</font>";
    }else{
      $action{$end_left} = "</font>";
    }

# product
    my($pstart, $pend);
    if($src eq "genomic" || $src eq "cdna"){
      $pstart = $from -1;
      $pend = $pstart + $size;
    }elsif($src eq "overlap"){
      $pstart = $end_left;
      $pend = $start_right;
    }elsif($src eq "snp"){
# placeholder
      $pstart = "0";
      $pend = "1";
    }
    if (exists($action{$pstart})){
      $action{$pstart} .= "<font color=$productColor>";
    }else{
      $action{$pstart} = "<font color=$productColor>";
    }
    if (exists($action{$pend})){
      $action{$pend} .= "</font>";
    }else{
      $action{$pend} = "</font>";
    }

# right primer
    if (exists($action{$start_right})){
      $action{$start_right} .= "<font color=$primerColor>";
    }else{
      $action{$start_right} = "<font color=$primerColor>";
    }
    if (exists($action{$end_right})){
      $action{$end_right} .= "</font>";
    }else{
      $action{$end_right} = "</font>";
    }

# The trimsize is the number of bases allowed as flanking
# next to the primers. There may not be enough sequence.
# 'A' means do not trim
  unless($trim eq 'A'){
      my$removeEnd = $end_right +$trim;
      unless($removeEnd > length($seq)){
        substr($seq, $removeEnd, (length($seq) - $removeEnd), "");
      } 
      my$remove = $start_left - $trim;
      if($remove <0 ){
        $remove="0";
      }
      if($remove){
        substr($seq, 0, $remove, "");
        my%newaction = ();
        foreach my$key(keys %action){
	  my$newkey = $key - $remove;
          $newaction{$newkey} = $action{$key};
        }
        %action = %newaction;
      }
    }
  }else{
# color the uppercase part
    while ($seq =~ m/[acgt][ACTG]/g){
      $action{pos($seq)-1} ="<font color=blue>"; 
    }
    while ($seq =~ m/[ACTG][actg]/g){
      $action{pos($seq)-1} ="</font>"; 
    }

  }
# from here on the routine is also done for regular cDNA
# add linebreaks
  for (my$i=$linesize; $i<= length($seq); $i+=$linesize){
    if(exists ($action{$i})){
      $action{$i} .= '<br>';
    }else{
      $action{$i} = '<br>';
    }
  }


# create output sequence with correct tags
  my$formatSeq="";
  my$pos = "0";
  foreach my$key(sort numerically (keys %action)){
	$formatSeq .= substr($seq, $pos, ($key-$pos));
	$formatSeq .= "$action{$key}";
	$pos=$key;
  }
  $formatSeq .= substr($seq, $pos, (length($seq)-$pos));

  return $formatSeq;
}

sub format_covered_cdna{

# used by find_overlap_primers
# takes in a list of matched regions, the original (complete) target,
# and the input sequence.
# The input sequence is lowercased outside the target
# and the matched regions are colored blue

  my$seq = pop @_;
  my$target = pop @_;
  my$linesize = pop @_;
  my@matched = @_;
  $seq = lc($seq);

# uppercase all matches
  foreach my$match (@matched){
#    $seq = uc($seq);
    $seq =~ s/$match/\U$match/i;
  }

my %action = ();
# now find the boundaries
  while ($seq =~ m/[acgt][ACTG]/g){
    $action{pos($seq)-1} ="<font color=blue>"; 
  }
  while ($seq =~ m/[ACTG][actg]/g){
    $action{pos($seq)-1} ="</font>"; 
  }

# add linebreaks
  for (my$i=$linesize; $i<= length($seq); $i+=$linesize){
    if(exists ($action{$i})){
      $action{$i} .= '<br>';
    }else{
      $action{$i} = '<br>';
    }
  }

# now convert sequence back to lowercase, and 
# uppercase target

  $seq = lc($seq);
  my($from, $size) = split /,/, $target;
  my$to = $size+$from -1;
  substr($seq, $from-1) =~ tr/actg/ACTG/;
  substr($seq, $to) =~ tr/ACTG/actg/;

# create output sequence with correct tags
  my$formatSeq="";
  my$pos = "0";
  foreach my$key(sort numerically (keys %action)){
        $formatSeq .= substr($seq, $pos, ($key-$pos));
        $formatSeq .= "$action{$key}";
        $pos=$key;
  }
  $formatSeq .= substr($seq, $pos, (length($seq)-$pos));

  return $formatSeq;
}




sub numerically {$a <=> $b}


############################################################
#                                                          #
#   subroutines used by all interfaces                     #
#                                                          #
############################################################

sub shared_parameters{

# used by the main program
# takes in the form data and creates a partial Primer3 input list
# each interface adds its own additional parameters before
# running Primer3

  my($q)= @_;
  my$PRIMER_OPT_SIZE=  $q-> param("PRIMER_OPT_SIZE") ||"20";
  my$PRIMER_MIN_SIZE=  $q-> param("PRIMER_MIN_SIZE") ||"18";
  my$PRIMER_MAX_SIZE=  $q-> param("PRIMER_MAX_SIZE") ||"23";
  my$PRIMER_OPT_TM=  $q-> param("PRIMER_OPT_TM") ||"60";
  my$PRIMER_MIN_TM=  $q-> param("PRIMER_MIN_TM") ||"55";
  my$PRIMER_MAX_TM=  $q-> param("PRIMER_MAX_TM") ||"65";
  my$PRIMER_MAX_DIFF_TM=  $q-> param("PRIMER_MAX_DIFF_TM") || "5";
  my$PRIMER_SELF_ANY=  $q-> param("PRIMER_SELF_ANY") ||"6.00";
  my$PRIMER_SELF_END=  $q-> param("PRIMER_SELF_END") ||"3.00";
  my$PRIMER_MAX_POLY_X=  $q-> param("PRIMER_MAX_POLY_X") ||"4";
  my$PRIMER_MAX_END_STABILITY= "9";
  my$PRIMER_MIN_GC=  $q-> param("PRIMER_MIN_GC") ||"30";
  my$PRIMER_GC_CLAMP=  $q-> param("PRIMER_GC_CLAMP") || "0";
  my$PRIMER_MAX_GC=  $q-> param("PRIMER_MAX_GC") ||"70";

  my$parameters = <<EOF;
PRIMER_OPT_SIZE=$PRIMER_OPT_SIZE
PRIMER_MIN_SIZE=$PRIMER_MIN_SIZE
PRIMER_MAX_SIZE=$PRIMER_MAX_SIZE
PRIMER_OPT_TM=$PRIMER_OPT_TM
PRIMER_MIN_TM=$PRIMER_MIN_TM
PRIMER_MAX_TM=$PRIMER_MAX_TM
PRIMER_MAX_DIFF_TM=$PRIMER_MAX_DIFF_TM
PRIMER_SELF_ANY=$PRIMER_SELF_ANY
PRIMER_SELF_END=$PRIMER_SELF_END
PRIMER_MAX_POLY_X=$PRIMER_MAX_POLY_X
PRIMER_MAX_END_STABILITY=$PRIMER_MAX_END_STABILITY
PRIMER_MIN_GC=$PRIMER_MIN_GC
PRIMER_MAX_GC=$PRIMER_MAX_GC
PRIMER_GC_CLAMP=$PRIMER_GC_CLAMP
=
EOF
;

  return $parameters;
}

sub get_filedata{

# subroutine for reading in files from form input

  my ($q) = @_;
  use constant UPLOAD_DIR => "/usr/local/apache/share/eur/htdocs/fgg/kgen/primer/output/";
  use constant BUFFER_SIZE =>16_384;
  use constant MAX_FILE_SIZE =>   5 * 1_048_576;
  use constant MAX_DIR_SIZE =>  10 * 1_048_576;
  use constant MAX_OPEN_TRIES =>100;
  
  $CGI::DISABLE_UPLOADS = 0;
  $CGI::POST_MAX = MAX_FILE_SIZE;
  
  $q->cgi_error and error( $q, "Error transferring file: ". $q->cgi_error );
  
  my $file = $q->param( "file" ) || error( $q, "No file received." );
  my $fh = $q->upload( "file" ) || error( $q, "upload failed: $file.");
  my $buffer = "";
  my $filedata="";
  
  if ( dir_size( UPLOAD_DIR ) + $ENV{CONTENT_LENGTH} > MAX_DIR_SIZE ) {
   error( $q, "Upload directory is full." );
  }
  
# read in contents of file
  while ( read($fh, $buffer, BUFFER_SIZE ) ) {
    $filedata="$filedata$buffer";
  }
  
# some error checking
  unless ($filedata =~ /ORIGIN/){
    error( $q, "There does not seem to be any sequence in the file")
  }
  
  $filedata =~ s/\r//g;
  return $filedata;
}

sub dir_size {

# used by sub get_filedata to list dir size

 my $dir = shift;
 my $dir_size = 0;

# loop through files and sum the sizes; doesn't descend down subdirs
 opendir DIR, $dir or die "Unable to open $dir: $!";
 while ( readdir DIR ) {
 $dir_size += -s "$dir/$_";
 }
 return $dir_size;
}

# the 4 subroutines below (get_annotation_and_dna, parse_annotation
# parse_features and revcom) are from BeginPerlBioinfo.pm,
# the library of subroutines that belongs to Beginning Perl for 
# Bioinformatics written by James Tisdall and published by 
# O'Reilly & Associates 
# (c) 2001 James Tisdall


# get_annotation_and_dna
#
#   - given filehandle to open GenBank library file, get next record

sub get_annotation_and_dna {

    my($record) = @_;

    my($annotation) = '';
    my($dna) = '';
    my($addit) = '';

    # Now separate the annotation from the sequence data
    ($annotation, $dna, $addit) = ($record =~ /^(LOCUS.*ORIGIN\s*\n)(.*)\/\/\n(.*)/s);

    # clean the sequence of any whitespace or / characters 
    #  (the / has to be written \/ in the character class, because
    #   / is a metacharacter, so it must be "escaped" with \)
    $dna =~ s/[\s\/\d]//g;

# temporary solution: genomic GenBank records do not contain any DNA anymore
# assume the user added a FASTA record to the GenBank file
    unless($dna){
      $dna = $addit;
      $dna =~ s/\>(.*?)\n//;
      $dna =~ s/[\s\/\d]//g;
    } 
    unless($dna){
      error ($q, "Inputfile does not seem to contain any sequence")
    }
    return($annotation, $dna)
}

# parse_annotation
#
#  given a GenBank annotation, returns a hash  with
#   keys: the field names
#   values: the fields

sub parse_annotation {

    my($annotation) = @_; 
    my(%results) = (  );

    while( $annotation =~ /^[A-Z].*\n(^\s.*\n)*/gm ) {
        my $value = $&;
        (my $key = $value) =~ s/^([A-Z]+).*/$1/s;
        $results{$key} = $value;
    }

    return %results;
}


# parse_features
#
#  extract the features from the FEATURES field of a GenBank record

sub parse_features {

    my($features) = @_;   # entire FEATURES field in a scalar variable

    # Declare and initialize variables
    my(@features) = ();   # used to store the individual features

    # Extract the features
    while( $features =~ /^ {5}\S.*\n(^ {21}\S.*\n)*/gm ) {

        my $feature = $&;
        push(@features, $feature);

    }

    return @features;
}

# revcom 
#
# A subroutine to compute the reverse complement of DNA sequence

sub revcom {

    my($dna) = @_;
    # First reverse the sequence
    my$revcom = reverse$dna;

    # Next, complement the sequence, dealing with upper and lower case
    # A->T, T->A, C->G, G->C
    $revcom =~ tr/ACGTacgt/TGCAtgca/;

    return $revcom;
}


sub run_primer3{

# used by all interfaces
# run_primer subroutine: copied from original
# primer3_www_results.cgi owned by the Whitehead Institute


  my($params_general) = shift @_;
  my($parameters_specific) = shift @_;
  my$PRIMER_BIN = '/srv/mblab/bin/pcr_pipeline/bin/primer3_core';
#  my$PRIMER_BIN = '/usr/local/apache/share/eur/htdocs/fgg/kgen/primer/bin/primer3_core';
  my $cline;
  my $cmd = "$PRIMER_BIN";
  my $primer3_pid;
  my @primers;
  my%results = ();

# keep taint-mode happy
  delete $ENV{PATH};
  delete $ENV{BASH_ENV};

  use IPC::Open3;
  use FileHandle;

  my$parameters = $parameters_specific.$params_general;
  my($childin, $childout) = (FileHandle->new, FileHandle->new);
  {
   local $^W = 0;
    $primer3_pid = open3($childin, $childout, $childout, $cmd);
  }

  if (!$primer3_pid) {
    print "Cannot excecute $cmd:<br>$!";
    print "wrapup\n";
    exit;    
  }

  print $childin $parameters;
  $childin->close;

  @primers = $childout->getlines;
  
  waitpid $primer3_pid, 0;
  foreach my$entry(@primers){
    my($name, $value)= split(/=/, $entry);
    chomp$value;
    $results{$name}=$value;
  }
  return %results;
}

sub print_primers{

# subroutine formats Primer3 output
# and extracts primers, primer sizes, gc content
# and melting temperature

  my(%results)= @_;
  (my$left_size = $results{PRIMER_LEFT}) =~ s/\d+\,//;
  (my$right_size = $results{PRIMER_RIGHT}) =~ s/\d+\,//;
  (my$left_seq = $results{PRIMER_LEFT_SEQUENCE}) =~ tr/actg/ACTG/;
  (my$right_seq = $results{PRIMER_RIGHT_SEQUENCE}) =~ tr/actg/ACTG/;
  my$left_repeat = find_repeat_in_primer($left_seq);
  my$right_repeat = find_repeat_in_primer($right_seq);
  my$left_warning = "";
  my$right_warning = "";
  if ($left_repeat){
    $left_warning = '<b>It looks like this primer contains a '."$left_repeat".' repeat!</b><br>';
  }
  if ($right_repeat){
    $right_warning = '<b>It looks like this primer contains a '."$right_repeat".' repeat!</b><br>';
  }

  my($info) = <<STOP
LEFT = $left_seq<br>
SIZE = $left_size<br>
TM = $results{PRIMER_LEFT_TM}<br>
GC% = $results{PRIMER_LEFT_GC_PERCENT}<br>
$left_warning
<br>
RIGHT = $right_seq<br>
SIZE = $right_size<br>
TM = $results{PRIMER_RIGHT_TM}<br>
GC% = $results{PRIMER_RIGHT_GC_PERCENT}<br>
$right_warning
<br>
PRODUCTSIZE = $results{PRIMER_PRODUCT_SIZE}<br>
STOP
;

  return $info;
}

sub find_repeat_in_primer{

# find_repeat_in_primer is called from print_primers
# It creates 6nt subsets of the inputstring
# and checks them for perfect repeats using
# find_repeat

  my($seq) = @_;
  for(my$i=0; $i < (length($seq) - 6); $i++){
    my$out = find_repeat(substr($seq, $i, 6));

# two copies of a sequence don't make a repeat
    if( $out ){
      unless($seq =~ /($out){3}/){
        $out="";
      }
    }
# one repeat is enough
    if($out){
      return $out;
    }
  }
}


sub find_repeat{

# find_repeat called from find_repeat_in_primer
# looks in subsets for perfect repeats
  my($string) = @_;
  my$count=length$string;
  my@ray=split('', $string);

  my$pos="-1";
  my$pos_string="-1";
  for(my$i=0; $i<$count; $i++){
    next if($i==0);
    my$nt=$ray[$i];
    my$flag;
    for(my$k=$pos+1; $k>-1; $k--){
      if($nt eq $ray[$k]){
        $pos=$k;
        $flag="1";
        last;  
      }
    }
    unless($flag){
      $pos="-1";
    }
    $pos_string.=$pos;
  }
# to get the unit of repeat, take the last $pos and remove that much from the string
# if it is a true repeat
  $pos++;
  my$replength=$count-$pos;
  if( ($count ne $replength) && (($count-$replength) % $replength) =="0"){
    my$repunit=substr($string, 0, $replength);
    return $repunit;
  }
}



############################################################
#                                                          #
#     subroutines used for genomic primers                 #
#                                                          #
############################################################

sub find_genomic_primers{

# called by main program
# create primers around each exon of a gene
# generate primer3 input parsed data:
# run primer3 with current parameters

  my($params_general) = shift @_;
  my($FLANKING_SIZE)=shift @_;
  my($PRIMER_PRODUCT_SIZE_RANGE) = shift @_;
  my(@exons) = @_;
  my$nr = "0";
  my$exon = "";
  my$exonsize;
  my($max) = $PRIMER_PRODUCT_SIZE_RANGE =~ /-(\d+)/;

# create cDNA from @exons 
  my$cDNA = join ("", @exons);
  $cDNA =~ s/[actgn]//g;
  $cDNA = format_sequence("0", "0", "60", "0", "0", "0", $cDNA);

# create primers for each exon
# generate primer3 input parsed data:
# run primer3 with current parameters

  foreach $exon(@exons){
    $nr++;
    $exonsize = length($exon) - 2*$max;
    push (@allout,  "<b>exon nr $nr</b>");
    if ($exonsize >$max){
      push (@allout, "exon $nr is too big; no primers designed (exon in blue)");

      $exon = format_sequence("0", "0", "60", "0", "0", "0", $exon);
      push (@allout, "$exon");
      next;
    }
    my$from = $max+1- $FLANKING_SIZE;
    my$size = $exonsize+ 2* $FLANKING_SIZE; 

    my$parameters_specific = <<ENDOF;
PRIMER_SEQUENCE_ID=exon $nr
PRIMER_PRODUCT_SIZE_RANGE=$PRIMER_PRODUCT_SIZE_RANGE
TARGET=$from,$size
SEQUENCE=$exon
PRIMER_NUM_RETURN=1
ENDOF
;
    my%results = run_primer3($params_general, $parameters_specific);
  
    unless(exists($results{PRIMER_LEFT})){
      push (@allout, "no primers found for exon $nr (exon in blue)");
      $exon = format_sequence("0", "0", "60", "0", "0", "0", $exon);
      push (@allout, "$exon");
      next;
    }
  
# send the primer info to the print subroutine
    my$onfo =  print_primers(%results);
    push (@allout, $onfo, " ");

# format output sequence
    $exon = format_sequence("genomic", \%results, "60", "blue", "red", "A", "0");
    push (@allout, $exon, "&nbsp");
  }
  push (@allout, "<hr>",'<a name="cDNA">and the cDNA is:</a>',
         "$cDNA");
  return @allout;
}


sub select_exons{

# used by main program, genomic primers loop
# runs subroutines to extract exons
# from GenBank file

  my$gbseq = shift @_;
  my$gene = shift@_;
  my$PRIMER_PRODUCT_SIZE_RANGE = shift@_;
  my$annotation;
  my$end_right;
  my$exon;
  my$exonnr;
  my$dna;
  my$length_left;
  my$length_right;
  my$start_left;
  my$start_right;
  my@exons = ();
  my@features = ();
  my@fields = ();
  my@primerseqs = ();
  my@upcase = ();
  my%fields = ();

  push (@allout, '<a href="#cDNA">Go to cDNA</a>');

# get the exon positions from the genbankfile

  ($annotation, $dna) = get_annotation_and_dna($gbseq);

  %fields = parse_annotation($annotation);

  @features = parse_features($fields{FEATURES});

  @exons = find_exons($gene, @features);

# the last entry is the orientation
  my$orientation = pop@exons;
  push (@allout, "The orientation of the gene was $orientation");

# the second last entry of @exons is a warning

  my$warning = pop@exons;
  push (@allout, "$warning");

# select the exon sequences and flanking sequence
# (exons in uppercase)

  @upcase = select_exons_and_flanking($dna,
    $PRIMER_PRODUCT_SIZE_RANGE, $orientation, @exons);

  return (@upcase);
}


sub find_exons{

# used by  sub select_exons
# find_exons: subroutine to select the begin and end location of
# each exon of each mRNA given for a gene.
# input is a gene name as it occurs in the GenBank file and an array 
# containing the features table of a genbank file as created by
# the parse_features subroutine
# output is an array with numbers separated by two periods
# this subroutine uses the extract_numbers subroutine

  my($gene) = shift(@_);
  my(@gbfeatures) = @_;
  my @allnumbers = ();
  my $flag = "";
  my $outsize= "";
  my $orientation = "normal";

  foreach my $feature(@gbfeatures){
    my($featurename) = ($feature =~ /(\S+)/);
    if ($flag && $featurename eq "mRNA" && $feature =~ /$gene/){
      if ($feature =~ /\(</){
        $outsize = "Warning! The sequence does not contain the start of the gene, skipped the first exon(s)";
      } elsif ($feature =~ /\.>/){
        $outsize = "Warning! The sequence does not contain the end of the gene, skipped the last exon(s)";
      }
      if($feature =~ /complement/){
        $orientation = "reverse";
      }
      my@numbers = extract_numbers ($feature);
      @allnumbers = (@allnumbers, @numbers);
    }elsif ($featurename =~ /gene/){
# genename sometimes flanked by quotes
      my($genename) = ($feature =~ /\/gene=\"{0,1}(.*)\"{0,1}/);
      $genename=~ s/\"//;
      if ($genename eq $gene){
        $flag = 1;
      }
    }
  }
  push (@allnumbers, $outsize, $orientation);
  return @allnumbers;
}


sub extract_numbers{

# used by sub find_exons
# gets exon start and end from GenBank feature

  my($feature) = @_;

  $feature =~ s/\s//g;

  my($numbers) = ($feature =~ /join\((.*?)\)/);
  my@numbers = split(",", $numbers);
  return @numbers;
}

sub select_exons_and_flanking{

# used by sub select_exons
# finds exons in genomic sequence,
# selects the exon with 200 bp flanking on each side
# and puts exon sequence in uppercase

  my($dna) = shift(@_);

  my($PRIMER_PRODUCT_SIZE_RANGE) = shift@_;
  my($orientation) = shift@_;
  my(@exons) = @_;
  my @exon_plus_flanking = ();
  my ($max) = $PRIMER_PRODUCT_SIZE_RANGE =~ /-(\d+)/;

  foreach my$numberset(@exons){
    my($begin, $end) = split (/\.\./, $numberset);
    my$length = $end - $begin +1;
    my$exon = substr($dna, $begin-1-$max, $length+2*$max);
    $exon =~ tr/ACTG/actg/;
    if(length $exon < ($max + $length)){
        error($q, "It looks like you did not use enough flanking sequence. There should be at least $max nt before the
          first and after the last exon");
    }
    substr($exon, $max, $length) =~ tr/actg/ACTG/;
    push(@exon_plus_flanking, $exon);
  }
  if($orientation eq "reverse"){
    my@rev = ();
     my$exn;
     while(@exon_plus_flanking){
       $exn = pop@exon_plus_flanking;
       $exn = revcom($exn);
       push (@rev, $exn);
     }
     @exon_plus_flanking = @rev;
  }

  return @exon_plus_flanking;
}

############################################################
#                                                          #
#   subroutines used for SNP primers                       #
#                                                          #
############################################################

sub find_snp_primers{

# called by main program
# create primers for each SNP
# generate primer3 input parsed data:
# run primer3 with current parameters

  my($params_general) = shift @_;
  my($PRIMER_PRODUCT_SIZE_RANGE) = shift @_;
  my(@snps) = @_;
  my$sz = scalar@snps;
  my$i = "0";
  my ($max) = $PRIMER_PRODUCT_SIZE_RANGE =~ /-(\d+)/;
  $max -= 5;
  my$excl;

  while ($i<$sz){
    my ($seq, $id, $letter, $location)  = 
      ($snps[$i], $snps[$i+1],$snps[$i+2],$snps[$i+3]);
    $i += 4;

# the crack_code subroutine translates Y into C/T etc

    my$letters = crack_code($letter);
    push (@allout,  "<b>$id = $letters ($letter)<\/b>", "position in sequence: $location");

# Primer3 does not recognise S, W, R etc...
# but don't want primer in SNP, so exclude these regions for primers

    my$nr_of_snps = scalar(my@nr = $seq =~ /[SWMKRY]/g);
    my@pos_SNP = ();
    while($seq =~ m/[SWMKRY]/g){
      push(@pos_SNP, pos($seq)); 
    }
    if (scalar@pos_SNP>1){
      $excl = join (",1 ", @pos_SNP);
      $excl .= ",1";
    } 

# Input must be ACG or T, but will be masked in Primer3

    (my$dna = $seq) =~ s/[WRCMKS]/T/g;

    my$parameters_specific = <<ENDOF;
PRIMER_SEQUENCE_ID=$id
PRIMER_PRODUCT_SIZE_RANGE=$PRIMER_PRODUCT_SIZE_RANGE
TARGET=$max,10
SEQUENCE=$dna
PRIMER_NUM_RETURN=1
ENDOF
;

    if ($excl){
      $parameters_specific.="EXCLUDED_REGION=$excl\n";
    }

# run Primer3 with current parameters

    my%results = run_primer3($params_general, $parameters_specific);
  
# see if there are primers in the Primer3 output

    unless(exists($results{PRIMER_LEFT})){
      push (@allout, "<b>no primers found</b>");
      $seq = format_sequence("0", "0", "60", "0", "0", "0", $seq);
      $seq =~ s/([KMYRSW])/<font color="blue">$1<\/font>/g;
      push (@allout, "$seq", "&nbsp");
      next;
    }
  
# send the primer info to the print subroutine
    my$onfo =  print_primers(%results);
    push (@allout,$onfo, " ");

# format the output sequence with the primers

    $seq = format_sequence("snp", \%results, 60, "blue", "black", "A", $seq);
# color the SNP
    $seq =~ s/([KMYRSW])/<font color="red">$1<\/font>/g;
    push (@allout, $seq, "&nbsp");

  }
  return @allout;
}


sub select_range{

# subroutine used by main program to select SNPs and their
# flanking sequence.

  my$gbseq = shift @_;
  my$PRIMER_PRODUCT_SIZE_RANGE = shift@_;
  my$end_right;
  my$exon;
  my$exonnr;
  my$length_left;
  my$length_right;
  my$start_left;
  my$start_right;
  my@fields = ();
  my@primerseqs = ();

# get the SNP information from the genbankfile

  my($annotation, $dna) = get_annotation_and_dna($gbseq);
  my%fields = parse_annotation($annotation);
  my@features = parse_features($fields{FEATURES});

# the find_variation subroutine selects only SNPs that
# are not marked 'ambiguous'

  my(@snps) = find_variation(@features);

# use the SNP information to get the sequence

  my@snpseq = select_SNPs_and_flanking($dna,
    $PRIMER_PRODUCT_SIZE_RANGE, @snps);

  return (@snpseq);
}

sub select_SNPs_and_flanking{

# subroutine used by sub select_range to
# get every SNP and flanking sequence

  my($dna) = shift(@_);
  my($PRIMER_PRODUCT_SIZE_RANGE) = shift@_;
  my(@snps) = @_;
  my @snp_plus_flanking = ();
  my ($max) = $PRIMER_PRODUCT_SIZE_RANGE =~ /-(\d+)/;
  my $sz = scalar@snps;
  my $i = "0";

# every 1st, 5th, 10th, etc entry of the list is a 
# SNP position

  while($i<$sz){
    my$number = $snps[$i];
    (my$al1 = $snps[$i+1]) =~ tr/[a-z]/[A-Z]/;
    (my$al2 = $snps[$i+2]) =~ tr/[a-z]/[A-Z]/;
    my$id = $snps[$i+3];

# find_code turns C/T into Y etc
    my$code = find_code($al1, $al2);
    substr($dna, $number-1, 1, "$code");
    my$startpoint = $number-1-$max;
    if($startpoint <0){
      $startpoint ="0";
    }
    my$snp = substr($dna, $startpoint, 2*$max);
    push(@snp_plus_flanking, $snp, $id, $code, $number-1);
    $i+=4;
  }
  return @snp_plus_flanking;
}

sub find_code{

# used by sub select_SNPs_and_flanking
# find_code is the reverse of crack_code:
# it inserts a single letter for every SNP:
# A/G becomes R etc.

  my$al1 = shift @_;
  my$al2 = shift @_;
  my$letter;

  if (($al1 eq "C" && $al2 eq "G")||($al1 eq "G" && $al2 eq "C")){
    $letter = "S";
  }
  elsif (($al1 eq "A" && $al2 eq "T")||($al1 eq "T" && $al2 eq "A")){
    $letter = "W";
  }
  elsif (($al1 eq "C" && $al2 eq "T")||($al1 eq "T" && $al2 eq "C")){
    $letter = "Y";
  }
  elsif (($al1 eq "A" && $al2 eq "G")||($al1 eq "G" && $al2 eq "A")){
    $letter = "R";
  }
  elsif (($al1 eq "A" && $al2 eq "C")||($al1 eq "C" && $al2 eq "A")){
    $letter = "M";
  }
  elsif (($al1 eq "G" && $al2 eq "T")||($al1 eq "T" && $al2 eq "G")){
    $letter = "K";
  }
  else{ $letter = "unknown";}
  return $letter;
}

sub crack_code{

# used by sub find_snp_primers
# reverse of find_code. Changes S into G/C etc.

  my($letter) = @_;
  my$alleles;

  if($letter eq "S"){
    $alleles = "G/C";
  }elsif($letter eq "W"){ 
    $alleles = "A/T";
  }elsif($letter eq "R"){ 
    $alleles = "G/A";
  }elsif($letter eq "Y"){ 
    $alleles = "T/C";
  }elsif($letter eq "K"){ 
    $alleles = "G/T";
  }elsif($letter eq "M"){ 
    $alleles = "A/C";
  }
  return $alleles;
}


sub find_variation{

# used by sub select_range
# takes a list of GenBank 'features' and select the SNPs
# skips the SNPs marked 'ambiguous'
# parses out the positions, the alleles and the position
# of the SNP

  my(@gbfeatures) = @_;
  my @allnumbers = ();
  foreach my $feature(@gbfeatures){
    if($feature =~ /ambiguous/){next}
    if($feature =~ 
      /variation\s*(\d+).*?(allele|replace)=\"(\w).*?(allele|replace)=\"(\w).*?(dbSNP\:\d+)\"/s){
      push(@allnumbers, $1, $3, $5, $6);
    }
  }
  return @allnumbers;
}

############################################################
#                                                          #
#     subroutines used for overlap primers                 #
#                                                          #
############################################################


sub find_overlap_primers{

# called by main program
# create primers in overlapping sets
# generate primer3 input parsed data:
# run primer3 with current parameters
# until complete target is covered

  my($params_general) = shift @_;
  my($PRIMER_PRODUCT_SIZE_RANGE) = shift @_;
  my($OVERLAP_SIZE_RANGE) = shift @_;
  my($TARGET) = shift @_;
  my($seq) = shift @_;
  my$nr = "1";
  my($from, $size) = split /,/, $TARGET;
  my($prod_size_lo, $prod_size_hi) = split/-/, $PRIMER_PRODUCT_SIZE_RANGE;
  my($overlap_lo, $overlap_hi) = split/-/, $OVERLAP_SIZE_RANGE;
  my$overlap = $overlap_hi - $overlap_lo;
  my@primers = ();
  my%results;

  my$outputseq = $seq;
  my$whole_target = $TARGET;

# remove 3' sequence
  my$cutoff = $from + $size + $prod_size_hi;

# keep track of end
  my$end = $from + $size;
  $seq = substr($seq, 0, $cutoff);

# the target of the first product is the given target

  my($new_start) = ($TARGET =~ /(^\d*)/); 

# the length of the target is $prod_size_lo

  $TARGET=$new_start.','.$prod_size_lo;

  my$parameters_specific;
  my$right;
  my$onfo;
  my$k;
  my$targetseq;
  my$noprimers_seq;
  my@matched_seqs = ();
  my$startpos = "0";
  my$endpos = "0";
  my$orisize = length($seq); 

# design primers while target is still part of $seq
  while ($end>0){
  print $q->p(" ");
    $parameters_specific = <<ENDOF;
PRIMER_SEQUENCE_ID=PCR product $nr
PRIMER_PRODUCT_SIZE_RANGE=$PRIMER_PRODUCT_SIZE_RANGE
PRIMER_NUM_RETURN=1
SEQUENCE=$seq
TARGET=$TARGET
ENDOF
;

    %results = run_primer3 ($params_general, $parameters_specific); 

# if there are no primers, cut off a 5' piece of the sequence
# (size is $prod_size_lo), keep track of $end and try again

    unless(exists($results{PRIMER_LEFT})){
      push (@allout, "<b>no primers found for target $nr</b>", "&nbsp");
      substr($seq, 0, $prod_size_lo,"");
      $end -= $prod_size_lo;
      $TARGET=$overlap.','.$prod_size_lo;
      $nr++;
      next;
    }
  
    push (@primers, $results{PRIMER_LEFT_SEQUENCE});
    $right = revcom($results{PRIMER_RIGHT_SEQUENCE});
    (my$pos) = $results{PRIMER_LEFT} =~ /(\d+)\,/;
    $startpos = $orisize - length($seq) + $pos +1;
    $endpos = $results{PRIMER_PRODUCT_SIZE};
    $endpos += $startpos - 1;
    push (@allout, "<b>target $nr</b> position in inputsequence: $startpos to $endpos");
    push (@primers, $right);

# send the primer info to the print subroutine
    $onfo =  print_primers(%results);
    push (@allout,$onfo);

# keep the matched seq (for later output)
    $targetseq = cutoff_product($seq, $results{PRIMER_LEFT}, $results{PRIMER_RIGHT});
    push(@matched_seqs, $targetseq);

# format the sequence that is covered (trim ends outside the product)

    my$printseq = format_sequence("overlap", \%results, "60", "blue", "black", "0", "0");
    push (@allout, $printseq, "&nbsp");

# cut off 5' end, keep track of $end 
    ($new_start) = ($results{PRIMER_RIGHT} =~ /(^\d*)/); 
    $new_start -= $overlap_hi;
    substr($seq, 0, $new_start,"");
    $end -= $new_start;
    $nr++;
    $TARGET=$overlap.','.$prod_size_lo;
  }

# End of loop. Format inputsequence to show which part is covered
# Color matched sequence in input and format 

  $outputseq = format_covered_cdna(@matched_seqs, "60", $whole_target, $outputseq);

  push (@allout, "<hr>","<b>your result</b> (TARGET IN UPPERCASE, 
<font color=\"blue\">covered 
regions are blue</font>. Colored regions include primers!):", 
$outputseq);

  return @allout;
}


sub cutoff_product{

# used by sub find_overlap_primers
# removes all sequence outside PCR product

  my($seq) = shift(@_);
  my($left) = shift(@_);
  my($right) = shift(@_);
      
  my($start_left, $length_left) = 
     split (/\,/, $left);
  my($end_right, $length_right) = 
     split (/\,/, $right);
  my $start_right = $end_right+1  - $length_right;

  $seq = substr($seq, 0, $end_right+1);
  substr($seq, 0, $start_left,"");   
  return uc($seq);
}


sub clean_sequence{

# used by main program in overlap_primers loop
# removes FASTA header line, spaces and numbers
# from inpusequence

  my($q) = @_;
  my$seq=$q->param("SEQUENCE") || error($q, "please fill in sequence!");
  if ($seq =~ />/){
    $seq =~ s/^>.*//;
  }
  if ($seq =~ /\d/){
    push (@allout, "Warning: numbers in input sequence were deleted!");
    $seq =~ s/\d//g;
  }
  $seq =~ s/\s//g;
  return $seq;
}



############################################################
#                                                          #
#     subroutines used for cDNA primers                    #
#                                                          #
############################################################


sub find_cDNA_primers{

# called by main program
# create primers around ORF for cDNA in each GenBank file
# if the ORF is the same size as the cDNA
# the PCR product is maximised.
# generate primer3 input parsed data:
# run primer3 with current parameters

  my($params_general) = shift @_;
  my($PRIMER_PRODUCT_SIZE_RANGE) = shift @_;
  my($from) = shift @_;
  my($to) = shift @_;
  my($seq) = shift @_;
  my$dnasize = length($seq);
  my$penalty = "0.05";
  my $size = $to+1 - $from;
  my($min, $max) = split /-/, $PRIMER_PRODUCT_SIZE_RANGE;
  my$parameters_specific;
  my%results = ();

  if ($size > $max){
    push (@allout, 
      "ORF size $size is larger than maximum size ($max), use overlappcr or increase maximum product size!", "Target in blue.");
    $seq = upcase_DNA($seq, $from-1, $size);
    $seq = format_sequence("0", "0", "60", "0", "0", "0", $seq);
    push (@allout, $seq);
    return @allout;

  } elsif ($size == $dnasize){
    push (@allout, 
      "ORF is same size as cDNA, creating maximum product possible");
    $parameters_specific = <<ENDOF;
PRIMER_PRODUCT_SIZE_RANGE=$PRIMER_PRODUCT_SIZE_RANGE
PRIMER_PRODUCT_OPT_SIZE=$size
PRIMER_PAIR_WT_PRODUCT_SIZE_GT=$penalty
PRIMER_PAIR_WT_PRODUCT_SIZE_LT=$penalty
ENDOF
;
  }else{
    $seq = upcase_DNA($seq, $from-1, $size);
    push (@allout, "ORF in uppercase. Size is $size.");
    $parameters_specific = <<ENDIF;
PRIMER_PRODUCT_SIZE_RANGE=$PRIMER_PRODUCT_SIZE_RANGE
TARGET=$from,$size
ENDIF
;
  }

  my $params_cDNA = <<END;
PRIMER_SEQUENCE_ID=komtnog
SEQUENCE=$seq
PRIMER_NUM_RETURN=1
END
;
  $params_cDNA .= $parameters_specific; 

  %results = run_primer3 ($params_general, $params_cDNA); 
  unless(exists($results{PRIMER_LEFT})){
    push (@allout, "no primers found (target in blue)");
    $seq = format_sequence("cdna", "0", "60", "0", "0", "0", $seq);
    push (@allout, "$seq");
    return @allout;
  }
  
# send the primer info to the print subroutine
  my$onfo =  print_primers(%results);
  push (@allout,$onfo, " ");


# format sequence for output
  $seq = format_sequence("cdna", \%results, "60", "blue", "red", "A", "0");

  push (@allout, $seq);
  return @allout;

}

sub upcase_DNA{

# used by sub find_cDNA primers
# puts sequence in uppercase letters
# based on start and end positions

  my($seq) = shift(@_);
  my($from) = shift(@_);
  my($size) = shift(@_);

  substr($seq, $from, $size) =~ tr/actg/ACTG/;
  return $seq;
}

sub get_cDNA_ID_and_sequence{

# used by main program, cDNA loop
# subroutine gets the cDNA ID, sequence and ORF position
# from a GenBank record

  my($record) = (@_); 
  my($annotation, $dna) = get_annotation_and_dna($record);
  my%feat;
  my %fields = parse_annotation($annotation);

# Extract the features from the FEATURES table
  my @features = parse_features($fields{FEATURES});

  foreach my$feature(@features){
# extract the name of the feature
    my($featurename) = ($feature =~ /^ {5}(\S+)/);
    $feat{$featurename} = $feature;
  }

# this needs some parsing to give the gene name only:
  my$gene;
  ($gene) = $feat{gene} =~ m/(gene=.*\n)/;
  $gene =~ s/gene="//;
  $gene =~ s/"\n//;

# extract the coding nucleotides from the CDS feature 

  my($coding_nts_from, $coding_nts_to) = ($feat{CDS}=~/(\d*)\.{2}(\d*)/);

  if($coding_nts_from eq "" | $coding_nts_to eq ""){
    print "WARNING: at least one gene $gene does not contain proper CDS 
feature and is not used for analysis\n"; 
    next;
  }
  return ($dna, $gene, $coding_nts_from, $coding_nts_to);
}



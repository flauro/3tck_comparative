#!/usr/bin/perl -w
# Federico Lauro
# This script is used to analyze the percentage identity between two sequences and the average nucleotide identity
# as described in Goris et al., IJSEM 57:81-91
# It is based on a quick parser (blast_parser_m8.pl)
# to convert a standard blast output to a <TAB> delimited file format
# With the following colums:
# Query_ID <TAB> Query_length <TAB> Hit_ID <TAB> HIT_description <TAB> evalue <TAB> score <TAB> Bitscore <TAB>
# Percent_identity <TAB> Best_HSP_start (posn. on query sequence) <TAB> Best_HSP_end (posn. on query sequence)


use strict;
use warnings;
use vars qw($USAGE);
use Getopt::Long;
use Bio::SearchIO;
use Bio::SeqIO;

my $counter=0;
my $seq_counter=0;
my $dna;
my $ANI_counter=0;
my $query_length=0;	# Stores the total query length
my $start;		# Stores the coord of the hit start
my $end;		# Stores the coord of the hit end
my $PCD=0;		# Stores the total length of the percentage conserved DNA
my $cutoff;		# Stores the 70% of the fragment length used as a cutoff for computing the ANI
my $aligned;		# Store the difference end-start
my $aligned_len;
my $rec_id;
my $ANI_sum=0;		# Store the total of all the $rec_id of the fragments that get over the 30% id

$USAGE = "\nANI.pl [-h] -q query -r reference -b blast_run -s size (default=1020 nt) -n undetermined -v verbose";

my ($query_file, $reference_file, $help, $blast_run, $undet, $size, $verbose) = (undef, undef, undef, undef, undef, '1020', '1');

&GetOptions('query|q=s'		     => \$query_file,
	    'reference|r=s'          => \$reference_file,
	    'blast_run|b'	     => \$blast_run,
	    'help|h'                 => \$help,
	    'undetermined|n'	     => \$undet,
	    'size|s=s'		     => \$size,
	    'verbose|v=s'            => \$verbose
	    );


if( $help ) {
    exec('perldoc', $0);
    die;
}

if( !defined $query_file ) {
    die($USAGE . "\n\tMust specify a query FASTA file\n\n");
}
if( !defined $reference_file ) {
    die($USAGE . "\n\tMust specify a reference FASTA file\n\n");
}

# Here i split the query genome in pieces $size nucleotides long
# and store them in a file called temp_query

$cutoff=$size*0.7;
open (TEMP, ">query_temp");

my $in = Bio::SeqIO->new(-file => "$query_file" , -format => 'FASTA');


	while ( my $seq = $in->next_seq()) {
                                        $dna=$seq->seq();
					$dna=~s/N//g if ( $undet );    # Removes all N's from a sequence if the -n flag is on
					$query_length += length ($dna);
                                        for (my $q=0; $q < (length ($dna) - $size+1); $q+=$size)
						{
						print TEMP ">",$seq->id,'_',$seq_counter,"\n";
						print TEMP substr($dna,$q,$size),"\n";
						++$seq_counter;
						}
					}


close TEMP;
print "\n\nSplit the query multiFASTA file into ", $seq_counter, " sequence fragments\n" if ($verbose > 0);

# Let the blastn begin
if ( !defined $blast_run && $verbose > 0) {system 'formatdb -p F -i '.$reference_file};
if ( !defined $blast_run) {system 'formatdb -p F -i '.$reference_file};
if ( !defined $blast_run) {system 'blastall -a 2 -p blastn -d '.$reference_file.' -o blast_temp -i query_temp -X 150 -q -1 -F F'};

    # Let the parsing of the BLAST results begin
    my $record = Bio::SearchIO->new(-format =>'blast',
                                    -file =>'blast_temp',
                                    );

X:                while (my $result = $record->next_result)
                         {
			# First reset or increase the counters
			  ++$counter;
			   while (my $hit = $result->next_hit) {

					# Check for hits and eventually exit the loop
					unless ($hit) {next X;}

		                          my $hsp = $hit->next_hsp;

 						unless ($hsp) {next X;}
		
						# Only parse the selected $max_hits number of hits
						$start=$hsp->query->start;
						$end=$hsp->query->end;
						$aligned=$end-$start+1;
						$aligned_len=$hsp->length;
						$rec_id=($hsp->percent_identity)*$aligned/$size;
						if ($hsp->percent_identity > 90) {$PCD += $aligned};
						if ($aligned > $cutoff && $rec_id > 30) {++$ANI_counter;
											$ANI_sum += $hsp->percent_identity}
						print  "\rAnalyzing ",$result->query_name;
						# print $rec_id,"\t";
						# print $aligned, "\t";
						# print $result->hit_name;
						# print $aligned_len,"\n"; 
		      				next X;

								}


		 	}


print STDERR "\n\nParsed ", $counter, " input sequences\n\n" if ($verbose > 0);

print "\n\nThe query length is ",$query_length,"\n";
print "\n\nThe Percentage conserved DNA between these 2 genomes is ",$PCD/$query_length*100,"\n";
print "\n\nThe ANI for these 2 sequences is ",$ANI_sum/$ANI_counter,"\n";


# Don't forget to remove query_temp sequence file; keep blast_temp as you might want to use it in subsequent runs with options -b n

exit;


__END__



#
# ANI
#
# Federico Lauro

# POD documentation - main docs before the code

=head1 NAME

ANI.pl - Script for calculating the ANI and percent conserved DNA between 2 genomes

=head1 SYNOPSIS

% ANI.pl [-h] -q query -r reference -b blast_run -s size (default=1020 nt) -v verbose

=head1 DESCRIPTION

This script will calculate the Percentage conserved DNA and Average Nucleotide Identity
between two genomes as described in Goris et al., IJSEM 57:81-91
Two genomes belong to the same species when ANI>95% and PCD>69% which
corresponds to a DNA-DNA hybridization (DDH) of ~ 70%.
Be sure to run the script both ways (with each genome as reference separately)
as the results may vary.

=head1 OPTIONS

-query (-q) multiFASTA of the query genome

-reference (-r) multiFASTA of the reference genome

-blast_run (-b) with the -b option the blastn portion of the script is skipped (requires a file named blast_temp)

-size (-s) Size of the fragments for splitting the query genome (default=1020nt)

-undetermined (-n) Removes all the N's from the query sequence before splitting

-help (-h) displays this help

-vervose (-v) amount of information displayed (0 = no messages)

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org               - General discussion
  http://bioperl.org/MailList.shtml   - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Federico Lauro

Email flauro@unsw.edu.au

=cut

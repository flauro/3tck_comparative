#!/usr/bin/perl
# Federico Lauro
# The script will take an input file of unassebled sequences and attempt to BLAST them on a scaffold and order
# them into a pseudomolecule WITH THE CONTIGS/SCAFFOLDS SEPARATED BY THE SPACER NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN
# To run faster this script uses MEGABLAST, but can be modified to use BLAST
# Unassembled queries will be binned into the file "unassembled.txt"
# Make sure that BLAST is installed and in your path or the script will not work
# Make sure you have the latest version of bioperl to parse the MEGABLAST results
# NOTE: THE FASTA PARSER IS HACKED QUICKLY SO IT WILL ONLY WORK WITH FASTA FILES FORMATTED WITH ONE SEQUENCE
# PER LINE



use Bio::SearchIO;

my %hits;
my %nH;
my $counter=0;
my $counter2=0;
my $ndx=0;
my %IofH;
my @sorted;
my @query_sequence;
my $flag;
my $matches;

print "Name of the file containing the fragments to order","\n";
$fragment_file=<stdin>;
chomp $fragment_file;
print "Name of the file containing the sequence scaffold","\n";
$scaffold_file=<stdin>;
chomp $scaffold_file;

# This prepares the scaffold database for blasting
system "formatdb","-i","$scaffold_file","-p","F","-o","T";

# This portion script explodes the  multiFASTA $fragment_file 
# into an array of multiple single FASTA files containing each one sequence

open (MULTIF, "$fragment_file");
while ($fasta_header=<MULTIF>)
        {push (@headers, $fasta_header);
        $sequence=<MULTIF>;
        push (@query_sequence, $sequence);
        ++$counter;
        }
close MULTIF;

print "\n", "COMPUTING", "\n";

# This part runs the searches of the fragments against the scaffold and saves the results
# In a new array

foreach $in_seq (@query_sequence) {
	open (TARGET, ">query.txt");
	print TARGET $headers[$counter2];
	print TARGET $in_seq;
	close TARGET;

	# I first store the query sequence in the array of hashes
	$hits[$counter2]{'query_seq'}=$in_seq;

	# Here the magic of blasting the fragments is performed and the results are saved in temp.txt
	system "megablast","-i","query.txt","-o","temp.txt","-d","$scaffold_file","-D","2","-v","1";

	# Let the parsing of the BLAST results begin
	my $record = Bio::SearchIO->new(-format => 'blast',
					-file =>'temp.txt');

		while (my $result = $record->next_result)
  			{my $hit = $result->next_hit;

			unless ($hit) {
				$hits[$counter2]{'match_query'}=undef;
				$hits[$counter2]{'hit_begin'}=undef;
				$hits[$counter2]{'query_strand'}=undef;
				$hits[$counter2]{'hit_strand'}=undef;
				last;}

				while (my $hsp = $hit->next_hsp)

					{$matches=$matches+abs ($hsp->query->end-$hsp->query->start);
					if ($flag=0) {
							$hits[$counter2]{'hit_begin'}=$hsp->hit->start;
							$hits[$counter2]{'query_strand'}=$hsp->query->strand;
							$hits[$counter2]{'hit_strand'}=$hsp->hit->strand;
							}
					++$flag;}

				# Here we store the information for each hit into the
				# array of hashes if over the cutoff
				if ($matches/($result->query_length) > 0.8)
				{
				$hits[$counter2]{'match_query'}=$result->query_name;
				}
			
			else {	$hits[$counter2]{'match_query'}='NONE';
				$hits[$counter2]{'hit_begin'}=undef;
				$hits[$counter2]{'query_strand'}=undef;
				$hits[$counter2]{'hit_strand'}=undef;
				}
		}
$flag=0;
$matches=0;
++$counter2;
print ".";
}

print "\n\n";

# Just a safety check
if ($counter ne $counter2) {print "Fede's code sucks", "\n";
				exit;}

# Let's start writing into the output files

open (UNASSEMBLED,">unassembled.txt");
open (PSEUDO, ">pseudomolecule.txt");
open (ORDER, ">contig_order.txt");

# This creates the FASTA header for the pseudomolecule
print PSEUDO "> Pseudomolecule of ",$counter," pieces from ",$hits[0]{'match_query'},"\n";

for (my $i = 0; $i < $counter; $i++)
{
# This part removes the contigs with no significant match to the scaffold
# These fragments will go in the file unassembled.txt

if ($hits[$i]{'match_query'} eq 'NONE') {print UNASSEMBLED $headers[$i];
				print UNASSEMBLED $query_sequence[$i];
				if ($hits[$i]{'query_seq'} ne $query_sequence[$i]) {
					print "Indexing problem: Fede's code REALLY sucks";
					exit;}
				delete $hits[$i]{'query_seq'};
				}
print $hits[$i]{'match_query'}," ";
print $hits[$i]{'hit_begin'}," ";
print $hits[$i]{'query_strand'}," ";
print $hits[$i]{'hit_strand'}," ";
print substr($hits[$i]{'query_seq'},0,5);
print "\n";
}

# This part reorders the array of hashes based on the hit coordinates
# and creates a pseudomolecule with the ordered sequences
# reverse-complemented where necessary
# It also stores the contig order in the file contig_order.txt
# Complemented contigs will have a (c) appended to the name 

foreach (@hits) {
$IofH{$_->{'hit_begin'}} = $ndx++;
}
foreach (sort {$a <=> $b} keys %IofH)
{
%nH = %{$hits[$IofH{$_}]};

if ($nH{'query_strand'} eq $nH{'hit_strand'})
	{chomp $nH{'query_seq'};
	print PSEUDO $nH{'query_seq'}, "NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN";
	print ORDER $nH{'match_query'}, "..";
	}
elsif ($nH{'query_strand'} ne $nH{'hit_strand'})
        {chomp $nH{'query_seq'};
	print PSEUDO revcom($nH{'query_seq'}), "NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN";
	print ORDER $nH{'match_query'}, "(c)", "..";
        }
}

close ORDER;
close PSEUDO;
close UNASSEMBLED;
exit;





##############################################
#                                            #
#  Subroutine for reverse-complementing      #
#                                            #
##############################################


sub revcom {

	my ($dna)=@_;

# First we reverse the DNA

	my $revcom=reverse $dna;

#  Then we find the complement of each base

	$revcom=~tr/ACGTacgt/TGCAtgca/;

	return $revcom;
}

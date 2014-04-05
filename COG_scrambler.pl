#!/usr/bin/perl -w
# Federico Lauro
# The script will take two lists of COG numbers and associate to each COG it's COG category
# It will use 2 datasets, a REF(erence) and a QUE(ry). First it wil compare them by taking subsamples
# and calculating the median of REF-QUE of the subsample comparison between the datasets.
# It will then repeat the calculation on a scrabled dataset
# It doesn't really matter what dataset is REF and what is QUE: the calculation will always
# be REF-QUERY so that when it does the statistics the categories overrepresented in the REF will result in positive values
# and those under-represented will be negative.
# It will output the comparison medians and the confidence levels for each category
# It has two usage modes; 
# MODE 0 (default) will just work with a list of COGs and not consider the number
# of times that COG appears on a datset.
# MODE 1 also takes into account a multeplicity which is the number of times
# a certain COG appears and the input file must be <TAB> delimited with the multeplicity
# after each COG name.
# The COG assignment list has been generated from the whoa file downloaded from NCBI 
# issuing the command 'grep "COG" > destination file'


use strict;
use warnings;
use POSIX qw (ceil floor);
use vars qw($USAGE);
use Getopt::Long;

my @COG_FILE;	  # This is the array of COGs that is used to parse out the association COG-category
my $cog_line;	  # This is the string that gets associated with the above array during the database preparation
my $input_line;	  #
my $count_file;	  #
my $multiplier;	  #
my @COGS_1;	  # Array of the COGs of REFERENCE dataset
my @COG_cat1;	  # Array of the COG categories associated with the above array
my @COGS_2;	  # Array of the COGs of QUERY dataset
my @COG_cat2;	  # Array of the COG categories associated with the above array
my @COGS_3;	  # Array of the COGs of REFERENCE+QUERY dataset (this is the array used for picking the scrambled sets)
my @COG_cat3;	  # Array of the COG categories associated with the above array
my $i;            # used to initialize the for-next loops
my $j;            # used to initialize the for-next loops
my $length;	  # used to check on the columns of the COG_list inputs
my $length1;      # stores the length of the REFERENCE array of COGs
my $length2;      # stores the length of the QUERY array of COGs
my $rand1;	  # stores the random pointer for the REFERENCE array
my $rand2;        # stores the random pointer for the QUERY array
my $top_conf;     # sets the top score of the bootstrap to be considered: ceil(($bootstrap*$cutoff)-1)
my $bottom_conf;  # sets the bottom score of the bootstrap to be considered: floor ($bootstrap*(1-cutoff))
my @variables;

# This initializes all the variables for the REFERENCE
my $len_r;

   my $f_E_r; my $f_G_r; my $f_D_r; my $f_N_r; my $f_M_r;
   my $f_B_r; my $f_H_r; my $f_Z_r; my $f_V_r; my $f_C_r;
   my $f_W_r; my $f_S_r; my $f_R_r; my $f_P_r; my $f_U_r;
   my $f_I_r; my $f_Y_r; my $f_F_r; my $f_O_r; my $f_A_r;
   my $f_L_r; my $f_Q_r; my $f_T_r; my $f_K_r; my $f_J_r;


# This initializes all the variables for the query
my $len_q;

   my $f_E_q; my $f_G_q; my $f_D_q; my $f_N_q; my $f_M_q;
   my $f_B_q; my $f_H_q; my $f_Z_q; my $f_V_q; my $f_C_q;
   my $f_W_q; my $f_S_q; my $f_R_q; my $f_P_q; my $f_U_q;
   my $f_I_q; my $f_Y_q; my $f_F_q; my $f_O_q; my $f_A_q;
   my $f_L_q; my $f_Q_q; my $f_T_q; my $f_K_q; my $f_J_q;



# This initializes the array where the results of REFERENCE-QUERY are stored
# These arrays are re-used for the bootstrapping too

   my @f_E_m; my @f_G_m; my @f_D_m; my @f_N_m; my @f_M_m;
   my @f_B_m; my @f_H_m; my @f_Z_m; my @f_V_m; my @f_C_m;
   my @f_W_m; my @f_S_m; my @f_R_m; my @f_P_m; my @f_U_m;
   my @f_I_m; my @f_Y_m; my @f_F_m; my @f_O_m; my @f_A_m;
   my @f_L_m; my @f_Q_m; my @f_T_m; my @f_K_m; my @f_J_m;




$USAGE = "\nCOG_scrambler.pl [-h] -r COG_list -q COG_list -c COG_categories -o output -z sample_size\n\nOther options\n-f cutoff [0.97]\n-s random_seed [64238]\n-b bootstrap_rep [100]\n-m mode [0=no mutliples]\n";

my ($input_file1, $input_file2, $COG_input, $output_file, $help, $size, $verbose, $cutoff, $bootstraps, $seed, $mode) = (undef, undef, undef, undef, undef, undef, '1', '0.97', '100', '64238', '0');

&GetOptions('category|c=s'           => \$COG_input,
	    'output|o=s'             => \$output_file,
	    'reference|r=s'          => \$input_file1,
	    'query|q=s'              => \$input_file2,
	    'size|z=s'		     => \$size,
	    'seed|s=s'		     => \$seed,
	    'cutoff|f=s'	     => \$cutoff,
	    'bootstraps|b=s'	     => \$bootstraps,
	    'mode|m=s'               => \$mode,
	    'help|h'                 => \$help,
	    'verbose|v=s'            => \$verbose
	    );


if( $help ) {
    exec('perldoc', $0);
    die;
}

if( !defined $COG_input ) {
    die($USAGE . "\n\tMust specify a COG category linkage file\n\n");
}
if( !defined $input_file1 ) {
    die($USAGE . "\n\tMust specify an input file for the REFERENCE\n\n");
}
if( !defined $input_file2 ) {
    die($USAGE . "\n\tMust specify an input file for the QUERY\n\n");
}
if( !defined $output_file ) {
    die($USAGE . "\n\tMust specify an output file\n\n");
}
if( !defined $size ) {
    die($USAGE . "\n\tMust specify the subsample size\n\n");
}



# This just checks that the files can be opened
unless (open (INPUT, "$input_file1")) {print "Could not open the input file for the REFERENCE dataset","\n";
					exit;}
unless (open (COG, "$COG_input")) {print "Could not open the COG category file","\n";
					exit;}
@COG_FILE=<COG>;

print STDERR "\n\nPreparing the COG dataset for the REFERENCE...\n\n" if ( $verbose > 0 );

while ($input_line=<INPUT>)
{
  # discard if blank line
  if ($input_line =~/^\s*$/) {next;}
  chomp $input_line;
  @variables=split(/\s+/,$input_line);

  # This sets up the parsing for the one column input
  if ($mode==0) {
     $length=@variables;
     # This just checks that the array contains two elements (COG and multeplicity)
     if ($length==2) {print STDERR "WARNING: There is more than one column in the REFERENCE input, maybe you should be using MODE 1???\r" if ( $verbose > 0 );}
     $multiplier=1;
		}

   # This sets up the parsing for the two column (COG+multiplicity) input
   if ($mode==1) {
     $length=@variables;
     # This just checks that the array contains two elements (COG and multeplicity)
     unless ($length==2) {print "There's something fishy here, are you sure this is the correct TAB-delimited file???\n";
			exit;}
     $multiplier=$variables[1];
	          }

 foreach $cog_line (@COG_FILE) {
		# Parse the matching line to the COG list into the output
		$cog_line=~/^\[\w+\]\s\w+\s?/;
		my $mat=$&;
		$mat=~s/\[\w+\]//;
		$mat=~s/\s//g;
		if ($variables[0] eq $mat) {

			# Then print to the array the COG (@COGS) and category (@COG_cat1) for the matching line of genome 1
			$cog_line=~/^\[\w+\]?/;
			my $cog_cats=$&;
			$cog_cats=~s/\[//;
			$cog_cats=~s/\]//;
			# Here we take care of the multeplicity adding a specific COG in the array as many times as its multiplicity
			# This will make things easier later
			for ($i=0; $i<$multiplier; $i += 1) {
				push (@COGS_1,$mat);
				push (@COG_cat1,$cog_cats);}
			# This should speed up the script by ending the loop whenever a match is found
			last;

			}
		}
}

close INPUT;

# This just checks that the file for dataset 2 (QUERY) can be opened
unless (open (INPUT, "$input_file2")) {print "Could not open the input file for the QUERY dataset","\n";
					exit;}

print STDERR "\n\nPreparing the COG dataset for the QUERY...\n\n" if ( $verbose > 0 );

while ($input_line=<INPUT>)
{
  # discard if blank line
  if ($input_line =~/^\s*$/) {next;}
  chomp $input_line;
  @variables=split(/\s+/,$input_line);

  # This sets up the parsing for the one column input
  if ($mode==0) {
     $length=@variables;
     # This just checks that the array contains two elements (COG and multeplicity)
     if ($length==2) {print STDERR "WARNING: There is more than one column in the QUERY input, maybe you should be using MODE 1???\r" if ( $verbose > 0 );}
     $multiplier=1;
		}

   # This sets up the parsing for the two column (COG+multiplicity) input
   if ($mode==1) {
     $length=@variables;
     # This just checks that the array contains two elements (COG and multeplicity)
     unless ($length==2) {print "There's something fishy here, are you sure this is the correct TAB-delimited file???\n";
			exit;}
     $multiplier=$variables[1];
	          }

 foreach $cog_line (@COG_FILE) {
		# Parse the matching line to the COG list into the output
		$cog_line=~/^\[\w+\]\s\w+\s?/;
		my $mat=$&;
		$mat=~s/\[\w+\]//;
		$mat=~s/\s//g;
		if ($variables[0] eq $mat) {

			# Then print to the array the COG (@COGS) and category (@COG_cat1) for the matching line of genome 1
			$cog_line=~/^\[\w+\]?/;
			my $cog_cats=$&;
			$cog_cats=~s/\[//;
			$cog_cats=~s/\]//;
			# Here we take care of the multeplicity adding a specific COG in the array as many times as its multiplicity
			# This will make things easier later
			for ($i=0; $i<$multiplier; $i += 1) {
				push (@COGS_2,$mat);
				push (@COG_cat2,$cog_cats);}
			# This should speed up the script by ending the loop whenever a match is found
			last;

			}
		}
}

close INPUT;
close COG;

# Retrieveing the length of the 2 arrays of COGs
$length1=@COGS_1;
$length2=@COGS_2;

print STDERR "\n\nREFERENCE dataset contains ",$length1," COGs" if ( $verbose > 0 );
print STDERR "\n\nQUERY dataset contains ",$length2," COGs\n\n\n" if ( $verbose > 0 );

#######################################################
#                                                     #
#   FIRST I SUBSAMPLE AND ANALYZE THE ORIGINAL        #
#                                                     #
#    DATASET (N) TIMES AND COMPARE                    #
#                                                     #
#######################################################

print STDERR "\n\nResampling the original datasets...\n\n" if ( $verbose > 0 );

srand ($seed);
# The outer loop does the bootstrap repetitions
for ($i=0; $i<$bootstraps; $i += 1) {

 print STDERR "\rResampling cycle ",$i+1 if ( $verbose > 0 );;

	# I have to reset the frequencies each time
	# This initializes all the variables for the REFERENCE
	$len_r=0;

	   $f_E_r=0; $f_G_r=0; $f_D_r=0; $f_N_r=0; $f_M_r=0;
	   $f_B_r=0; $f_H_r=0; $f_Z_r=0; $f_V_r=0; $f_C_r=0;
	   $f_W_r=0; $f_S_r=0; $f_R_r=0; $f_P_r=0; $f_U_r=0;
	   $f_I_r=0; $f_Y_r=0; $f_F_r=0; $f_O_r=0; $f_A_r=0;
	   $f_L_r=0; $f_Q_r=0; $f_T_r=0; $f_K_r=0; $f_J_r=0;


	# This initializes all the variables for the query
	$len_q=0;

	   $f_E_q=0; $f_G_q=0; $f_D_q=0; $f_N_q=0; $f_M_q=0;
	   $f_B_q=0; $f_H_q=0; $f_Z_q=0; $f_V_q=0; $f_C_q=0;
	   $f_W_q=0; $f_S_q=0; $f_R_q=0; $f_P_q=0; $f_U_q=0;
	   $f_I_q=0; $f_Y_q=0; $f_F_q=0; $f_O_q=0; $f_A_q=0;
	   $f_L_q=0; $f_Q_q=0; $f_T_q=0; $f_K_q=0; $f_J_q=0;


		# The inner loop cycles through the subsample size
		for ($j=0; $j<$size; $j += 1) {
		$rand1=int(rand $length1);
		$rand2=int(rand $length2);
		COG_count ($COG_cat1[$rand1],$COG_cat2[$rand2]);
        	# At this stage, storing the COG names has no purpose
		# print $COGS_1[$rand1],"\t",$COGS_2[$rand2],"\n";
	}

  # Here I save the REFERENCE-QUERY % differences at the end of the result array
  push(@f_E_m,($f_E_r/$len_r-$f_E_q/$len_q)*100);
  push(@f_G_m,($f_G_r/$len_r-$f_G_q/$len_q)*100);
  push(@f_D_m,($f_D_r/$len_r-$f_D_q/$len_q)*100);
  push(@f_N_m,($f_N_r/$len_r-$f_N_q/$len_q)*100);
  push(@f_M_m,($f_M_r/$len_r-$f_M_q/$len_q)*100);
  push(@f_B_m,($f_B_r/$len_r-$f_B_q/$len_q)*100);
  push(@f_H_m,($f_H_r/$len_r-$f_H_q/$len_q)*100);
  push(@f_Z_m,($f_Z_r/$len_r-$f_Z_q/$len_q)*100);
  push(@f_V_m,($f_V_r/$len_r-$f_V_q/$len_q)*100);
  push(@f_C_m,($f_C_r/$len_r-$f_C_q/$len_q)*100);
  push(@f_W_m,($f_W_r/$len_r-$f_W_q/$len_q)*100);
  push(@f_S_m,($f_S_r/$len_r-$f_S_q/$len_q)*100);
  push(@f_R_m,($f_R_r/$len_r-$f_R_q/$len_q)*100);
  push(@f_P_m,($f_P_r/$len_r-$f_P_q/$len_q)*100);
  push(@f_U_m,($f_U_r/$len_r-$f_U_q/$len_q)*100);
  push(@f_I_m,($f_I_r/$len_r-$f_I_q/$len_q)*100);
  push(@f_Y_m,($f_Y_r/$len_r-$f_Y_q/$len_q)*100);
  push(@f_F_m,($f_F_r/$len_r-$f_F_q/$len_q)*100);
  push(@f_O_m,($f_O_r/$len_r-$f_O_q/$len_q)*100);
  push(@f_A_m,($f_A_r/$len_r-$f_A_q/$len_q)*100);
  push(@f_L_m,($f_L_r/$len_r-$f_L_q/$len_q)*100);
  push(@f_Q_m,($f_Q_r/$len_r-$f_Q_q/$len_q)*100);
  push(@f_T_m,($f_T_r/$len_r-$f_T_q/$len_q)*100);
  push(@f_K_m,($f_K_r/$len_r-$f_K_q/$len_q)*100);
  push(@f_J_m,($f_J_r/$len_r-$f_J_q/$len_q)*100);

}


# After the resampling of the original dataset, I calculate the median
# for each category of the differences in the original dataset and save them in
# strings
my $median_E = median(@f_E_m);
my $median_G = median(@f_G_m);
my $median_D = median(@f_D_m);
my $median_N = median(@f_N_m);
my $median_M = median(@f_M_m);
my $median_B = median(@f_B_m);
my $median_H = median(@f_H_m);
my $median_Z = median(@f_Z_m);
my $median_V = median(@f_V_m);
my $median_C = median(@f_C_m);
my $median_W = median(@f_W_m);
my $median_S = median(@f_S_m);
my $median_R = median(@f_R_m);
my $median_P = median(@f_P_m);
my $median_U = median(@f_U_m);
my $median_I = median(@f_I_m);
my $median_Y = median(@f_Y_m);
my $median_F = median(@f_F_m);
my $median_O = median(@f_O_m);
my $median_A = median(@f_A_m);
my $median_L = median(@f_L_m);
my $median_Q = median(@f_Q_m);
my $median_T = median(@f_T_m);
my $median_K = median(@f_K_m);
my $median_J = median(@f_J_m);


#######################################################
#                                                     #
#               THEN                                  #
#                                                     #
#       LET THE SCRAMBLING BEGIN......                #
#                                                     #
#                                                     #
#######################################################

# Re-initialize the results array for bootstrapping
   @f_E_m = (); @f_G_m = (); @f_D_m = (); @f_N_m = (); @f_M_m = ();
   @f_B_m = (); @f_H_m = (); @f_Z_m = (); @f_V_m = (); @f_C_m = ();
   @f_W_m = (); @f_S_m = (); @f_R_m = (); @f_P_m = (); @f_U_m = ();
   @f_I_m = (); @f_Y_m = (); @f_F_m = (); @f_O_m = (); @f_A_m = ();
   @f_L_m = (); @f_Q_m = (); @f_T_m = (); @f_K_m = (); @f_J_m = ();

# Merge the original sample arrays
@COGS_3 = (@COGS_1,@COGS_2);
@COG_cat3 = (@COG_cat1,@COG_cat2);
$length1=@COGS_3;

print STDERR "\n\nBootstrapping...\n\n" if ( $verbose > 0 );


# The outer loop does the bootstrap repetitions
for ($i=0; $i<$bootstraps; $i += 1) {

print STDERR "\rBootstrap cycle ",$i+1 if ( $verbose > 0 );

	# I have to reset the frequencies each time
	# This initializes all the variables for the REFERENCE
	$len_r=0;

	   $f_E_r=0; $f_G_r=0; $f_D_r=0; $f_N_r=0; $f_M_r=0;
	   $f_B_r=0; $f_H_r=0; $f_Z_r=0; $f_V_r=0; $f_C_r=0;
	   $f_W_r=0; $f_S_r=0; $f_R_r=0; $f_P_r=0; $f_U_r=0;
	   $f_I_r=0; $f_Y_r=0; $f_F_r=0; $f_O_r=0; $f_A_r=0;
	   $f_L_r=0; $f_Q_r=0; $f_T_r=0; $f_K_r=0; $f_J_r=0;


	# This initializes all the variables for the query
	$len_q=0;

	   $f_E_q=0; $f_G_q=0; $f_D_q=0; $f_N_q=0; $f_M_q=0;
	   $f_B_q=0; $f_H_q=0; $f_Z_q=0; $f_V_q=0; $f_C_q=0;
	   $f_W_q=0; $f_S_q=0; $f_R_q=0; $f_P_q=0; $f_U_q=0;
	   $f_I_q=0; $f_Y_q=0; $f_F_q=0; $f_O_q=0; $f_A_q=0;
	   $f_L_q=0; $f_Q_q=0; $f_T_q=0; $f_K_q=0; $f_J_q=0;


		# The inner loop cycles through the subsample size. This time I resample both datasets from
		# the same bin combined
		for ($j=0; $j<$size; $j += 1) {
		$rand1=int(rand $length1);
		$rand2=int(rand $length1);
		COG_count ($COG_cat3[$rand1],$COG_cat3[$rand2]);
        	# At this stage, storing the COG names has no purpose
#		print $COGS_1[$rand1],"\t",$COGS_2[$rand2],"\n";
	}

  # Here I save the REFERENCE-QUERY % differences of each category at the end of the result array
  push(@f_E_m,($f_E_r/$len_r-$f_E_q/$len_q)*100);
  push(@f_G_m,($f_G_r/$len_r-$f_G_q/$len_q)*100);
  push(@f_D_m,($f_D_r/$len_r-$f_D_q/$len_q)*100);
  push(@f_N_m,($f_N_r/$len_r-$f_N_q/$len_q)*100);
  push(@f_M_m,($f_M_r/$len_r-$f_M_q/$len_q)*100);
  push(@f_B_m,($f_B_r/$len_r-$f_B_q/$len_q)*100);
  push(@f_H_m,($f_H_r/$len_r-$f_H_q/$len_q)*100);
  push(@f_Z_m,($f_Z_r/$len_r-$f_Z_q/$len_q)*100);
  push(@f_V_m,($f_V_r/$len_r-$f_V_q/$len_q)*100);
  push(@f_C_m,($f_C_r/$len_r-$f_C_q/$len_q)*100);
  push(@f_W_m,($f_W_r/$len_r-$f_W_q/$len_q)*100);
  push(@f_S_m,($f_S_r/$len_r-$f_S_q/$len_q)*100);
  push(@f_R_m,($f_R_r/$len_r-$f_R_q/$len_q)*100);
  push(@f_P_m,($f_P_r/$len_r-$f_P_q/$len_q)*100);
  push(@f_U_m,($f_U_r/$len_r-$f_U_q/$len_q)*100);
  push(@f_I_m,($f_I_r/$len_r-$f_I_q/$len_q)*100);
  push(@f_Y_m,($f_Y_r/$len_r-$f_Y_q/$len_q)*100);
  push(@f_F_m,($f_F_r/$len_r-$f_F_q/$len_q)*100);
  push(@f_O_m,($f_O_r/$len_r-$f_O_q/$len_q)*100);
  push(@f_A_m,($f_A_r/$len_r-$f_A_q/$len_q)*100);
  push(@f_L_m,($f_L_r/$len_r-$f_L_q/$len_q)*100);
  push(@f_Q_m,($f_Q_r/$len_r-$f_Q_q/$len_q)*100);
  push(@f_T_m,($f_T_r/$len_r-$f_T_q/$len_q)*100);
  push(@f_K_m,($f_K_r/$len_r-$f_K_q/$len_q)*100);
  push(@f_J_m,($f_J_r/$len_r-$f_J_q/$len_q)*100);

}


print STDERR "\n\nCalculating the confidence intervals and writing to the output...\n" if ( $verbose > 0 );

# Here I sort the results of the boot strapping results in ascending order
@f_E_m = sort {$a <=> $b} @f_E_m;
@f_G_m = sort {$a <=> $b} @f_G_m;
@f_D_m = sort {$a <=> $b} @f_D_m;
@f_N_m = sort {$a <=> $b} @f_N_m;
@f_M_m = sort {$a <=> $b} @f_M_m;
@f_B_m = sort {$a <=> $b} @f_B_m;
@f_H_m = sort {$a <=> $b} @f_H_m;
@f_Z_m = sort {$a <=> $b} @f_Z_m;
@f_V_m = sort {$a <=> $b} @f_V_m;
@f_C_m = sort {$a <=> $b} @f_C_m;
@f_W_m = sort {$a <=> $b} @f_W_m;
@f_S_m = sort {$a <=> $b} @f_S_m;
@f_R_m = sort {$a <=> $b} @f_R_m;
@f_P_m = sort {$a <=> $b} @f_P_m;
@f_U_m = sort {$a <=> $b} @f_U_m;
@f_I_m = sort {$a <=> $b} @f_I_m;
@f_Y_m = sort {$a <=> $b} @f_Y_m;
@f_F_m = sort {$a <=> $b} @f_F_m;
@f_O_m = sort {$a <=> $b} @f_O_m;
@f_A_m = sort {$a <=> $b} @f_A_m;
@f_L_m = sort {$a <=> $b} @f_L_m;
@f_Q_m = sort {$a <=> $b} @f_Q_m;
@f_T_m = sort {$a <=> $b} @f_T_m;
@f_K_m = sort {$a <=> $b} @f_K_m;
@f_J_m = sort {$a <=> $b} @f_J_m;


# Select the ranges for the cutoffs
$bottom_conf = floor ($bootstraps*(1-$cutoff)-1);
$top_conf = ceil(($bootstraps*$cutoff));
# This just prevents an out-of-range pointer with small number of bootstraps
if ($bottom_conf < 0) {$bottom_conf=0}
if ($top_conf > ($bootstraps-1)) {$top_conf=$bootstraps-1}



unless (open (COUNTS, ">$output_file")) {print "Could not open the output file","\n";
					exit;}

print COUNTS "REFERENCE dataset: ",$input_file1, "\n";
print COUNTS "QUERY dataset: ",$input_file2, "\n";
print COUNTS "Subsample size: ",$size,"\n";
print COUNTS "Bootstrap replicates: ",$bootstraps,"\n";
print COUNTS "Confidence level: ",$cutoff,"\n";
print COUNTS "Category\tMedian\tBottom confidence\tTop confidence\n";
print COUNTS "E","\t",$median_E,"\t",$f_E_m[$bottom_conf],"\t",$f_E_m[$top_conf],"\n";
print COUNTS "G","\t",$median_G,"\t",$f_G_m[$bottom_conf],"\t",$f_G_m[$top_conf],"\n";
print COUNTS "D","\t",$median_D,"\t",$f_D_m[$bottom_conf],"\t",$f_D_m[$top_conf],"\n";
print COUNTS "N","\t",$median_N,"\t",$f_N_m[$bottom_conf],"\t",$f_N_m[$top_conf],"\n";
print COUNTS "M","\t",$median_M,"\t",$f_M_m[$bottom_conf],"\t",$f_M_m[$top_conf],"\n";
print COUNTS "B","\t",$median_B,"\t",$f_B_m[$bottom_conf],"\t",$f_B_m[$top_conf],"\n";
print COUNTS "H","\t",$median_H,"\t",$f_H_m[$bottom_conf],"\t",$f_H_m[$top_conf],"\n";
print COUNTS "Z","\t",$median_Z,"\t",$f_Z_m[$bottom_conf],"\t",$f_Z_m[$top_conf],"\n";
print COUNTS "V","\t",$median_V,"\t",$f_V_m[$bottom_conf],"\t",$f_V_m[$top_conf],"\n";
print COUNTS "C","\t",$median_C,"\t",$f_C_m[$bottom_conf],"\t",$f_C_m[$top_conf],"\n";
print COUNTS "W","\t",$median_W,"\t",$f_W_m[$bottom_conf],"\t",$f_W_m[$top_conf],"\n";
print COUNTS "S","\t",$median_S,"\t",$f_S_m[$bottom_conf],"\t",$f_S_m[$top_conf],"\n";
print COUNTS "R","\t",$median_R,"\t",$f_R_m[$bottom_conf],"\t",$f_R_m[$top_conf],"\n";
print COUNTS "P","\t",$median_P,"\t",$f_P_m[$bottom_conf],"\t",$f_P_m[$top_conf],"\n";
print COUNTS "U","\t",$median_U,"\t",$f_U_m[$bottom_conf],"\t",$f_U_m[$top_conf],"\n";
print COUNTS "I","\t",$median_I,"\t",$f_I_m[$bottom_conf],"\t",$f_I_m[$top_conf],"\n";
print COUNTS "Y","\t",$median_Y,"\t",$f_Y_m[$bottom_conf],"\t",$f_Y_m[$top_conf],"\n";
print COUNTS "F","\t",$median_F,"\t",$f_F_m[$bottom_conf],"\t",$f_F_m[$top_conf],"\n";
print COUNTS "O","\t",$median_O,"\t",$f_O_m[$bottom_conf],"\t",$f_O_m[$top_conf],"\n";
print COUNTS "A","\t",$median_A,"\t",$f_A_m[$bottom_conf],"\t",$f_A_m[$top_conf],"\n";
print COUNTS "L","\t",$median_L,"\t",$f_L_m[$bottom_conf],"\t",$f_L_m[$top_conf],"\n";
print COUNTS "Q","\t",$median_Q,"\t",$f_Q_m[$bottom_conf],"\t",$f_Q_m[$top_conf],"\n";
print COUNTS "T","\t",$median_T,"\t",$f_T_m[$bottom_conf],"\t",$f_T_m[$top_conf],"\n";
print COUNTS "K","\t",$median_K,"\t",$f_K_m[$bottom_conf],"\t",$f_K_m[$top_conf],"\n";
print COUNTS "J","\t",$median_J,"\t",$f_J_m[$bottom_conf],"\t",$f_J_m[$top_conf],"\n";


close COUNTS;


print STDERR "\n\nDONE\n" if ( $verbose > 0 );

exit;


####################################################################
#                                                                  #
#      Subroutine for counting the categories associated with      #
#                                                                  #
#                    a specific COG                                #
#                                                                  #
# NOTE: The multiplier (number of times a certain COG appears      #
# in a sample) has already been taken into account by making the   #
# arrays for the datasets                                          #
#                                                                  #
####################################################################

sub COG_count {

my (@sub_input) = @_;
my $cogs_1=$sub_input [0];
my $cogs_2=$sub_input [1];

# Increase the category counts for the REFERENCE
for (my $ii=0; $ii<(length($cogs_1)); $ii += 1)
	{$len_r=$len_r+1;
  		if(substr($cogs_1,$ii,1)=~/E/) {$f_E_r=$f_E_r+1}
	      	if(substr($cogs_1,$ii,1)=~/G/) {$f_G_r=$f_G_r+1}
  		if(substr($cogs_1,$ii,1)=~/D/) {$f_D_r=$f_D_r+1}
  		if(substr($cogs_1,$ii,1)=~/N/) {$f_N_r=$f_N_r+1}
  		if(substr($cogs_1,$ii,1)=~/M/) {$f_M_r=$f_M_r+1}
  		if(substr($cogs_1,$ii,1)=~/B/) {$f_B_r=$f_B_r+1}
  		if(substr($cogs_1,$ii,1)=~/H/) {$f_H_r=$f_H_r+1}
  		if(substr($cogs_1,$ii,1)=~/Z/) {$f_Z_r=$f_Z_r+1}
  		if(substr($cogs_1,$ii,1)=~/V/) {$f_V_r=$f_V_r+1}
  		if(substr($cogs_1,$ii,1)=~/C/) {$f_C_r=$f_C_r+1}
  		if(substr($cogs_1,$ii,1)=~/W/) {$f_W_r=$f_W_r+1}
  		if(substr($cogs_1,$ii,1)=~/S/) {$f_S_r=$f_S_r+1}
  		if(substr($cogs_1,$ii,1)=~/R/) {$f_R_r=$f_R_r+1}
  		if(substr($cogs_1,$ii,1)=~/P/) {$f_P_r=$f_P_r+1}
  		if(substr($cogs_1,$ii,1)=~/U/) {$f_U_r=$f_U_r+1}
  		if(substr($cogs_1,$ii,1)=~/I/) {$f_I_r=$f_I_r+1}
  		if(substr($cogs_1,$ii,1)=~/Y/) {$f_Y_r=$f_Y_r+1}
  		if(substr($cogs_1,$ii,1)=~/F/) {$f_F_r=$f_F_r+1}
  		if(substr($cogs_1,$ii,1)=~/O/) {$f_O_r=$f_O_r+1}
  		if(substr($cogs_1,$ii,1)=~/A/) {$f_A_r=$f_A_r+1}
		if(substr($cogs_1,$ii,1)=~/L/) {$f_L_r=$f_L_r+1}
		if(substr($cogs_1,$ii,1)=~/Q/) {$f_Q_r=$f_Q_r+1}
		if(substr($cogs_1,$ii,1)=~/T/) {$f_T_r=$f_T_r+1}
		if(substr($cogs_1,$ii,1)=~/K/) {$f_K_r=$f_K_r+1}
		if(substr($cogs_1,$ii,1)=~/J/) {$f_J_r=$f_J_r+1}

		}

# Increase the category counts for the QUERY
for (my $jj=0; $jj<(length($cogs_2)); $jj += 1)
	{$len_q=$len_q+1;
  		if(substr($cogs_2,$jj,1)=~/E/) {$f_E_q=$f_E_q+1}
	      	if(substr($cogs_2,$jj,1)=~/G/) {$f_G_q=$f_G_q+1}
  		if(substr($cogs_2,$jj,1)=~/D/) {$f_D_q=$f_D_q+1}
  		if(substr($cogs_2,$jj,1)=~/N/) {$f_N_q=$f_N_q+1}
  		if(substr($cogs_2,$jj,1)=~/M/) {$f_M_q=$f_M_q+1}
  		if(substr($cogs_2,$jj,1)=~/B/) {$f_B_q=$f_B_q+1}
  		if(substr($cogs_2,$jj,1)=~/H/) {$f_H_q=$f_H_q+1}
  		if(substr($cogs_2,$jj,1)=~/Z/) {$f_Z_q=$f_Z_q+1}
  		if(substr($cogs_2,$jj,1)=~/V/) {$f_V_q=$f_V_q+1}
  		if(substr($cogs_2,$jj,1)=~/C/) {$f_C_q=$f_C_q+1}
  		if(substr($cogs_2,$jj,1)=~/W/) {$f_W_q=$f_W_q+1}
  		if(substr($cogs_2,$jj,1)=~/S/) {$f_S_q=$f_S_q+1}
  		if(substr($cogs_2,$jj,1)=~/R/) {$f_R_q=$f_R_q+1}
  		if(substr($cogs_2,$jj,1)=~/P/) {$f_P_q=$f_P_q+1}
  		if(substr($cogs_2,$jj,1)=~/U/) {$f_U_q=$f_U_q+1}
  		if(substr($cogs_2,$jj,1)=~/I/) {$f_I_q=$f_I_q+1}
  		if(substr($cogs_2,$jj,1)=~/Y/) {$f_Y_q=$f_Y_q+1}
  		if(substr($cogs_2,$jj,1)=~/F/) {$f_F_q=$f_F_q+1}
  		if(substr($cogs_2,$jj,1)=~/O/) {$f_O_q=$f_O_q+1}
  		if(substr($cogs_2,$jj,1)=~/A/) {$f_A_q=$f_A_q+1}
		if(substr($cogs_2,$jj,1)=~/L/) {$f_L_q=$f_L_q+1}
		if(substr($cogs_2,$jj,1)=~/Q/) {$f_Q_q=$f_Q_q+1}
		if(substr($cogs_2,$jj,1)=~/T/) {$f_T_q=$f_T_q+1}
		if(substr($cogs_2,$jj,1)=~/K/) {$f_K_q=$f_K_q+1}
		if(substr($cogs_2,$jj,1)=~/J/) {$f_J_q=$f_J_q+1}

		}

 }


####################################################################
#                                                                  #
#      Subroutine for calculating the mean of an array after       #
#                                                                  #
#                      ordering it                                 #
#                                                                  #
####################################################################

sub median {

my $median;
my (@median_input) = @_;
my $count = @median_input+1;
@median_input = sort {$a <=> $b} @median_input;

if ($count % 2) {$median = $median_input[int($count/2)]}
	else
		{$median = ($median_input[$count/2] + $median_input[$count/2 - 1])/2}

return $median;
}




__END__



#
# COG_scrambler
#
# Federico Lauro

# POD documentation - main docs before the code

=head1 NAME

COG_scrambler.pl - script for comparing two lists of COGs to find the COG-categories over- and under-represnted in one dateset vs. the other

=head1 SYNOPSIS

% COG_scrambler.pl [-h] -r COG_list -q COG_list -c COG_categories -o output_prefix -z sample_size [-f -s -b -m]

=head1 DESCRIPTION

The script will take two lists of COG numbers and associate to each COG it's COG category
It will use 2 datasets, a REF(erence) and a QUE(ry). First it wil compare them by taking subsamples
and calculating the median of REF-QUE of the subsample comparison between the datasets.
It will then repeat the calculation on a scrabled dataset
It doesn't really matter what dataset is REF and what is QUE: the calculation will always
be REF-QUERY so that when it does the statistics the categories overrepresented in the REF will result in positive values
and those under-represented will be negative.
It will output the comparison medians and the confidence levels for each category

The COG_linkage file can been generated from the whoa file downloaded from NCBI issuing the command 'grep "COG" > destination file'

=head1 OPTIONS

-cutoff (-f) sets the cutoff for the confidence levels (default=0.97)

-size (-z) sets the size of each subsample to be compared (required)

-seed (-s) sets the random seed for resampling (default=64238)

-category (-c) name of the file with the COG category linkage file (required)

-reference (-r) name of the REFERENCE dataset (required)

-query (-q) name of the QUERY dataset (required)

-output (-o) sets the output filename (required)

-bootstraps (-b) number of bootstrap replicates to perform (default=100)

-mode (-m)

	MODE 0 (default) will just work with a list of COGs and not consider the number of times that COG appears on a datset.

	MODE 1 also takes into account a multeplicity which is the number of times a certain COG appears

-verbose (-v) use 0 for not displaying progress messages on screen (default=1)

-help (-h) displays this help page


=head1 AUTHOR - Federico Lauro

Email flauro@unsw.edu.au

=cut

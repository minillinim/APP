#!/usr/bin/perl
###############################################################################
#
#    app_normalise_by_tables.pl
#    
#    Find a centroid OTU table with a given number of individuals
#
#    Copyright (C) 2011 Michael Imelfort and Paul Dennis
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;

#CPAN modules
use Statistics::R;
use File::Basename;
use Data::Dumper;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright
printAtStart();

my $options = checkParams();

######################################################################
# CODE HERE
######################################################################
print "Checking if everything makes sense...\n";
# get a working directory
my $global_wd = "working";
if(exists $options->{'working'}) {$global_wd = $options->{'working'}; }
$global_wd .= "/";
`mkdir -p $global_wd`;

# output files
my $global_dist_fh = undef;
if(exists $options->{'dist'}) { open $global_dist_fh, ">", $options->{'dist'} or die "**ERROR: Could not open file: ".$options->{'dist'}." $!\n"; }

my $log_fn = $options->{'table'}.".log";
if(exists $options->{'log'}) { $log_fn = $options->{'log'}; } 
open my $global_log_fh, ">", $log_fn or die "**ERROR: Could not open file: $log_fn $!\n";

my $global_out_fn = $options->{'table'}.".normalised";
if(exists $options->{'out'}) { $global_out_fn = $options->{'out'}; } 

# work out how many samples we're dealing with
my $global_num_samples = find_number_of_samples($options->{'table'});

# work out the number of reps we need to do
my $global_norm_num_reps = 1000;
if(exists $options->{'reps'}) { $global_norm_num_reps = $options->{'reps'}; }

# make rarefactions
print "Making $global_norm_num_reps multiple rarefactions to $options->{'norm'} individuals...\n";
my $cmd_str = "multiple_rarefactions_even_depth.py -i ".$options->{'table'}." -o $global_wd -d ".$options->{'norm'}." -n $global_norm_num_reps --lineages_included --k";
`$cmd_str`;

# find centroid
print "Finding centroid table...\n";
my $centroid_index = find_centroid_table();
my $cp_str = "cp $global_wd"."/rarefaction_".$options->{'norm'}."_$centroid_index".".txt $global_out_fn";
`$cp_str`;

# clean up
print "Cleaning up...\n";
if(exists $options->{'dist'}) {close $global_dist_fh; }
close $global_log_fh;

######################################################################
# CUSTOM SUBS
######################################################################
sub find_number_of_samples
{
	#-----
	# Given an OTU table, determine the number of samples it has in it!
	#
	my ($table_fn) = @_;
	open my $t_fh, "<", $table_fn or die "**ERROR: Cannot open OTU table: $table_fn $! \n";
	my $line_count = 0;
	while(<$t_fh>)
	{
		# Take the fourth line in the file
		$line_count++;
		if($line_count < 4)
		{
			next;
		}
		else
	    {
	    	chomp $_;
	    	my @otu_fields = split /\t/, $_;
	    	return $#otu_fields - 1;
	    	last;
	    }
	}
	close $t_fh;
}

sub find_centroid_table
{
    #-----
    # Find a representative set of $global_norm_sample_size sequences (do this in RRRRRR...)
    #
    #### Create a communication bridge with R and start R
	my $R_instance = Statistics::R->new();
    $R_instance->start();

    my $sampl_p1 = $global_num_samples + 1;
    
    # read in the list of distance matricies
    $R_instance->run(qq`library(foreign);`);
    $R_instance->run(qq`a<-list.files("$global_wd", "*.txt");`);
    
    # work out how many there are and allocate an array
    $R_instance->run(qq`len_a <- length(a);`);
    $R_instance->run(qq`big_frame <- array(0,dim=c($global_num_samples,$global_num_samples,len_a));`);
    
    # load each file individually into a big frame
    my $r_str = "for (i in c(1:len_a)) { j <- i - 1; name <- paste(\"$global_wd\",\"rarefaction_".$options->{'norm'}."\",\"_\",j,\".txt\",sep=\"\"); u<-read.table(name,sep=\"\\t\",row.names=1); u[,$sampl_p1]<-NULL; big_frame[,,i]<-as.matrix(dist(t(u), upper=TRUE, diag=TRUE)); i<-i+1; }";
    $R_instance->run($r_str);    

    # find the average matrix
    $R_instance->run(qq`ave <- big_frame[,,1];`);
    $R_instance->run(qq`for (i in c(2:len_a)) { ave <- ave + big_frame[,,i]; }`);
    $R_instance->run(qq`ave <- ave/len_a;`);
    
    # find the euclidean distance of each matrix from the average
    $R_instance->run(qq`dist<-array(0,dim=c(len_a));`);
    $R_instance->run(qq`for (i in c(1:len_a)) { dist[i] <- sqrt(sum(big_frame[,,i]-ave)^2); }`);
    
    # find the min value
    $R_instance->run(qq`min_index <- which.min(dist);`);
    my $centroid_otu_index = $R_instance->get('min_index');
    
    # make stats on the distances
    # and log what we did
    $R_instance->run(qq`max_dist <- max(dist);`);
    $R_instance->run(qq`min_dist <- min(dist);`);
    $R_instance->run(qq`range_dist <- max_dist - min_dist;`);
    $R_instance->run(qq`mean_dist <- mean(dist);`);
    $R_instance->run(qq`median_dist <- median(dist);`);
    
    print $global_log_fh "---------------------------------------------------\n";
    print $global_log_fh "  Centroid OTU table based normalised statistics\n";
    print $global_log_fh "---------------------------------------------------\n";
    print $global_log_fh "Max dist:\t".$R_instance->get('max_dist')."\n";
    print $global_log_fh "Min dist:\t".$R_instance->get('min_dist')."\n";
    print $global_log_fh "Range:\t".$R_instance->get('range_dist')."\n";
    print $global_log_fh "Mean:\t".$R_instance->get('mean_dist')."\n";
    print $global_log_fh "Median:\t".$R_instance->get('median_dist')."\n";
  
    if(2 < $global_num_samples)
    {
        $R_instance->run(qq`library(permute);`);
        $R_instance->run(qq`library(vegan);`);
        $R_instance->run(qq`mantel.otu <- mantel(ave,big_frame[,,min_index]);`);
        $R_instance->run(qq`m_stat <- mantel.otu\$statistic;`);
        $R_instance->run(qq`m_sig <- mantel.otu\$signif;`);
        print $global_log_fh "Mantel P stat:\t".$R_instance->get('m_sig')."\n";
        print $global_log_fh "Mantel R stat:\t".$R_instance->get('m_stat')."\n";
    }
    else
    {
        print $global_log_fh "Too few samples to perform a mantel test.\n";
    }
    
    # print all the distances to a file so we can make purdy pictures from them later
    if(exists $options->{'dist'})
    {
	    my $num_tables = $R_instance->get('len_a');
	    foreach my $counter (1..$num_tables)
	    {
	        print $global_dist_fh $R_instance->get("dist[$counter]")."\n"         
	    }
    }   

    # let the user know the result
    return $centroid_otu_index - 1;     
}

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "out|o:s", "table|t:s", "norm|n:i", "log|l:s", "dist|d:s", "reps|r:i", "working|w:s");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsosy items
    #if(!exists $options{''} ) { print "**ERROR: \n"; exec("pod2usage $0"); }
    if(!exists $options{'table'} ) { print "**ERROR: You need to tell me which OTU table to normalise\n"; exec("pod2usage $0"); }
    if(!exists $options{'norm'} ) { print "**ERROR: You need to specify the number of sequences to normalise to\n"; exec("pod2usage $0"); }

    return \%options;
}

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2011 Michael Imelfort and Paul Dennis
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    app_normalise_by_tables.pl
    
=head1 COPYRIGHT

   copyright (C) 2011 Michael Imelfort and Paul Dennis

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

    Find a centroid OTU table with a given number of individuals

=head1 SYNOPSIS

    app_normalise_by_tables.pl -table|t OTU_TABLE -norm|n NORMALISATION_SIZE -log|l LOG_FILE
    
    Normalise a set of OTU tables

      -table -t OTU_TABLE          OTU table to normalise on
      -norm -n NORMALISATION_SIZE  Number of sequences to normalise to
      [-log -l LOG_FILE]           File to store results of mantel tests etc... (default: OTU_TABLE.log)
      [-out -o OUT_FILE]           Output notmalised file (default: OTU_TABLE.normalised)
      [-reps -r NUM_REPS]          Number of reps to take (default: 1000)
      [-working -w WORKING_DIR]    Place to put the multitude of files which will be created (default: working) 
      [-dist -d DIST_FILE]         File to store OTU table distances
      [-help -h]                   Displays basic usage information

=cut


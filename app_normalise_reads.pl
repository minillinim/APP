#!/usr/bin/perl
################################################################################
#
#    app_normalise_reads.pl
#    Version: 0.1
#
#    Use this script to normailise a pyrotag data file
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
################################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;

#CPAN modules
use List::Util 'shuffle';

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
open my $in_fh, "<", $options->{'in'} or die $!;
my $out_fn = $options->{'in'}.$options->{'num'}.".rand";
if(exists $options->{'out'}) { $out_fn = $options->{'out'}; }
open my $out_fh, ">", $out_fn or die $!;

my %global_sample_head_counts = ();
my %global_sample_head_line_nums = ();
my %global_sample_head_line_nums_rand = ();


# Get the total number of reads for each sample
# We assume the file looks like this:
#
# >BOB_dfjgdsf
# atcagatcagcatcgac
# dlsflkdjsflkjdslkjfds
# >SALLY_gsdfshsdk
# sjfdkjhfkdsfkds
# >BOB_reiuyewiu
# ...

my $line_num = 0;
while(<$in_fh>)
{
    $line_num++;
    if($_ =~ /^>/)
    {
        my @f1 = split /_/, $_;
        my @f2 = split />/, $f1[0];
        my $sample_name =$f2[1];
        if(exists($global_sample_head_counts{$sample_name}))
        {
            $global_sample_head_counts{$sample_name}++;
            
            push @{$global_sample_head_line_nums{$sample_name}}, $line_num;
        }
        else
        {
            # first time we've seen this sample
            my @tmp_array = ();
            push @tmp_array, $line_num;
            $global_sample_head_line_nums{$sample_name} = \@tmp_array;

            $global_sample_head_counts{$sample_name} = 1;
        }
    }
}
close $in_fh;

# check to see that there are at least num reads for each sample
# and choose which ones to keep and which to kill...
foreach my $key (keys %global_sample_head_counts)
{
    if($global_sample_head_counts{$key} < $options->{'num'})
    {
        print "Warning! $key has only $global_sample_head_counts{$key} reads.\n";
        print "This sample wil be ignored when normalising to: ".$options->{'num'}."\n";
        $global_sample_head_counts{$key} = 0;
    }
    else
    {
        # randomise the line nums
        my @tmp_array = shuffle @{$global_sample_head_line_nums{$key}};
        # get the first num entries and put them on the list
        foreach my $i (1..$options->{'num'})
        {
            $global_sample_head_line_nums_rand{$tmp_array[$i]} = 1;
        }
    }
}

#### Go through the file one more time and select the reads

open $in_fh, "<", $options->{'in'};
$line_num = 0;
my $in_fasta = -1;
my $seq = "";
my $header = "";
my $Y = 0;
my $N = 0;
while(<$in_fh>)
{
    $line_num++;
    if($_ =~ /^>/)
    {
        if(-1 != $in_fasta)
        {
            # we have previously been "in fasta"
            if(exists $global_sample_head_line_nums_rand{$in_fasta})
            {
                # this is a keeper
                print $out_fh "$header$seq\n";
            }
        }
        $seq = "";
        $header = $_;
        $in_fasta = $line_num;
    }
    else
    {
        chomp $_;
        $seq .= $_;
    }
}
$line_num++;
if(-1 != $in_fasta)
{
    if(exists $global_sample_head_line_nums_rand{$in_fasta})
    {
        # this is a keeper
        print $out_fh "$header$seq\n";
    }
}
close $in_fh;
close $out_fh;

######################################################################
# CUSTOM SUBS
######################################################################


######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "in|i:s", "num|n:s", "out|o:s" );
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
    if(!exists $options{'in'} ) { print "**ERROR: you MUST give a input file\n"; exec("pod2usage $0"); }
    if(!exists $options{'num'} ) { print "**ERROR: you MUST provide the number of reads to normalise to\n"; exec("pod2usage $0"); }
    #if(!exists $options{''} ) { print "**ERROR: \n"; exec("pod2usage $0"); }

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

    app_normalise_reads.pl

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

   Normalise the number of reads from each sample in a 454 pyrotag run

=head1 SYNOPSIS

    app_normalise_reads.pl  [-help|h]

      -in -i INFILE                Input file (must be shuffled fasta)
      -num -n NORMALISE_NUM        Number of reads to select
      [-out -o OUTFILE]            File to write to [default: INFILE_NORMALISE_NUM.rand]
      [-help -h]                   Displays basic usage information
      
=cut



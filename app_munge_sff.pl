#!/usr/bin/perl
###############################################################################
#
#    app_munge_sff.pl
#    
#    Splits a sff file into job specific folders and creates app_config files
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

#locally-written modules
#load the pretty names for the fields
use AppConfig;

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
#### GLOBALS
my $global_raw_sff = $options->{'prefix'}.".sff";
my $global_raw_fna = $options->{'prefix'}.".fna";
my $global_raw_qual = $options->{'prefix'}.".qual";
my $global_raw_sff_gz = $options->{'prefix'}.".sff.gz";
my $global_raw_fna_gz = $options->{'prefix'}.".fna.gz";
my $global_raw_qual_gz = $options->{'prefix'}.".qual.gz";
my $global_raw_mapping = $options->{'prefix'}.".pdbm";

my $global_can_remove_fna = 1;

# store all the unique job IDs
my %global_jobIDs = ();
my %global_jobFNAs = ();
my %global_jobQUALs = ();
my %global_jobCONFs = ();

# does the mapping file exist?
if(-e $global_raw_mapping)
{
    open my $map_fh, "<", $global_raw_mapping or die "Could not open $global_raw_mapping\n";
    while(<$map_fh>)
    {
        next if($_ =~ "^#");
        chomp $_;
        my @fields = split /\t/, $_;
        my @job_fields = split /\./, $fields[0];
        
        # push this string onto the hash so we can make a conf out of it
        if(! exists $global_jobIDs{$job_fields[0]})
        {
            my @tmp_SAMPs = ();
            $global_jobIDs{$job_fields[0]} = \@tmp_SAMPs;
        }
        push @{$global_jobIDs{$job_fields[0]}}, "$job_fields[1]\t$fields[$FNA{BarcodeSequence}]\t$fields[$FNA{LinkerPrimerSequence}]\t$fields[$FNA{Description}]";
    }
    close $map_fh; 
}
else
{
    print "**ERROR: mapping file must exist and be called: $global_raw_mapping\n";
    print "If you have a mapping file, please rename it and try again...\n";
    die;
}

my $global_split_args = "split_libraries.py -b variable_length -m $global_raw_mapping -f $global_raw_fna -a 2 -H 10 -M 1";

#### Check to see if sffinfo needs to be called
print "Checking whether we need to split up the sff file\n";
if(!(-e $global_raw_fna))
{
    print "Converting files using sffinfo...\t\t";
    # call sffinfo
    `sffinfo -q $global_raw_sff  > $global_raw_qual`;
    `sffinfo -s $global_raw_sff > $global_raw_fna`;
    print "Done\n";
}
else
{
    print "\n\t*** $global_raw_fna already exists - so using this.\n";
    print "\t*** If this is not what you want, delete or rename $global_raw_fna.\n\n";

    # we don;t want to remove these
    $global_can_remove_fna = 0;
}

if(-e $global_raw_qual)
{
    $global_split_args .= " -q $global_raw_qual ";
}

#### Split the fna and qual file
# first split it.
print "Splitting sequences by MID...\t\t";
`$global_split_args`;
print "Done\n";

#### Make the output directories and file handles and make the confs
print "Creating config files and job folders...\t\t";
foreach my $job_IDs (keys %global_jobIDs)
{
    # make the directory
    my $job_dir = "$APP_BYJOB/$job_IDs";
    `mkdir -p $job_dir`;
    
    # make the FNA and QUAL files
    open my $tmp_fna_fh, ">", "$job_dir/$job_IDs.fna";
    $global_jobFNAs{$job_IDs} = \$tmp_fna_fh;
    open my $tmp_qual_fh, ">", "$job_dir/$job_IDs.qual";
    $global_jobQUALs{$job_IDs} = \$tmp_qual_fh;
    
    # make the config files
    open my $tmp_conf_fh, ">", "$job_dir/app_".$job_IDs.".config";
    print $tmp_conf_fh "$FNA_HEADER\n";
    foreach my $sample_line (@{$global_jobIDs{$job_IDs}})
    {
        print $tmp_conf_fh $sample_line.$FNA_LINE_FINISHER;
    }
    print $tmp_conf_fh "$FNA_FOOTER\n";
    close $tmp_conf_fh;
}
print "Done\n";

#### Create the job specific fasta files
print "Splitting the fna file into multiple folders...\t\t";
split_fna_by_job();
print "Done\n";

#### Close all the files and clean up
foreach my $job_IDs (keys %global_jobIDs)
{
    close ${$global_jobFNAs{$job_IDs}};
    close ${$global_jobQUALs{$job_IDs}};
}

if(exists $options->{'cleanup'})
{
    # remove all the files we made
    if(1 == $global_can_remove_fna)
    {
        `rm $global_raw_fna`;
        `rm $global_raw_qual`;
    }
    `rm seqs.fna`;
    `rm histograms.txt`;
    `rm split_library_log.txt`;
}

######################################################################
# CUSTOM SUBS
######################################################################
sub split_fna_by_job
{
    #-----
    # split a full run's worth of reads into individual job folders
    #
    # first we make a lookup of readIDs to job IDs
    my %reads_2_jobs_hash = ();

    # open the raw file
    open my $fna_fh, "<" , "seqs.fna" or die $!;
    while(<$fna_fh>)
    {
        my $header = $_;
        
        # discard the sequence
        <$fna_fh>;
        
        # split on the 'dot'
        my @tmp_1 = split /\./, $header;
        # then split on ">"
        my @tmp_2 = split />/, $tmp_1[0];
        # then splt on spaces
        my @tmp_3 = split / /, $tmp_1[1];

        $reads_2_jobs_hash{$tmp_3[1]} = $tmp_2[1];
    }
    close $fna_fh;
    
    # now open the original fna and qual
    open $fna_fh, "<", $global_raw_fna or die $!;
    open my $qual_fh, "<", $global_raw_qual or die $!;
    my $seq = "";
    my $qual = "";
    my $header = "";
    my $qual_raw = "";
    my $in_fasta = 0;
    while(<$fna_fh>)
    {
        $qual_raw = <$qual_fh>;
        if($_ =~ /^>/)
        {
            if(0 != $in_fasta)
            {
                # not the first time we've been here
                # see if this guy is in the list and write away if so
                # we need to get the seqID
                # split on ">"
                my @tmp_1 = split />/, $header;
                # then splt on spaces
                my @tmp_2 = split / /, $tmp_1[1];
                if(exists $reads_2_jobs_hash{$tmp_2[0]})
                {
                    # this guy is in our lists!
                    my $fna_out_fh = ${$global_jobFNAs{$reads_2_jobs_hash{$tmp_2[0]}}};
                    my $qual_out_fh = ${$global_jobQUALs{$reads_2_jobs_hash{$tmp_2[0]}}};
                    print $fna_out_fh $header.$seq;
                    print $qual_out_fh $header.$qual;
                }
            }
            else
            {
                $in_fasta = 1;
            }
            # reset these guys
            $header = $_;
            $seq = "";
            $qual = "";
        }
        else
        {
            $seq .= $_;
            $qual .= $qual_raw;
        }
    }
    if(0 != $in_fasta)
    {
        # not the first time we've been here
        # see if this guy is in the list and write away if so
        # we need to get the seqID
        # split on ">"
        my @tmp_1 = split />/, $header;
        # then splt on spaces
        my @tmp_2 = split / /, $tmp_1[1];
        if(exists $reads_2_jobs_hash{$tmp_2[0]})
        {
            # this guy is in out lists!
            my $fna_out_fh = ${$global_jobFNAs{$reads_2_jobs_hash{$tmp_2[0]}}};
            my $qual_out_fh = ${$global_jobQUALs{$reads_2_jobs_hash{$tmp_2[0]}}};
            print $fna_out_fh $header.$seq;
            print $qual_out_fh $header.$qual;
        }
    }
    close $fna_fh;
    close $qual_fh;
}

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "prefix|p:s", "cleanup:s");
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
    if(!exists $options{'prefix'} ) { print "**ERROR: you MUST give a sff file\n"; exec("pod2usage $0"); }
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

    app_munge_sff.pl

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

   Splits a sff file into job specific folders and creates app_config files

=head1 SYNOPSIS

    app_munge_sff.pl -p|prefix SFF_PREFIX [-cleanup] [-help|h]

      -p SFF_PREFIX                Prefix of the sff, mapping, qual files etc..
      [-cleanup]                   Remove all temp files made.
      [-help -h]                   Displays basic usage information
         
=cut


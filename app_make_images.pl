#!/usr/bin/perl
###############################################################################
#
#    app_make_images.pl
#    
#    Make purdy looking images using results from APP
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
use File::Basename;
use Statistics::R;
use Data::Dumper;

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

#### Create a communication bridge with R and start R
my $global_R_instance = Statistics::R->new();
$global_R_instance->start();

getWorkingDirs($options->{'config'});

makeImageDirs(0);

# sample IDs
my %global_samp_ID_list = ();
my $global_num_samples = 0;

# analyses we need to do
my %global_otu_table_list = ();
my %global_heat_map_list = ();

# nice names for the OTUs
my %global_otu_nice_names_list = ();

# we need to modify these files
my $global_tn_modified = modify_otu($tn_otu_table_file, $global_TB_results_dir);
$global_otu_nice_names_list{$global_tn_modified} = "$global_TB_results_dir/$tn_prefix"."_hellinger_otu_table";
$global_otu_nice_names_list{$tn_otu_table_file} = "$global_TB_results_dir/$tn_prefix"."_otu_table";

my $global_sn_modified;

my %allowable_trans_types = ('weighted_unifrac' => 1,
                             'unweighted_unifrac' => 1,
                             );
#                             'euclidean' => 1,
#                             'hellinger' => 1,
#                             'binary_euclidean' => 1,
#                             'chord' => 1
#                             );

# there are a number of different ways to normalise
# by default don'r normalise
my $global_norm_style = "NONE";

print "Checking if all the config checks out...\t\t";
parse_config_images();
print "All good!\n";

# make a list of all the analyses we'll need to do
#### TABLE NORM IMAGES

$global_otu_table_list{ $global_tn_modified } = $global_TB_results_dir;

# always do these guys
my @otu_search_dirs = ("$global_TB_results_dir/beta_diversity/$nn_prefix/", "$global_TB_results_dir/beta_diversity/$tn_prefix/");
my @hm_search_dirs = ("$global_TB_results_dir/breakdown_by_taxonomy");

#### SEQ NORM IMAGES IF NEEDED
if($global_norm_style eq "SEQ")
{
    # make the seq processing dir if we need to
    makeImageDirs(1);
    
    # add OTU table files
    $global_sn_modified = modify_otu($sn_otu_table_file, $global_SB_results_dir);
    $global_otu_nice_names_list{$global_sn_modified} = "$global_SB_results_dir/$sn_prefix"."_hellinger_otu_table";
    $global_otu_nice_names_list{$sn_otu_table_file} = "$global_SB_results_dir/$sn_prefix"."_otu_table";
    $global_otu_table_list{ $global_sn_modified } = $global_SB_results_dir;
    push @otu_search_dirs, "$global_SB_results_dir/beta_diversity/";
    
    # now do heatmaps
    push @hm_search_dirs, "$global_SB_results_dir/breakdown_by_taxonomy";
}

# now process each of the search dirs.
foreach my $base_dir (@otu_search_dirs)
{
    my @matrix_types = `ls -1 $base_dir`;
    
    foreach my $mt (@matrix_types)
    {
        # check that this is a valid directory
        chomp $mt;
        next if(!exists($allowable_trans_types{$mt}));
        
        my $working_dir = $base_dir.$mt."/";
        my $otu_table = `ls -1 $working_dir | grep -e "otu_table.txt\$"`;
        chomp $otu_table;
        $global_otu_table_list{$working_dir.$otu_table} = $working_dir;
        my $otu_nice = substr($otu_table, 0, (length($otu_table) - 14));
        $global_otu_nice_names_list{$working_dir.$otu_table} = $working_dir.$otu_nice."_otu_table";
    }
}

# now make the OTU images -> no hellinger
print "Making OTU PCA (etc...) images\n";
foreach my $this_otu (keys %global_otu_table_list)
{
    make_otu_images($this_otu, $global_otu_table_list{$this_otu}, 0);
}

# now do the two hellinger types on the normalised otu tables
make_otu_images($global_tn_modified, $global_TB_results_dir, 1);
if($global_norm_style eq "SEQ")
{
    make_otu_images($global_sn_modified, $global_SB_results_dir, 1);
}
print "\n";

# process heatmap data
print "Making taxa heatmaps\n";
foreach my $base_dir (@hm_search_dirs)
{
    my $cmd = "ls -1 $base_dir | grep -e \"txt\$\"";
    my @tax_tables = `$cmd`;
    
    foreach my $tt (@tax_tables)
    {
        chomp $tt;
        my $tt_m = modify_otu_hm($base_dir."/".$tt);
        make_heatmap_images($tt_m);        
    }
}
print "\nDone!\n";
######################################################################
# CUSTOM SUBS
######################################################################

sub modify_otu
{
    #-----
    # copy and modify an otu table, remove execess crap
    # used when making pca plots
    #
    my ($otu_table, $working_dir) = @_;
    my $ret_table = $otu_table.".modified";
    if(! -e  $ret_table)
    {
        `sed -e "s/#rar.*//" -e "s/#OTU ID//" -e "/^\$/d" -e "s/\\t[^\\t]*\$//" $otu_table > $ret_table`;
    }
    return $ret_table; 
}

sub modify_otu_hm
{
    #-----
    # copy and modify an otu table, remove execess crap and change numbers to percentages
    # used whe making tax heatmaps
    #
    my ($otu_table) = @_;
    my $ret_table = $otu_table.".modified";
    if(1)#! -e  $ret_table)
    {
        # create the file
        open my $r_fh, ">", $ret_table or die $!;
        open my $o_fh, "<", $otu_table or die $!;
        while(<$o_fh>)
        {
            if($_ =~/^Taxon/) { $_ =~ s/Taxon//; print $r_fh $_; }
            else
            {
                chomp $_;
                my @line_fields = split /\t/, $_;
                print $r_fh $line_fields[0];
                foreach my $i (1..$#line_fields)
                {
                    # we need to go from 0.00122448979592 to 0.0012
                    my $print_num = sprintf("%.4f", $line_fields[$i]);
                    print $r_fh "\t$print_num";
                }
                print $r_fh "\n";
            }
        }
        close $r_fh;
        close $o_fh;
        
    }
    return $ret_table; 
}

sub make_heatmap_images
{
    #-----
    # run some R and make a few heatmaps
    #
    my ($otu_table) = @_;
    
    print ".";
    
    my $nice_dir = dirname($otu_table);
    my $full_name = basename($otu_table);
    my $neat_name = substr($full_name, 0, index($full_name, '.'));
    my $hm_svg = $neat_name.".svg";
    my $hm_pdf = $neat_name.".pdf";
    
    # shorten file names by setting the wd
    $global_R_instance->run(qq`start_dir <- getwd()`);
    $global_R_instance->run(qq`setwd("$nice_dir")`);
    
    $global_R_instance->run(qq`library (lattice)`);
    $global_R_instance->run(qq`time<-read.table("$otu_table",header=TRUE,row.names=1,sep="\\t")`);
    $global_R_instance->run(qq`t(time)->time`);
    
    $global_R_instance->run(qq`pdf(file='$hm_pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='letter', pointsize=12)`);
    $global_R_instance->run(qq`levelplot (time)`);
    $global_R_instance->run(qq`col.1<-colorRampPalette(c("white","grey","yellow","red"))(30)`);
    $global_R_instance->run(qq`levelplot (time, col.regions=col.1)`);
    $global_R_instance->run(qq`dev.off()`);

    $global_R_instance->run(qq`svg(filename='$hm_svg', height=6, width=6)`); 
    $global_R_instance->run(qq`levelplot (time)`);
    $global_R_instance->run(qq`col.1<-colorRampPalette(c("white","grey","yellow","red"))(30)`);
    $global_R_instance->run(qq`levelplot (time, col.regions=col.1)`);
    $global_R_instance->run(qq`dev.off()`);
        
    # put back the working directory before we leave
    $global_R_instance->run(qq`setwd(start_dir)`);
}

sub make_otu_images
{
    #-----
    # run some R magic and make a few images
    # 
    my($table, $working_dir, $hellinger) = @_;
    
    print ".";
    
    # before we start, get the current working directory
    $global_R_instance->run(qq`start_dir <- getwd()`);
    
    # make sure there a a blank line at the end of the file
    #`echo "\n" >> $table`;
    # remove any blank lines
    #`sed -i -e "/^\$/d" $table`;
    
    $global_R_instance->run(qq`library(vegan)`);
    $global_R_instance->run(qq`read.table("$table",header=TRUE,row.names=1)->tb.tmp`);
    if(1 == $hellinger)
    {
        $global_R_instance->run(qq`t(tb.tmp)->tb`);
    }
    else
    {   
        $global_R_instance->run(qq`t(sqrt(tb.tmp))->tb`);
    }
    
    # get the pretty name for this table
    my $nice_name = $global_otu_nice_names_list{$table};
    my $nice_dir = dirname($nice_name);
    
    # shorten file names by setting the wd
    $global_R_instance->run(qq`setwd("$nice_dir")`);
    
    # get the name and title
    $nice_name = basename($nice_name);
    my $nice_title = $nice_name;
    
    # fix this guy up
    $nice_name =~ s/unweighted_unifrac/uu/;
    $nice_name =~ s/weighted_unifrac/wu/;
    $nice_name =~ s/non_normalised/nn/;
    $nice_name =~ s/sequence_normalised/sn/;
    $nice_name =~ s/table_normalised/tn/;
   
    # PCA - get associated data and figures
    
    $global_R_instance->run(qq`rda(tb)->tb.pca`);
    
    my $species_file = $nice_name.".pca_species.scores.txt";
    my $sites_file = $nice_name.".pca_sites.scores.txt";
    $global_R_instance->run(qq`tmp_sc <- scores(tb.pca)`);
    $global_R_instance->run(qq`write.table(tmp_sc\$sites, "$sites_file")`);
    $global_R_instance->run(qq`write.table(tmp_sc\$species, "$species_file")`);
    
    my $pca_pdf = $nice_name.".pca.pdf";
    my $pca_svg = $nice_name.".pca.svg";
    
    my $pca_title = $nice_title.".pca";
    
    $global_R_instance->run(qq`pdf(file='$pca_pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='letter', pointsize=12)`);
    $global_R_instance->run(qq`plot(tb.pca,type='t',scaling=3,main='$pca_title')`);
    $global_R_instance->run(qq`dev.off()`);
    $global_R_instance->run(qq`svg(filename='$pca_svg', height=6, width=6)`); 
    $global_R_instance->run(qq`plot(tb.pca,main='$pca_title')`);
    $global_R_instance->run(qq`dev.off()`);
    
    # NMDS
    
    $global_R_instance->run(qq`metaMDS(tb,distance='euclidean')->tb.nmds`);

    my $nmds_file = $nice_name.".nmds.txt";
    my $nmds_scores_file = $nice_name.".nmds.scores.txt";
    #$global_R_instance->run(qq`write.table(tb.nmds, "$nmds_scores_file")`);
    $global_R_instance->run(qq`tmp_sc <- scores(tb.nmds)`);
    $global_R_instance->run(qq`write.table(tmp_sc, "$nmds_scores_file")`);

    my $nmds_pdf = $nice_name.".nmds.pdf";
    my $nmds_svg = $nice_name.".nmds.svg";
    
    my $nmds_title = $nice_title.".nmds";
    
    $global_R_instance->run(qq`pdf(file='$nmds_pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='letter', pointsize=12)`);
    $global_R_instance->run(qq`plot(tb.nmds,type='t',main='$nmds_title')`);
    $global_R_instance->run(qq`dev.off()`);
    $global_R_instance->run(qq`svg(filename='$nmds_svg', height=6, width=6)`);
    $global_R_instance->run(qq`plot(tb.nmds,type='t',main='$nmds_title')`);
    $global_R_instance->run(qq`dev.off()`);
        
    # Hierarchical clustering
    
    my $clust_pdf = $nice_name.".clust.pdf";
    my $clust_svg = $nice_name.".clust.svg";
    my $clust_title = $nice_title.": Hierarchical complete linkage clustering";

    $global_R_instance->run(qq`tb_clust <- hclust((vegdist(tb,method='euclidean')),method='complete')`);
    $global_R_instance->run(qq`pdf(file='$clust_pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='letter', pointsize=12)`);
    $global_R_instance->run(qq`plot(tb_clust,main='$nice_title')`);
    $global_R_instance->run(qq`dev.off()`);
    $global_R_instance->run(qq`svg(filename='$clust_svg', height=6, width=6)`); 
    $global_R_instance->run(qq`plot(tb_clust,main='$nice_title')`);
    $global_R_instance->run(qq`dev.off()`);
    
    # put back the working directory before we leave
    $global_R_instance->run(qq`setwd(start_dir)`);
    
    # done!
}

sub parse_config_images
{
    #-----
    # parse the app config file and produce a qiime mappings file
    #
    open my $conf_fh, "<", $options->{'config'} or die $!;
    while(<$conf_fh>)
    {
        next if($_ =~ /^#/);
        last if($_ =~ /^@/);
        chomp $_;
        my @fields = split /\t/, $_;
        
        # save the MID for later
        $global_samp_ID_list{$fields[$FNA{'SampleID'}]} = $fields[$FNA{'USE'}];
        
        # we need to find the globally maximal number of sequences for any USED sample
        if("1" eq int $fields[$FNA{'USE'}])
        {
            $global_num_samples++;
        }
    }
    
    # user options section
    while(<$conf_fh>)
    {
        chomp $_;
        my @fields = split /=/, $_;
        if($#fields > 0)
        {
            if($fields[0] eq "NORMALISE")
            {
                # is this guy set?
                if($fields[1] ne "")
                {
                    chomp $fields[1];
                    my @norm_fields = split /,/, $fields[1];
                    $global_norm_style = $norm_fields[0];

                    # check the normailisation style is sane
                    if($global_norm_style ne "SEQ")
                    {
                        die "You must specify 'SEQ' as normalisation methods (if you specify anything)\n";
                    }
                }
            }
        }
    }    
    close $conf_fh;
}

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "config|c:s");
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
    if(!exists $options{'config'} ) { print "**ERROR: you MUST give a config file\n"; exec("pod2usage $0"); }
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

    app_make_images.pl
    
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

    Make purdy looking images using results from APP

=head1 SYNOPSIS

    app_make_images.pl -c|config CONFIG_FILE [-help|h]
    
    Make purdy looking images using results from APP

      -c CONFIG_FILE               app config file to be processed
      [-help -h]                   Displays basic usage information

=cut


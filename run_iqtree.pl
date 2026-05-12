#!/usr/bin/perl -w

use strict;
use Parallel::ForkManager;

my $pm = 80;
my $set = $pm/4;
my $manager = new Parallel::ForkManager( $set );
my $dir = $ARGV[0]; #directory of your alignments
my $out = $ARGV[1]; #comma separated list of outgroups
chdir("$dir");
my @files=<*.*>;
chdir("../");

my $path = `pwd`;
chomp($path);

my $seed = int(rand(10000000));

my $time = `date`;
$time =~ s/\s+/_/g;
$time = substr($time, 0, -1);
`mkdir Trees_$time`;
chdir ("Trees_$time");
my $workingpath = $path . "/Trees_" . $time;

my $outgroup;
my @outgroups = split (/,/, $out);

my %outgroup_files;

OUTGROUP: for my $file (sort {$b cmp $a} @files){
        
		for my $tout (@outgroups){
			my @real_out = `grep "$tout" ../$dir/$file`;
			if(@real_out){
					my $true_out = substr($real_out[0], 1);
					$outgroup_files{$file}=$true_out;
					next OUTGROUP;
			}
		}
        
        
}


for my $file (keys %outgroup_files){

        $manager->start and next;
	`sbatch /mrm/bin/iqtree2.srun ../$dir/$file $file $outgroup_files{$file}`;
        #`/home/mmckain/bin/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -T 4 -f a -x 111123 -p 1233435 -# 500 -m GTRGAMMA -s ../$dir/$file -n $file.rax -o $outgroup_files{$file} > raxml_runt.out`;

        $manager->finish;
}

#!/usr/bin/perl -w
#use strict;
use Parallel::ForkManager;
my $proc=22;
my $pm = Parallel::ForkManager->new($proc);
my $infile    = $ARGV[0];
my $infiledup = $ARGV[1];

#my $nsp = $ARGV[2];
my $percent = $ARGV[2];
my $cutoff  = $ARGV[3];
my %gclust  = ();
my %groups  = ();
my ( $dup1, $dup2, %dups, %edups );
my ( $gid, $sid );
my %event;
&id_clusters();
&get_goodclusters();

my $time = `date`;
$time =~ s/\s+/_/g;
$time =~ s/\:/-/g;
$time =~ s/_$//g;
open my $log, ">", "Duplication_Pipeline_" . $time . ".log";  

`mkdir clusters`;
`mkdir alignments`;
`mkdir pal2nal\_alignments`;
`mkdir cleaned_alignments`;



my $j=1;
my $i=1;
for $gid (sort keys %gclust){
    if($j <= 1000){
	open my $setfile, ">>", "$i.set.txt";
	print $setfile "$gid:\t";
	for $sid (@{$gclust{$gid}}){
	    print $setfile "$sid\t";
	}
	print $setfile "\n";
	$j++;
    }
    else{
    	$j = 1;
    	$i++;
    	open my $setfile, ">>", "$i.set.txt";
	print $setfile "$gid:\t";
	for $sid (@{$gclust{$gid}}){
	    print $setfile "$sid\t";
	}
	print $setfile "\n";
	$j++;
    }
}

my @tasks = (1..$i);
TASKS:
for my $task (@tasks){

        $pm->start and next TASKS;
        `perl /grps2/mrmckain/bin/WGD_Identification_Pipeline/DuplicationAlignment_part2.pl $task.set.txt $percent $cutoff`;
	print $log "Running $task.set.txt...\n";

        $pm->finish;
}
$pm->wait_all_children;



###########
sub id_clusters {
    open my $DUP, "<", $infiledup;
    while (<$DUP>) {
        chomp;
        my ( $dup1, $dup2, $event ) = split /\s+/;
        if ( exists $dups{$dup1} ) {
            $dup1 .= "xx" . rand(100000000);
           # print "dup1 exists now it is $dup1\n";
        }
        #$event = substr($event, index($event, "-")+1);
        $dups{$dup1}        = $dup2;
        $edups{$dup1}{$dup2} = $event;
        
    }

    close($DUP);
}
#################
sub get_goodclusters {

open my $collapse_perblock, ">", "genefamiles_collapsed_by_event.txt";    
open my $orthogroups, "<", $ARGV[0];
while(<$orthogroups>){
	chomp;
	 my @current = split /\s+/;
        foreach my $cur (@current) {
            chomp;
            if ( $cur =~ /:$/ ) {
                $cur = substr( $cur, 0, -1 );
                $gid = $cur;
                next;
            }
            else{
                $groups{$gid}{$cur}=1;
            }
        }
    }
 my %goodgroups; 
 my %collapse;

 for my $orthid (sort keys %groups){
 #	print "$orthid\n";
	$collapse{$orthid}=();
 }
 
 my $ndup; 
=cut 
for my $dup ( sort keys %dups ) {
	my ($d1id, $d2id);
        if ( $dup =~ /xx/ ) {
            my $pos = index( $dup, "xx" );
            $ndup = substr( $dup, 0, $pos );
        }
        else { 
        	$ndup = $dup; 
		print "$ndup\n";
        }
        for my $key ( sort keys %groups ) {
		
            if (exists $groups{$key}{$ndup}) {
                    $d1id = $key;
		    $goodgroups{$d1id}=1;
		    print "key found\n";
                }
	    if (exists $groups{$key}{$dups{$ndup}}){
                    $d2id = $key;
   		    $goodgroups{$d2id}=1;
		    print "key found\n";
            }
        }
        
	
	if($d1id && $d2id){
	    if ($d1id eq $d2id){
        	next;
	    }
	    
	    else{
        	my $size1 = scalar keys %{$groups{$d1id}};
        	my $size2 = scalar keys %{$groups{$d2id}};
        	
        	if($size1 > $size2){
		    
        		for my $seqsid (sort keys %{$groups{$d2id}}){
			    $groups{$d1id}{$seqsid}=1;
        		}
			print $collapse_perblock "$d1id\t$d2id\t$ndup\t$dups{$ndup}\t$edups{$ndup}{$dups{$ndup}}\n";
        		delete $groups{$d2id};
			
        		for my $collid (sort keys %{$collapse{$d2id}}){
			    $collapse{$d1id}{$collid}=1;
        		}
			$collapse{$d1id}{$d2id}=1;

        		delete $collapse{$d2id};
		    }
        	
        	if($size2 > $size1){
		    
		    for my $seqsid (sort keys %{$groups{$d1id}}){
			$groups{$d2id}{$seqsid}=1;
		    }
		    print $collapse_perblock "$d2id\t$d1id\t$ndup\t$dups{$ndup}\t$edups{$ndup}{$dups{$ndup}}\n";
		    delete $groups{$d1id};
		    
		    for my $collid (sort keys %{$collapse{$d1id}}){
			$collapse{$d2id}{$collid}=1;
		    }
		    $collapse{$d2id}{$d1id}=1;
		    delete $collapse{$d1id};
        	}
	    }
	}
    }
=cut
for my $glid (sort keys %groups){
    #if(exists $goodgroups{$glid}){
	
	for my $gliseqs (sort keys %{$groups{$glid}}){
	    push( @{ $gclust{$glid} }, $gliseqs );
	}
   # }
}
	
open my $collapseoutfile, ">", $infiledup . ".condensed_clusters.txt";
for my $collgid (sort keys %collapse){
    print $collapseoutfile "$collgid:\t";
    for my $collseq (sort keys %{$collapse{$collgid}}){
	print $collapseoutfile "$collseq\t";
    }
    print $collapseoutfile "\n";
}
close $collapseoutfile;    		


}
#############



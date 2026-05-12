#!/usr/bin/perl -w
use strict;

my $percent = $ARGV[1];
my $cutoff  = $ARGV[2];
my %gclust  = ();
my ( $gid, $sid );
my $time = `date`;
$time =~ s/\s+/_/g;
$time =~ s/\:/-/g;
$time =~ s/_$//g;
open my $log, ">", "Duplication_Pipeline_Part2" . $time . ".log";

#my %skiptaxa;

#open my $skiptaxa, "<", "skip_22gtaxa.txt";
#while(<$skiptaxa>){
#	chomp;
#	$skiptaxa{$_}=1;
#}

open my $orthogroups, "<", $ARGV[0];
while(<$orthogroups>){
	chomp;
	 my @current = split /\s+/;
        for my $cur (@current) {
            chomp;
            if ( $cur =~ /:$/ ) {
                $cur = substr( $cur, 0, -1 );
                $gid = $cur;
                next;
            }
            else{
		my $tack=0;
		#for my $skip (sort keys %skiptaxa){
		#if($cur =~ /$skip/){
		#	$tack=1;
		#	}
		#}
		if($tack==0){
			push( @{$gclust{$gid}}, $cur );
		
		}
            }
        }
    }

&get_cdna_aa();
&run_mafft();
&run_pal2nal();
&alignment_cleaner();




##############
sub get_cdna_aa {
    
    my (%concatpep, %concatcdna);
    
    open my $cdnafileseq, "<", "concat.cdna";
    my $concatid;
    while(<$cdnafileseq>){
    	chomp;
	$_=~s/\s*$//g;
    	if(/>/){
		
    		$concatid=substr($_,1);
    	}
    	else{
    		$concatcdna{$concatid}.=$_;
   	}
    }
    close $cdnafileseq;
    
    open my $pepfileseq, "<", "concat.pep";
    while(<$pepfileseq>){
    	chomp;
	$_=~s/\s*$//g;
    	if(/>/){
    		$concatid=substr($_,1);
    	}
    	else{
    		$concatpep{$concatid}.=$_;
    	}
    }
    close $pepfileseq;
    
  
    for $gid ( sort keys %gclust ) {
    #	unless($gid == "9937"){
#		next;
#	}
    	print $log "Working on $gid...\n";
	open my $temp_pep, ">", "clusters/$gid.pep";
    	open my $temp_cdna, ">", "clusters/$gid.fsa";
       
        for $sid ( @{ $gclust{$gid} } ) {
	    print $log "Working on $sid...\n";
	    if(exists $concatpep{$sid} && $concatcdna{$sid}){
            print $temp_pep ">$sid\n$concatpep{$sid}\n";
            print $temp_cdna ">$sid\n$concatcdna{$sid}\n";
		}
	    else{
		if(!exists $concatpep{$sid}){
			print $log "$sid is not in the peptide file!\n";
		}
		if(!exists $concatcdna{$sid}){
                        print $log "$sid is not in the cDNA file!\n";
                }
              }
      }      
    }
   

}

##############
sub run_mafft {

   
   
    for $gid ( sort keys %gclust ) {
        print $log "Running mafft alignment on $gid...\n";
	   system "/grps2/mrmckain/bin/mafft-7.490-with-extensions/bin/mafft --quiet --auto clusters/$gid.pep > alignments/$gid.pep.align"; # pipe the muscle output 
    }
    
    
}

###########
sub run_pal2nal {
   
    
    for $gid ( sort keys %gclust ) {
	print $log "Running pal2nal for $gid...\n";
	    system "perl /mrm/bin/pal2nal.pl alignments/$gid.pep.align clusters/$gid.fsa -output fasta > pal2nal\_alignments/$gid.p2n.fasta\n";
    }
}

###########
sub alignment_cleaner {
    my $slen;
    my $cid;
    my $seq;
    my %alignlen;
    
    for $gid ( sort keys %gclust ) {
        #unless($gid == "1450"){
	#	next;
#	}
	print $log "Running alignment_cleaner for $gid...\n";
	my %gens     = ();
        my %newalign = ();
        my $pseqid;
        my $pseqs;
        open my $p2nfile, "<", "pal2nal_alignments/$gid.p2n.fasta";
        while (<$p2nfile>) {
        	chomp;
        	if(/>/){
        		if($pseqid){
        			$slen = length($pseqs);
            		my @current = split( //, $pseqs );
            		my $count = 0;
            		foreach my $sq (@current) {
                		$count++;
                		$gens{$pseqid}{$count} = $sq;
            		}
			$pseqs=();
			}
            		$pseqid=$_;
            	}
            
            else{
            	$pseqs.=$_;
            }	
        }
        $slen = length($pseqs);
        my @current = split( //, $pseqs );
        my $count = 0;
        foreach my $sq (@current) {
        	$count++;
            $gens{$pseqid}{$count} = $sq;
        }
       

        my $numk = scalar keys %gens;
        print "$numk\n";
        for ( my $i = 1 ; $i <= $slen ; $i++ ) {
            my $badid    = 0;
            my $keycount = 0;
            foreach my $keys ( sort keys(%gens) ) {
                $keycount++;
                if ( $gens{$keys}{$i} eq "-" ) {
                    $badid++;
		    
                }
                if (   ( $keycount == $numk )
                    && ( ( 1 - $badid / $numk ) >= $percent ) ){
                    foreach my $kys ( sort keys(%gens) ) {
                        if ( exists $newalign{$kys} ) {
                            $newalign{$kys} = $newalign{$kys} . $gens{$kys}{$i};
                        }
                        else { $newalign{$kys} = $gens{$kys}{$i}; }
                    }
                }
	
            }
        }

        foreach my $ids ( sort keys(%newalign) ) {
            my $stringlength = length( $newalign{$ids} );
#            print "$newalign{$ids}\n";
            my @cuent = split( //, $newalign{$ids} );
            my $oldseq = ();
            foreach my $sqq (@cuent) {
                if ( $sqq ne "-" ) {
                    $oldseq .= $sqq;
                }
            }
            my $oldlength = length($oldseq);
            if ( ( $oldlength / $stringlength ) < $cutoff ) {
                delete $newalign{$ids};
            }
        }

        open my $OUT1, ">",  "cleaned_alignments/$gid\_pairs.cleaned.fasta";
        foreach my $ky ( sort keys(%newalign) ) {
            my $kys;
            print $OUT1 "$ky\n$newalign{$ky}\n";
        }
    }
}


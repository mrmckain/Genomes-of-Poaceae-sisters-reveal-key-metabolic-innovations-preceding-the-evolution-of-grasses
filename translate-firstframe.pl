#!/usr/bin/perl -w
use strict;

my %seqs;
open my $cdnafile, "<", $ARGV[0];
my $sid;
while(<$cdnafile>){
        chomp;
        if(/>/){
                $sid = $_;
                next;
        }
        $seqs{$sid} .= $_;
}
close $cdnafile;
my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G', 'GCN'=>'A', 'CGN'=>'R', 'GGN'=>'G', 'CCN'=>'P', 'TCN'=>'S', 'ACN'=>'T', 'GTN'=>'V');

open my $OUT, ">",  $ARGV[1] . ".pep";
open my $OUT2, ">", $ARGV[1] .".fixed";
open my $OUT3, ">", $ARGV[1] . ".stops.txt";
open my $stopfile, ">", $ARGV[1] . ".seqswithstops.txt";
foreach my $seqid (sort keys %seqs){
        my $protein;    
        my $cdna;
        my $codon;
        for(my $i=0;$i<(length($seqs{$seqid})-2);$i+=3){
                $codon=substr($seqs{$seqid},$i,3);
                $codon= uc $codon;
                if (exists $g{$codon}){
                        if($g{$codon} eq "\_"){
                                                              print $stopfile "$seqid\n";
							      print $OUT3 "$seqid\n";
                        }
                        else{   
                                $protein .= $g{$codon};
                        }
                }
                else{
                        print "Bad codon: $codon\n";
                }
                unless(exists $g{$codon}){
                        #if ($codon =~ /N|X/){
                                next; 
				#$protein .= "X";
                        #}
                }
                if (exists $g{$codon} && $g{$codon} eq "\_"){
                        print "STOP is codon: $codon\n";
			print $OUT3 "$codon\n";
			
                }
                else{
                        $cdna .= $codon;
                }
        }
        print $OUT "$seqid\n$protein\n";
        print $OUT2 "$seqid\n$cdna\n";
}
close $OUT;
close $OUT2;
close $stopfile;


#!/usr/bin/perl -w

#usage： perl flt_bam.pl  /data/frr/Qu_methylation/all_bam /data/frr/Qu_methylation/flt_bam

use strict;
use File::Basename;

my $bam_dir = shift;   #/data/frr/Qu_methylation/all_bam
opendir (DIR, $bam_dir) or die "can't open the directory!";
my @files = grep/\.bam$/,readdir DIR;
closedir(DIR);
my $path = shift;      #/data/frr/Qu_methylation/flt_bam

foreach my $file (@files){
    open(IN, "samtools view -h $file |") or die "Cannot open $file:$!";
	my $basename=basename($file, '.bam'); 
	my $outfile=$path.'/'.$basename.'_flt.bam';
	open(OUT, "| samtools view -b -@ 50 - > $outfile") or die "Cannot write into $outfile:$!";
		while(<IN>){
			if(/^\@/){
				print OUT $_;
			}else{
				my ($flag, $mate_chr)=(split/\s+/)[1,6];  #my ($flag, $mq, $mate_chr)=(split/\s+/)[1,4,6];
				if($mate_chr eq '=' and ($flag & 2 and ! ($flag & 4) and ! ($flag & 8))){
					print OUT $_;
				}
			}
		}
	}	
		close(IN);
		close(OUT);	

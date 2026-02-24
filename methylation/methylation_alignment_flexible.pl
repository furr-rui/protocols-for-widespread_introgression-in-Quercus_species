#!/usr/bin/perl -w
# usage: 
# a.mkdir v #make a directory to store cleandata
# b. perl methylation_alignment_flexible.pl my_reference.fa cleandata

use strict;
our $ref_file = shift;
#my $out_path=shift;
our $out_path='/data/frr/Qu_methylation/all_bam';
our $dir = shift;
opendir (DIR, $dir) or die "can't open the directory!";
my @files = grep/\.fq.gz$/,readdir DIR;
closedir(DIR);
#print join(",",@files),"\n";
#print "$files[0]\n";
#my ($tag,$order) = $files[0] =~/^(.+)_([1|2]).clean.fq.gz$/;
#print "$tag\t$order\n";
my %lib=();
for(my $i=0;$i<=$#files;$i++){
     my ($tag,$order) = $files[$i] =~/^(.+)_([1|2]).clean.fq.gz$/;      #小括号里建议再详细写
     $lib{$tag}=1;                                                          
    #print "$tag\t$lib{$tag}\n";
 }
    #print "$tag\t$order\n";   
foreach (keys %lib){
	my $fq1 = $_.'_1.clean.fq.gz';                                 #.表示字符串拼接
	my $fq2 = $_.'_2.clean.fq.gz';
	#print "$_\t$lib{$_}\n";
	if(-e $fq1 and -e $fq2){                                              #-e判断文件或目录名存不存在
       # print "$fq1\t$fq2\n";
	system("biscuit align -@ 50 -R \"\@RG\\tID:$_\\tSM:$_\" $ref_file ./$fq1 ./$fq2 | dupsifter $ref_file | samtools sort -m 4g -@ 30 -o $out_path/$_\'_srt.markdup.bam\' -O BAM -");
	}else{
	print STDERR "found error";
    }
}

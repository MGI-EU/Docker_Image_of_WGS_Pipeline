#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Basename;

die "Usage:perl $0 <in.picard.bam>" unless @ARGV==1;
my $bam = shift;
my $clean_reads=0;
my $clean_bases=0;
my $mapped_reads=0;
my $mapped_bases=0;
my $mapping_rate=0;
my $dup_reads=0;
my $dup_rate=0;
my $mismatch_bases=0;
my $mismatch_rate=0;
my $uniq_reads=0;
my $uniq_bases=0;
my $uniq_rate=0;

open BAM,"samtools view  $bam | " or die;
while(<BAM>)
{
		chomp;
		my $info=$_;
		my @info=split /\s+/,$info;
		$clean_reads++;
		$clean_bases+=length($info[9]);
		unless($info[1]=~/u/)
		{
			if($info[4] >= 1)
			{
				if($info=~/XC:i:(\d+)/){$uniq_reads++;$uniq_bases+=$1;}
				else{$uniq_reads++;$uniq_bases+=length($info[9]);}
			}
		  	if($info=~/XC:i:(\d+)/){$mapped_reads++;$mapped_bases+=$1;}
			else{$mapped_reads++;$mapped_bases+=length($info[9]);}
			if($info[1]=~/d/){$dup_reads++;}
			if($info=~/XM:i:(\d+)/){$mismatch_bases+=$1;}
		}
}
$mismatch_rate=$mismatch_bases/$mapped_bases;
$mapping_rate=$mapped_reads/$clean_reads;
$dup_rate=$dup_reads/$mapped_reads;
$uniq_rate=$uniq_reads/$mapped_reads;

my $name=basename($bam);
my $sample=(split /\./,$name)[0];
print "Sample\t$sample\n";
print "Clean reads\t$clean_reads\n";
print "Clean bases(bp)\t$clean_bases\n";
print "Mapped reads\t$mapped_reads\n";
print "Mapped bases(bp)\t$mapped_bases\n";
print "Mapping rate\t",sprintf("%.2f%%",100*$mapping_rate),"\n";
print "Uniq reads\t$uniq_reads\n";
print "Uniq bases(bp)\t$uniq_bases\n";
print "Unique rate\t",sprintf("%.2f%%",100*$uniq_rate),"\n";
print "Duplicate reads\t$dup_reads\n";
print "Duplicate rate\t",sprintf("%.2f%%",100*$dup_rate),"\n";
print "Mismatch bases(bp)\t$mismatch_bases\n";
print "Mismatch rate\t",sprintf("%.2f%%",100*$mismatch_rate),"\n";

close $bam;

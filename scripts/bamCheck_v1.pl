#!/usr/bin perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin); ## bin path for picard and samtools

my $usage=<<'USAGE';
    #modidy from /hwfssz1/ST_CANCER/CGR/SHARE/CancerPipeline/bin/bamCheck_v1.pl	
    #BGISEQ-500 data could not run picard ValidateSamFile (PL=COMPLETE)
    This program is designed to check the bam files for cancer pipline
                        version 1.0 : 2015-11-16
                             author : luoshuzhen@genomics.cn  
    Options:
         -a STR  fq.stat/bam file for check  [require]
         -b STR  bam file for save [require]
         -r STR  remove the file from previous step,0 for remove ,1 not remove [defaut 0]
         -h Help Information
USAGE

my ($checkFile,$saveFile,$rm);

GetOptions(
    'a=s' => \$checkFile,
    'b=s' => \$saveFile,
    'r=s' => \$rm,
    'h|?' => sub{die "$usage\n";},
);

die "$usage\n" unless (defined($checkFile) && defined($saveFile));

$rm ||=0;

my($cleanReads,$bamReads);
my($F1,$F2,$var1,$var2);
my $flag="bam";

if($checkFile=~/\.bam$/){ 
    $checkFile =~ s/\.bam$//;
    open IA, "$checkFile.stat" || die $!; ## .stat file is from samtools flagstat
}else{
    $flag="clean";
    open IA, $checkFile || die $!;
}
$saveFile =~ s/\.bam$//;
open IB, "$saveFile.stat" || die $!;

while( defined($F1=<IA>) and defined($F2=<IB>)){
    chomp($F1);
    chomp($F2);
    if($F1 =~ /^Number of Reads/i){ $cleanReads=(split /\s+/,$F1)[-1];
    }elsif($F1 =~ /in total/i){ ($cleanReads)= $F1 =~ /^(\d+)/;}
    if($F2 =~ /in total/i){ ($bamReads)= $F2 =~ /^(\d+)/;
    }
}
close IA;
close IB;

my $n1=(split /\//,$checkFile)[-1]; ## get the file name
my $n2=(split /\//,$saveFile)[-1];
print "Type\t$n1\t$n2\n";
print "cleanReads\t$cleanReads\t$bamReads\n";

if($bamReads eq $cleanReads ){
    if($flag eq "bam"){
        if(-e "$checkFile.Val.txt"){ $var1=`grep "No errors found" $checkFile.Val.txt`;}
        else{$var1="No errors found";} # BGISEQ-500 data could not run picard ValidateSamFile
	if(-e "$saveFile.Val.txt"){ $var2=`grep "No errors found" $saveFile.Val.txt`;}
	else{$var2="No errors found";} # BGISEQ-500 data could not run picard ValidateSamFile
        chomp($var1);
        chomp($var2);
        if($var1 && $var2){
            chomp($var1);
            chomp($var2);
            print "ValidateFile\t$var1\t$var2\n";
            if($rm==0){
                system"rm -rf $checkFile.bam";
	            print "$checkFile.bam had removed !\n";
            }
        }else{
	        print "bam file maybe incompleted or does not check ,you can remove it after checking\n";
        }
    }else{
        if(-e "$saveFile.Val.txt"){ $var2=`grep "No errors found" $saveFile.Val.txt`;}
	else{$var2="No errors found";} # BGISEQ-500 data could not run picard ValidateSamFile
        if($var2){print "$saveFile.bam has no errors\n";}
    }
}else{
    print "reads maybe inconsistent,you must check it\n";
}



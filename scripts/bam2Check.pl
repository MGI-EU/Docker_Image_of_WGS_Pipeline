#!/usr/bin perl -w
use strict;

print STDERR "BGISEQ-500 data could not run picard ValidateSamFile (PL=COMPLETE)";
die "Usage:\tperl $0 <bam.list> <check.bam> \n"if(@ARGV!=2);
my $bamList=$ARGV[0]; ## bam list for remove
my $checkfile=$ARGV[1];

my $flag=0;
my $rebam='';

## calculate total reads for merge.bam
$checkfile =~ s/\.bam$//;
open CH,"$checkfile.stat" || die $!;
my ($checkReads,$checkvar);
print "Type\tcleanReads\tValidateFile\n";
while(<CH>){
    chomp;
    if($_ =~ /in total/i){ ($checkReads)= $_ =~ /^(\d+)/; }
}
close CH;
my $name=(split /\//,$checkfile)[-1];
if(-e "$checkfile.Val.txt"){ $checkvar=`grep "No errors found" $checkfile.Val.txt`;}
else{$checkvar="No errors found"} #BGISEQ-500 data could not run picard ValidateSamFile
if($checkvar){
    chomp($checkvar);
    print "$name\t$checkReads\t$checkvar\n";
}else{ 
    print "$name\t$checkReads\n";
}

## calculate reads for per bam and count the sum
open IN,$bamList || die $!;
my $totalReads=0;
my $num=0;
while(my $bam=<IN>){
    chomp($bam);
    $num +=1;
    if(defined $bam){
        $bam =~ s/\.bam$//;
        my $name=(split /\//,$bam)[-1];
        open IA, "$bam.stat" || die $!;
        my $bamRead;
        while(<IA>){
            chomp;
            if($_ =~ /in total/i){ ($bamRead)= $_ =~ /^(\d+)/; }
        }
        close IA;
        $totalReads += $bamRead;
        my $var;
        if(-e "$bam.Val.txt"){ $var=`grep "No errors found" $bam.Val.txt`;}
	else{$var="No errors found"} #BGISEQ-500 data could not run picard ValidateSamFile
        if($var) { $flag += 1; }
        print "$name\t$bamRead\t$var";
        $rebam .= "$bam.bam ";
    }
}
close IN;
chop($rebam);

if( $totalReads == $checkReads && defined $checkvar){
    if($num==$flag){
        print "totalReads\t$totalReads\tNO errors\n";
        system("rm -rf $rebam");
        print "$rebam has removed";
    }else{ print "bam file maybe incompleted or does not check ,you can remove it after checking\n";}
}else{
    print "reads maybe inconsistent,you must check it\n";
}



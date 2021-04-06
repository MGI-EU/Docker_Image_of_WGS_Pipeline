#!/usr/bin/perl -w

use strict;

die "Usage:\tperl $0 <in.list> <in.config> <outdir> \n"if(@ARGV!=3);

my($List,$config,$outdir)=@ARGV;
my(%samples,%LIBS,%SortBAM);
my(%hash_sampair); 

my($BIN,$DATAPATH,$LIBPATH,$reference,$FileName,$platform,$center,$java8); ## basic rm $queue
my($SOAPnuke,$SOAPnukePara,$fqcheck,$bamcheck,$sortcheck,$bwacheck,$mkdupcheck,$coveragecheck,$bwa,$bwaMemPara,$samtools);# fq clean and alignment 
my($picard,$mkdupPara); ## mark duplication
my ($run_deepvariant,$container,$deepvariant_SIF,$deepvariant,$deepvariant_version,$deepvariant_mode,$deepvariant_bed,$deepvariant_core,$deepvariant_mem,$deepvariant_bam,$GNUdir,$gpu); ## SNP calling
my ($CLEAN_job_name,$CLEAN_nodes,$CLEAN_cpus_per_task,$CLEAN_mem,$CLEAN_time);
my ($ALIGN_job_name,$ALIGN_nodes,$ALIGN_cpus_per_task,$ALIGN_mem,$ALIGN_time);
my ($CKBWA_job_name,$CKBWA_nodes,$CKBWA_cpus_per_task,$CKBWA_mem,$CKBWA_time);
my ($SORT_job_name,$SORT_nodes,$SORT_cpus_per_task,$SORT_mem,$SORT_time);
my ($CKBAM_job_name,$CKBAM_nodes,$CKBAM_cpus_per_task,$CKBAM_mem,$CKBAM_time);
my ($MKDUP_job_name,$MKDUP_nodes,$MKDUP_cpus_per_task,$MKDUP_mem,$MKDUP_time);
my ($CKDUP_job_name,$CKDUP_nodes,$CKDUP_cpus_per_task,$CKDUP_mem,$CKDUP_time);
my ($COVER_job_name,$COVER_nodes,$COVER_cpus_per_task,$COVER_mem,$COVER_time);
my ($DEEP_job_name,$DEEP_nodes,$DEEP_cpus_per_task,$DEEP_mem,$DEEP_time);

################################ Read configuration file ###################################

open CFG,"$config" or die $!;
while(<CFG>){
    chomp;
    if($_=~/^\#/){next;}
    ## basic parameters
    if($_=~/^FileName=(.*?);/){$FileName=$1;}
    if($_=~/^reference=(.*?);/){$reference=$1;}
    if($_=~/^DATAPATH=(.*?);/){$DATAPATH=$1;}
    if($_=~/^BIN=(.*?);/){$BIN=$1;}
    if($_=~/^LIBPATH=(.*?);/){$LIBPATH=$1;}
    if($_=~/^SeqPlatform=(.*?);/){$platform=$1;}
    if($_=~/^SeqCenter=(.*?);/){$center=$1;}
    ##fastq clean
    if($_=~/^SOAPnuke=(.*?);/){$SOAPnuke=$1;}
    if($_=~/^SOAPnukePara=(.*?);/){$SOAPnukePara=$1;}
    if($_=~/^fqcheck=(.*?);/){$fqcheck=$1;}
    if($_=~/^CLEAN_job_name=(.*?);/){$CLEAN_job_name=$1;}
    if($_=~/^CLEAN_nodes=(.*?);/){$CLEAN_nodes=$1;}
    if($_=~/^CLEAN_cpus_per_task=(.*?);/){$CLEAN_cpus_per_task=$1;}
    if($_=~/^CLEAN_mem=(.*?);/){$CLEAN_mem=$1;}
    if($_=~/^CLEAN_time=(.*?);/){$CLEAN_time=$1;}

    ##bwa mem
    if($_=~/^bwa=(.*?);/){$bwa=$1;}
    if($_=~/^bwaMemPara=(.*?);/){$bwaMemPara=$1;}
    if($_=~/^samtools=(.*?);/){$samtools=$1;}
    if($_=~/^BAMcheck=(.*?);/){$bamcheck=$1;}
    if($_=~/^ALIGN_job_name=(.*?);/){$ALIGN_job_name=$1;}
    if($_=~/^ALIGN_nodes=(.*?);/){$ALIGN_nodes=$1;}
    if($_=~/^ALIGN_cpus_per_task=(.*?);/){$ALIGN_cpus_per_task=$1;}
    if($_=~/^ALIGN_mem=(.*?);/){$ALIGN_mem=$1;}
    if($_=~/^ALIGN_time=(.*?);/){$ALIGN_time=$1;}
    if($_=~/^BwaCheck=(.*?);/){$bwacheck=$1;} 

    if($_=~/^CKBWA_job_name=(.*?);/){$CKBWA_job_name=$1;}
    if($_=~/^CKBWA_nodes=(.*?);/){$CKBWA_nodes=$1;}
    if($_=~/^CKBWA_cpus_per_task=(.*?);/){$CKBWA_cpus_per_task=$1;}
    if($_=~/^CKBWA_mem=(.*?);/){$CKBWA_mem=$1;}
    if($_=~/^CKBWA_time=(.*?);/){$CKBWA_time=$1;}
	
    ##sort and merge bams
    if($_=~/^SortCheck=(.*?);/){$sortcheck=$1;}
    if($_=~/^SORT_job_name=(.*?);/){$SORT_job_name=$1;}
    if($_=~/^SORT_nodes=(.*?);/){$SORT_nodes=$1;}
    if($_=~/^SORT_cpus_per_task=(.*?);/){$SORT_cpus_per_task=$1;}
    if($_=~/^SORT_mem=(.*?);/){$SORT_mem=$1;}
    if($_=~/^SORT_time=(.*?);/){$SORT_time=$1;}

    if($_=~/^CKBAM_job_name=(.*?);/){$CKBAM_job_name=$1;}
    if($_=~/^CKBAM_nodes=(.*?);/){$CKBAM_nodes=$1;}
    if($_=~/^CKBAM_cpus_per_task=(.*?);/){$CKBAM_cpus_per_task=$1;}
    if($_=~/^CKBAM_mem=(.*?);/){$CKBAM_mem=$1;}
    if($_=~/^CKBAM_time=(.*?);/){$CKBAM_time=$1;}
    if($_=~/^GNUdir=(.*?);/){$GNUdir=$1;}
    if($_=~/^CoverageCheck=(.*?);/){$coveragecheck=$1;}

    ## bam mkdup 
    if($_=~/^MkupCheck=(.*?);/){$mkdupcheck=$1;}
    if($_=~/^picard=(.*?);/){$picard=$1;}
    if($_=~/^mkdupPara=(.*?);/){$mkdupPara=$1;}
    if($_=~/^MKDUP_job_name=(.*?);/){$MKDUP_job_name=$1;}
    if($_=~/^MKDUP_nodes=(.*?);/){$MKDUP_nodes=$1;}
    if($_=~/^MKDUP_cpus_per_task=(.*?);/){$MKDUP_cpus_per_task=$1;}
    if($_=~/^MKDUP_mem=(.*?);/){$MKDUP_mem=$1;}
    if($_=~/^MKDUP_time=(.*?);/){$MKDUP_time=$1;}

    if($_=~/^CKDUP_job_name=(.*?);/){$CKDUP_job_name=$1;}
    if($_=~/^CKDUP_nodes=(.*?);/){$CKDUP_nodes=$1;}
    if($_=~/^CKDUP_cpus_per_task=(.*?);/){$CKDUP_cpus_per_task=$1;}
    if($_=~/^CKDUP_mem=(.*?);/){$CKDUP_mem=$1;}
    if($_=~/^CKDUP_time=(.*?);/){$CKDUP_time=$1;}

    if($_=~/^COVER_job_name=(.*?);/){$COVER_job_name=$1;}
    if($_=~/^COVER_nodes=(.*?);/){$COVER_nodes=$1;}
    if($_=~/^COVER_cpus_per_task=(.*?);/){$COVER_cpus_per_task=$1;}
    if($_=~/^COVER_mem=(.*?);/){$COVER_mem=$1;}
    if($_=~/^COVER_time=(.*?);/){$COVER_time=$1;}

    ##deepvariant
    if($_=~/^run_deepvariant=(.*?);/){$run_deepvariant=$1;}
    if($_=~/^deepvariant_bam=(.*?);/){$deepvariant_bam=$1;}
    if($_=~/^container=(.*?);/){$container=$1;}
    if($_=~/^deepvariant=(.*?);/){$deepvariant=$1;}
    if($_=~/^deepvariant_SIF=(.*?);/){$deepvariant_SIF=$1;}
    if($_=~/^deepvariant_version=(.*?);/){$deepvariant_version=$1;}
    if($_=~/^deepvariant_mode=(.*?);/){$deepvariant_mode=$1;}
    if($_=~/^deepvariant_bed=(.*?);/){$deepvariant_bed=$1;}
    if($_=~/^deepvariant_core=(.*?);/){$deepvariant_core=$1;}
    if($_=~/^deepvariant_mem=(.*?);/){$deepvariant_mem=$1;}
    if($_=~/^gpu=(.*?);/){$gpu=$1;}
    
    if($_=~/^DEEP_job_name=(.*?);/){$DEEP_job_name=$1;}
    if($_=~/^DEEP_nodes=(.*?);/){$DEEP_nodes=$1;}
    if($_=~/^DEEP_cpus_per_task=(.*?);/){$DEEP_cpus_per_task=$1;}
    if($_=~/^DEEP_mem=(.*?);/){$DEEP_mem=$1;}
    if($_=~/^DEEP_time=(.*?);/){$DEEP_time=$1;}

    ## other tools
    if($_=~/^java8=(.*?);/){$java8=$1;}
}

################################### main parameters #####################################
$BIN ||= "/opt/pipeline/bin";
$DATAPATH ||= "/opt/pipeline/database";
$LIBPATH ||= "/opt/pipeline/lib";
#$BIN ||= "/home/cheng/Desktop/t2/bin";
#$DATAPATH ||= "/home/cheng/Desktop/t2/database";
#$LIBPATH ||= "/home/cheng/Desktop/t2/lib";
die "The bin directory $BIN does not exist!\n" if(!-e $BIN);
die "The database directory $DATAPATH does not exist!\n" if(!-e $DATAPATH);
die "The library directory $LIBPATH does not exist!\n" if(!-e $LIBPATH);

$reference ||="$DATAPATH/reference/hg19/hg19.fa";
die "The reference $reference does not exist!\n" if(!-e $reference);

## tools
$SOAPnuke ||="$BIN/SOAPnuke1.5.6";
$fqcheck ||="$BIN/fqcheck33";
$bwa ||="$BIN/bwa";
$samtools ||="$BIN/samtools"; 
$picard ||="$BIN/picard-2.10.9.jar";
$java8 ||="$BIN/java";
$container ||="$BIN/singularity";
$GNUdir ||= "$BIN/gnuplot-4.6.7/PostScript";

## clean ,bwa-mem and picard parameters
$SOAPnukePara ||= "-l 5 -q 0.5 -n 0.1 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG";
$bwaMemPara ||= "-t 4 -k 32";
$mkdupPara ||= "REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true PROGRAM_RECORD_ID=null";

## deepvariant
$run_deepvariant ||= "true";
$deepvariant ||= "/opt/deepvariant/bin/run_deepvariant";
$deepvariant_version ||= "0.10.0";
$deepvariant_mode ||= "WGS";
$deepvariant_bed ||= "";
$deepvariant_core ||= 4;
$deepvariant_mem ||= "$DATAPATH/singularity_resource.toml";
$deepvariant_SIF ||= "$BIN/deepvariant_0.10.0.sif";
$deepvariant_bam ||= "mkdup";
$gpu ||= "false";

## others
$center ||= "BGI";
$platform ||= "BGISEQ-500";
$bwacheck ||= "true";
$sortcheck ||= "true";
$bamcheck ||= "true";
$mkdupcheck ||= "false";
$coveragecheck ||= "true";
############################### unify chromosome names ####################################

my @chrs=();
for (my $i=1; $i<=22; $i++) {push @chrs,"chr$i";}
push @chrs,"chrX";
push @chrs,"chrY";
push @chrs,"chrM";

################################ Process sequence ##########################################

open LL,$List or die $!;
`mkdir -p $outdir/QCcheck`;
open INFOLIST,">$outdir/QCcheck/result.info.list" or die $!;
while(<LL>){
    chomp;
    my($q,$id);
    my ($samp,$qual,$readsId,$lib,$lane,$lane_barcode,$fqprefix)=split /\s+/, $_;
    if($qual==33){ $q=2; }else{ $q=1; }
    if($readsId =~ /old/i){ $id=0 ; }elsif($readsId =~ /new/i){ $id=1 ; }else{ die "please offer reads id : old or new !";}
    `mkdir -p $outdir/$samp`;
    `mkdir -p $outdir/$samp/$lane_barcode`;
    `mkdir -p $outdir/$samp/$lane_barcode/clean`;
    `mkdir -p $outdir/$samp/$lane_barcode/bwa`;
    `mkdir -p $outdir/$samp/$lane_barcode/shell`;
    `mkdir -p $outdir/$samp/$lane_barcode/tmp`;
    $LIBS{$samp}{$lane_barcode}=1;
    &peLaneClean($samp,$fqprefix,$lane_barcode,$q,$id);
    &Alignment($samp,$lane,$lib,$lane_barcode);
    &GATKbrecal($samp,$lane_barcode,$lib);
    print INFOLIST "$samp\t$lib\t$lane_barcode\t$outdir/$samp/$lane_barcode/clean/$lane_barcode.clean.stat\t$outdir/$samp/$lane_barcode/bwa/$lane_barcode.stat\t$outdir/$samp/mkdup-realn/$samp.mkdup.stat\t$outdir/$samp/Coverage/information.xls\n";
}
close LL;
close INFOLIST;

&MergeBam(); ## Merge sorted bams
if ($deepvariant_bam =~ /mkdup/i) {&Mkdup()}; ##  Marking duplications
if ($coveragecheck=~/t/i) {&coverageAdepth()}; ## calculate coverage and depth information
if ($run_deepvariant=~/t/i) {&deepVariant()}; ## SNP calling by deepvariant

## check bam list for the pipline
if($bamcheck=~/true/i){
    open INFOLIST,">$outdir/QCcheck/check.bam.list" or die $!;
    foreach my $samp (sort keys %samples){
        print INFOLIST "$samp\t";
        foreach my $lib (sort keys %{$samples{$samp}}){
            foreach my $lane_barcode (sort keys %{$samples{$samp}{$lib}}) {
                print INFOLIST "$outdir/$samp/$lane_barcode/shell/clean_Check.txt\t$outdir/$samp/$lane_barcode/shell/bwa_Check.txt\t$outdir/$samp/$lane_barcode/shell/brecal_Check.txt\t";
            }
        }
        print INFOLIST "$outdir/$samp/shell/merge_Check.txt\t$outdir/$samp/shell/mkdup_Check.txt\t";
    }
    close INFOLIST;
}


############################### monitor to control the pipeline #######################

&edgeList();  ##creat scrips to monitor

#########################################################################################

## date
sub date{
    my $dat='';
    my $tmp=`date`;
    my @tmp=split /\s+/,$tmp;
    $dat.=$tmp[2];
    my @time=split /:/,$tmp[3];
    foreach (@time){$dat.=$_;}
    return $dat;
}

## echo sign for monitor  
sub echostring{
    my $sh=shift;
    my $ostr="echo ==========end at : `date` ========== && \\\n";
    $ostr.="echo Still_waters_run_deep 1>&2 && \\\n";
    $ostr.="echo Still_waters_run_deep > $sh.sign\n";
    return $ostr;
}

## fastq Clean by SOAPnuke
sub peLaneClean{
    my ($samp,$fqprefix,$lane_barcode,$q,$id)=@_;
    my ($fq1,$fq2);
    my @fqfiles=`ls $fqprefix\_*.fq*`;
    foreach(@fqfiles){
        chomp;
        if($_=~/\_1.fq$/ || $_=~/\_1.fq.gz$/){$fq1=$_;}
        if($_=~/\_2.fq$/ || $_=~/\_2.fq.gz$/){$fq2=$_;}
    }
    die "Do not exists $fqprefix*.fq.gz\n" unless ($fq1 && $fq2);

    my $cleandir="$outdir/$samp/$lane_barcode/clean";
    my $fqout1="$cleandir/$lane_barcode\_1.clean.fq.gz";
    my $fqout2="$cleandir/$lane_barcode\_2.clean.fq.gz";
    my $cleanCmd="#!/bin/bash\n";
    $cleanCmd.="#SBATCH --job-name=$CLEAN_job_name\n#SBATCH --time=$CLEAN_time\n#SBATCH --mem=$CLEAN_mem\n#SBATCH --cpus-per-task=$CLEAN_cpus_per_task\n#SBATCH --nodes=$CLEAN_nodes\necho ==========start at : `date` ========== && \\\n";
    $cleanCmd.="if [ -e \"$fqout1\" ];then rm -rf $fqout1;fi && \\\n";
    $cleanCmd.="if [ -e \"$fqout2\" ];then rm -rf $fqout2;fi && \\\n";
    $cleanCmd.="$SOAPnuke filter -1 $fq1 -2 $fq2";
    $cleanCmd.=" $SOAPnukePara -Q $q -G --seqType $id -o $cleandir -C $fqout1 -D $fqout2 && \\\n"; 
    $cleanCmd.="gzip -t $fqout1 $fqout2 && \\\n";
    $cleanCmd.="perl $BIN/soapnuke_stat.pl $cleandir/Basic_Statistics_of_Sequencing_Quality.txt $cleandir/Statistics_of_Filtered_Reads.txt > $cleandir/$lane_barcode.clean.stat && \\\n";
    $cleanCmd.="rm -rf $cleandir/*.txt && \\\n";
    $cleanCmd.="$fqcheck -r $fqout1 -c $cleandir/$lane_barcode\_1.clean.fqcheck && \\\n";
    $cleanCmd.="$fqcheck -r $fqout2 -c $cleandir/$lane_barcode\_2.clean.fqcheck && \\\n";
    if ($GNUdir){$cleanCmd.="export GNUPLOT_PS_DIR=$GNUdir && \\\n";}
    $cleanCmd.="perl $BIN/fqcheck_distribute.pl $cleandir/$lane_barcode\_1.clean.fqcheck $cleandir/$lane_barcode\_2.clean.fqcheck -o $cleandir/$lane_barcode.clean. && \\\n";
    open CSH,">$outdir/$samp/$lane_barcode/shell/clean.sh";
    $cleanCmd.=echostring("$outdir/$samp/$lane_barcode/shell/clean.sh");
    print CSH $cleanCmd;
    close CSH;
}

## Alignment by bwa mem
sub Alignment{
    my ($samp,$lane,$lib,$lane_barcode)=@_;
    my $rgid = $lane_barcode;
	$rgid=~s/_read//;
	my ($fq1,$fq2);
    my $cleandir="$outdir/$samp/$lane_barcode/clean";
    my $bwadir="$outdir/$samp/$lane_barcode/bwa";
    $fq1="$cleandir/$lane_barcode\_1.clean.fq.gz";
    $fq2="$cleandir/$lane_barcode\_2.clean.fq.gz";
    
    # mapping reads to reference
    open BMSH,">$outdir/$samp/$lane_barcode/shell/bwa_mem.sh" or die $!;
    print BMSH "#!/bin/bash\n";
    print BMSH "#SBATCH --job-name=$ALIGN_job_name\n#SBATCH --time=$ALIGN_time\n#SBATCH --mem=$ALIGN_mem\n#SBATCH --cpus-per-task=$ALIGN_cpus_per_task\n#SBATCH --nodes=$ALIGN_nodes\necho ==========start at : `date` ========== && \\\n";
    print BMSH "$bwa mem $bwaMemPara -M -R \'\@RG\\tID:$rgid\\tSM:$samp\\tLB:$lib\\tPU:$lane\\tPL:$platform\\tCN:$center\' $reference $fq1 $fq2 | $samtools view -b -S -F 256 -t $reference.fai -o $bwadir/$lane_barcode.bam - && \\\n";
    print BMSH echostring("$outdir/$samp/$lane_barcode/shell/bwa_mem.sh");
    close BMSH;
    
    # Validation for bwa
    if ($bwacheck=~/true/){
        open BMSH,">$outdir/$samp/$lane_barcode/shell/bwa_check.sh" or die $!;
        print BMSH "#!/bin/bash\n";
        print BMSH "#SBATCH --job-name=$CKBWA_job_name\n#SBATCH --time=$CKBWA_time\n#SBATCH --mem=$CKBWA_mem\n#SBATCH --cpus-per-task=$CKBWA_cpus_per_task\n#SBATCH --nodes=$CKBWA_nodes\necho ==========start at : `date` ========== && \\\n";
        print BMSH "$samtools flagstat $bwadir/$lane_barcode.bam >$bwadir/$lane_barcode.stat && \\\n";
        if($bamcheck=~/true/i){
            print BMSH "perl $BIN/bamCheck_v1.pl -a $bwadir/$lane_barcode.stat -b $bwadir/$lane_barcode.bam >$outdir/$samp/$lane_barcode/shell/clean_Check.txt && \\\n"; 
        }
        print BMSH "rm -rf $fq1 $fq2 && \\\n";
        print BMSH echostring("$outdir/$samp/$lane_barcode/shell/bwa_check.sh");
        close BMSH;
    }
    
    #sort bam file 
    open BMSH,">$outdir/$samp/$lane_barcode/shell/bam_sort.sh" or die $!;
    print BMSH "#!/bin/bash\n";
    print BMSH "#SBATCH --job-name=$SORT_job_name\n#SBATCH --time=$SORT_time\n#SBATCH --mem=$SORT_mem\n#SBATCH --cpus-per-task=$SORT_cpus_per_task\n#SBATCH --nodes=$SORT_nodes\necho ==========start at : `date` ========== && \\\n";
    print BMSH "$java8 -Xmx6g -XX:MaxMetaspaceSize=512m -XX:-UseGCOverheadLimit -XX:ParallelGCThreads=1 -jar $picard SortSam SORT_ORDER=coordinate I=$bwadir/$lane_barcode.bam O=$bwadir/$lane_barcode.sort.bam CREATE_INDEX=true && \\\nln -s $bwadir/$lane_barcode.sort.bai $bwadir/$lane_barcode.sort.bam.bai && \\\n";
    print BMSH "rm -rf $bwadir/tmp && \\\n";
    print BMSH echostring("$outdir/$samp/$lane_barcode/shell/bam_sort.sh");
    close BMSH;

    # Validation for sort
    if ($sortcheck=~/true/){
        open BMSH,">$outdir/$samp/$lane_barcode/shell/sort_check.sh" or die $!;
        print BMSH "#!/bin/bash\n";
        print BMSH "#SBATCH --job-name=$CKBAM_job_name\n#SBATCH --time=$CKBAM_time\n#SBATCH --mem=$CKBAM_mem\n#SBATCH --cpus-per-task=$CKBAM_cpus_per_task\n#SBATCH --nodes=$CKBAM_nodes\necho ==========start at : `date` ========== && \\\n";
        print BMSH "$samtools flagstat $bwadir/$lane_barcode.sort.bam >$bwadir/$lane_barcode.sort.stat && \\\n";
        if($bamcheck=~/true/i){
            print BMSH "perl $BIN/bamCheck_v1.pl -a $bwadir/$lane_barcode.bam -b $bwadir/$lane_barcode.sort.bam >$outdir/$samp/$lane_barcode/shell/bwa_Check.txt && \\\n";
        }else{
            print BMSH "rm -rf $bwadir/$lane_barcode.bam && \\\n";
        }
        print BMSH echostring("$outdir/$samp/$lane_barcode/shell/sort_check.sh");
        close BMSH;
    }
}

## Base Recalibratio by gatk
sub GATKbrecal{
    my ($samp,$lane_barcode,$lib)=@_;
    my $bwadir="$outdir/$samp/$lane_barcode/bwa";
    my $bam="$bwadir/$lane_barcode.sort.bam";
    $samples{$samp}{$lib}{$lane_barcode}="$bam"; 
}

sub MergeBam{
    foreach my $samp (sort keys %LIBS) {
	`mkdir -p $outdir/$samp/merge `;
	`mkdir -p $outdir/$samp/shell`;
	open STBM,">$outdir/$samp/shell/bam_merge.sh" or die $!;
	my @keys = %{$LIBS{$samp}};
	my $sortCmd;
	my $lane_barcode = $keys[0];
	my $bwadir = "$outdir/$samp/$lane_barcode/bwa";
	#print join "\t",@keys,"\n";
	if ($#keys > 1){
		$sortCmd= "#!/bin/bash\n";
		$sortCmd.= "#SBATCH --job-name=$SORT_job_name\n#SBATCH --time=$SORT_time\n#SBATCH --mem=$SORT_mem\n#SBATCH --cpus-per-task=$SORT_cpus_per_task\n#SBATCH --nodes=$SORT_nodes\necho ==========start at : `date` ========== && \\\n";
		$sortCmd.= "$java8 -Xmx6g -XX:MaxMetaspaceSize=512m -XX:-UseGCOverheadLimit -XX:ParallelGCThreads=1 -jar $picard MergeSamFiles AS=true";
		foreach my $lib (sort keys %{$samples{$samp}}){
			foreach my $lane_barcode (sort keys %{$samples{$samp}{$lib}}) {
				$sortCmd.=" I=$samples{$samp}{$lib}{$lane_barcode}";
			}
		}
		$sortCmd.= " O=$outdir/$samp/merge/$samp.sort.bam && \\\n$samtools index $outdir/$samp/merge/$samp.sort.bam && \\\n";
	}else{
		$sortCmd="#!/bin/bash\n";
		$sortCmd.= "#SBATCH --job-name=$SORT_job_name\n#SBATCH --time=$SORT_time\n#SBATCH --mem=128\n#SBATCH --cpus-per-task=1\n#SBATCH --nodes=1\necho ==========start at : `date` ========== && \\\n";
		$sortCmd.="ln -s $bwadir/$lane_barcode.sort.bam $outdir/$samp/merge/$samp.sort.bam && \\\nln -s $bwadir/$lane_barcode.sort.bam.bai $outdir/$samp/merge/$samp.sort.bam.bai && \\\n";
	}
	print STBM "$sortCmd";
	print STBM echostring("$outdir/$samp/shell/bam_merge.sh");
	close STBM;
        $SortBAM{$samp}="$outdir/$samp/merge/$samp.sort.bam";
    }
}

## mark duplication
sub Mkdup{
    foreach my $samp (sort keys %samples) {
        `mkdir -p $outdir/$samp/mkdup-realn`;
        `mkdir -p $outdir/$samp/tmp`;        
        my $brecalBam='';
        # use picard to mark duplication 
        open DUPSH,">$outdir/$samp/shell/mkdup.$samp.sh" or die $!;
        my $mkdupCmd="#!/bin/bash\n";
	$mkdupCmd.="#SBATCH --job-name=$MKDUP_job_name\n#SBATCH --time=$MKDUP_time\n#SBATCH --mem=$MKDUP_mem\n#SBATCH --cpus-per-task=$MKDUP_cpus_per_task\n#SBATCH --nodes=$MKDUP_nodes\necho ==========start at : `date` ========== && \\\n";
        $mkdupCmd.="$java8 -Xmx6g -XX:MaxMetaspaceSize=512m -XX:-UseGCOverheadLimit -XX:ParallelGCThreads=1 -jar $picard MarkDuplicates CREATE_INDEX=true ";
        $mkdupCmd.="I=$outdir/$samp/merge/$samp.sort.bam O=$outdir/$samp/mkdup-realn/$samp.mkdup.bam METRICS_FILE=$outdir/$samp/mkdup-realn/$samp.mkdup.bam.met TMP_DIR=$outdir/$samp/tmp $mkdupPara && \\\n";
        $mkdupCmd.="rm -rf $outdir/$samp/tmp && \\\n";
	$mkdupCmd.=echostring("$outdir/$samp/shell/mkdup.$samp.sh");
        print DUPSH $mkdupCmd;
        close DUPSH;

        # Validation for mark duplication
	if ($mkdupcheck=~/true/){
            chop($brecalBam);
            open DUPSH,">$outdir/$samp/shell/mkdup_check.$samp.sh" or die $!;
            print DUPSH "#!/bin/bash\n";
            print DUPSH "#SBATCH --job-name=$CKDUP_job_name\n#SBATCH --time=$CKDUP_time\n#SBATCH --mem=$CKDUP_mem\n#SBATCH --cpus-per-task=$CKDUP_cpus_per_task\n#SBATCH --nodes=$CKDUP_nodes\necho ==========start at : `date` ========== && \\\n";
	    print DUPSH "$samtools flagstat $outdir/$samp/mkdup-realn/$samp.mkdup.bam >$outdir/$samp/mkdup-realn/$samp.mkdup.stat &&\\\n";
	        if($bamcheck=~/true/i){
	            open REMOVE,">$outdir/$samp/mkdup-realn/mkdup.bam.list" or die $!; ## bam list for check and remove
	            my @F=split/\s+/,$brecalBam;
	            foreach(@F){ print REMOVE "$_\n"; }
	            close REMOVE;
                print DUPSH "perl $BIN/bam2Check.pl $outdir/$samp/mkdup-realn/mkdup.bam.list $outdir/$samp/mkdup-realn/$samp.mkdup.bam >$outdir/$samp/shell/mkdup_check.txt && \\\n";
                print DUPSH "rm -rf $outdir/$samp/mkdup-realn/mkdup.bam.list && \\\n";
	        }else{ 
                print DUPSH "rm -rf $brecalBam && \\\n"; 
                }
            print DUPSH echostring("$outdir/$samp/shell/mkdup_check.$samp.sh");
            close DUPSH;
	}
    }
}

## calculate the coverage and depth information for each sample
sub coverageAdepth{
    foreach my $samp (sort keys %SortBAM) { 
        my $bam=$SortBAM{$samp};
        `mkdir -p $outdir/$samp/Coverage`;
        open COVDEPSH,">$outdir/$samp/shell/covdep.$samp.sh" or die $!;
        print COVDEPSH "#!/bin/bash\n";
        print COVDEPSH "#SBATCH --job-name=$COVER_job_name\n#SBATCH --time=$COVER_time\n#SBATCH --mem=$COVER_mem\n#SBATCH --cpus-per-task=$COVER_cpus_per_task\n#SBATCH --nodes=$COVER_nodes\necho ==========start at : `date` ========== && \\\n";
        print COVDEPSH "perl $BIN/picard.pl $bam >$outdir/$samp/Coverage/information.xls && \\\n";
        print COVDEPSH "perl $BIN/depthV2.0.pl $bam $outdir/$samp/Coverage >>$outdir/$samp/Coverage/information.xls && \\\n";
        print COVDEPSH echostring("$outdir/$samp/shell/covdep.$samp.sh");
        close COVDEPSH;
    }
}
## generate the execution script 
sub edgeList{
    `mkdir -p $outdir/shell_run/`;
    open EDGE1,">$outdir/shell_run/clean.sh" or die $!;
    open EDGE2,">$outdir/shell_run/bwa.sh" or die $!;
    if ($bwacheck=~/true/){ open EDGE3,">$outdir/shell_run/bwaCheck.sh" or die $!;}
    open EDGE4,">$outdir/shell_run/sortBam.sh" or die $!;
    open EDGE10,">$outdir/shell_run/mergeBam.sh" or die $!;
    if ($sortcheck=~/true/) { open EDGE5,">$outdir/shell_run/sortCheck.sh" or die $!;}
    if ($deepvariant_bam=~/mkdup/) {open EDGE6,">$outdir/shell_run/mkdup.sh" or die $!;}
    if ($mkdupcheck=~/true/){open EDGE7,">$outdir/shell_run/mkdupCheck.sh" or die $!;}
    open EDGE8,">$outdir/shell_run/covdep.sh" or die $!;
    open EDGE9,">$outdir/shell_run/deepvariant.sh" or die $!;
    my %flagsh;
	
    foreach my $samp(sort keys %samples){
        ## sequence processing
        foreach my $lib(sort keys %{$samples{$samp}}){
            foreach my $lane_barcode(sort keys %{$samples{$samp}{$lib}}){
                print EDGE1 "sh $outdir/$samp/$lane_barcode/shell/clean.sh\n";
                print EDGE2 "sh $outdir/$samp/$lane_barcode/shell/bwa_mem.sh\n";
                if ($bwacheck=~/true/) {print EDGE3 "sh $outdir/$samp/$lane_barcode/shell/bwa_check.sh\n";}
                print EDGE4 "sh $outdir/$samp/$lane_barcode/shell/bam_sort.sh\n";
                if ($sortcheck=~/true/){print EDGE5 "sh $outdir/$samp/$lane_barcode/shell/sort_check.sh\n";}
            }
        }	
		print EDGE10 "sh $outdir/$samp/shell/bam_merge.sh\n";
		if ($deepvariant_bam=~/mkdup/) {print EDGE6 "sh $outdir/$samp/shell/mkdup.$samp.sh\n";}
		if ($mkdupcheck=~/true/) {print EDGE7 "sh $outdir/$samp/shell/mkdup_check.$samp.sh\n";}
		print EDGE8 "sh $outdir/$samp/shell/covdep.$samp.sh\n";
		print EDGE9 "sh $outdir/$samp/shell/deepVariant.$samp.sh\n";
    }
    close EDGE1;close EDGE2;close EDGE4;close EDGE8;close EDGE9;close EDGE10;
    if ($bwacheck=~/true/){close EDGE3;}
    if ($sortcheck=~/true/){close EDGE5;}
    if ($deepvariant_bam=~/mkdup/){close EDGE6;}
    if ($mkdupcheck=~/true/){close EDGE7;}

    my @r = ('a'..'z','A'..'Z');
    $FileName ||= join '', map { $r[int rand @r] } 0..6;
    open RUNSH,">$outdir/shell_run/run.$FileName.sh" or die $!;
    print RUNSH "## cleaning raw fq files\nsh $outdir/shell_run/clean.sh\n\n";
    print RUNSH "## aligning using BWA\nsh $outdir/shell_run/bwa.sh\n\n";
    if ($bwacheck=~/true/){print RUNSH "## checking the results of BWA\nsh $outdir/shell_run/bwaCheck.sh &\n\n";}
    print RUNSH "## sorting bam\nsh $outdir/shell_run/sortBam.sh\n\n";
    if ($sortcheck=~/true/){print RUNSH "## checking sorted BAM\nsh $outdir/shell_run/sortCheck.sh &\n\n";}
    print RUNSH "## merging bam\nsh $outdir/shell_run/mergeBam.sh\n\n";
    if ($deepvariant_bam=~/mkdup/){print RUNSH "## marking duplicates in BAMs\nsh $outdir/shell_run/mkdup.sh\n\n";}
    print RUNSH "## Bam coverage statistics\nsh $outdir/shell_run/covdep.sh \n\n";
    if ($mkdupcheck=~/true/){ print RUNSH "## checking duplicates-marked BAM\nsh $outdir/shell_run/mkdupCheck.sh \n\n";}
    print RUNSH "## SNP calling by Deepvariant\n#sh $outdir/shell_run/deepvariant.sh\n";  
    close RUNSH;
}

sub deepVariant{
	foreach my $samp (sort keys %samples){
		open DV, ">$outdir/$samp/shell/deepVariant.$samp.sh" or die $!;
		`mkdir -p $outdir/$samp/SNP`;
		my ($ref_dir,$ref_fa)=$reference=~/(.*)\/(.*)/;
		my $bamdir = $deepvariant_bam=~/sort/ ? "merge" : "mkdup-realn";
		my $dp_bam = $bamdir =~ /merge/ ? "$samp.sort.bam" : "$samp.mkdup.bam";
		my $BED = ($deepvariant_bed) ? "--regions=$deepvariant_bed" : "";
		print DV "#!/bin/bash\n";
		# cpu version
		if ($gpu=~/false/i){
			print DV "#SBATCH --job-name=$DEEP_job_name\n#SBATCH --time=$DEEP_time\n#SBATCH --mem=$DEEP_mem\n#SBATCH --cpus-per-task=$DEEP_cpus_per_task\n#SBATCH --nodes=$DEEP_nodes\necho ==========start at : `date` ========== && \\\n";
			print DV "INPUT_DIR=\"$outdir/$samp/$bamdir\" && \\\nOUTPUT_DIR=\"$outdir/$samp/SNP\" && \\\nREF_DIR=\"$ref_dir\" && \\\n";
			if ($container=~/singularity/i){
				print DV "#sudo $container pull docker://google/deepvariant:\"$deepvariant_version\" && \\\n";
				print DV "sudo $container run --apply-cgroups $deepvariant_mem -B \${INPUT_DIR}:/input -B \${OUTPUT_DIR}:/output -B \${REF_DIR}:/reference  $deepvariant_SIF $deepvariant --model_type=$deepvariant_mode --ref=/reference/$ref_fa --reads=/input/$dp_bam $BED --output_vcf=/output/$samp.vcf.gz --output_gvcf=/output/$samp.g.vcf.gz --num_shards=$deepvariant_core && \\\n";
			}else{
				my $docker_mem=$DEEP_mem."M";
				print DV "service docker start && \\\n";
				print DV "sudo $container pull google/deepvariant:\"$deepvariant_version\" && \\\n";
				print DV "sudo $container run -m $docker_mem -c $deepvariant_core  --privileged=true -v \"\${INPUT_DIR}\":\"/input\":rw -v \"\${OUTPUT_DIR}\":\"/output\":rw -v \"\${REF_DIR}\":\"/reference\":rw  google/deepvariant:\"$deepvariant_version\" $deepvariant --model_type=$deepvariant_mode --ref=/reference/$ref_fa --reads=/input/$dp_bam $BED --output_vcf=/output/$samp.vcf.gz --output_gvcf=/output/$samp.g.vcf.gz --num_shards=$deepvariant_core && \\\n";
			}
		}else{ ## GPU version
			print DV "#SBATCH --job-name=$DEEP_job_name\n#SBATCH --time=$DEEP_time\n#SBATCH --mem=$DEEP_mem\n#SBATCH --cpus-per-task=$DEEP_cpus_per_task\n#SBATCH --nodes=$DEEP_nodes\necho ==========start at : `date` ========== && \\\n";
			print DV "INPUT_DIR=\"$outdir/$samp/$bamdir\" && \\\nOUTPUT_DIR=\"$outdir/$samp/SNP\" && \\\nREF_DIR=\"$ref_dir\" && \\\n";
			if ($container=~/singularity/i){
				print DV "#sudo $container pull docker://google/deepvariant:\"$deepvariant_version-gpu\" && \\\n";
				print DV "sudo $container run --nv --apply-cgroups $deepvariant_mem -B \${INPUT_DIR}:/input -B \${OUTPUT_DIR}:/output -B \${REF_DIR}:/reference  $deepvariant_SIF [ =docker://google/deepvariant:\"\${BIN_VERSION}-gpu\" ] $deepvariant --model_type=$deepvariant_mode --ref=/reference/$ref_fa --reads=/input/$dp_bam $BED --output_vcf=/output/$samp.vcf.gz --output_gvcf=/output/$samp.g.vcf.gz --num_shards=$deepvariant_core && \\\n";
			}else{
				my $docker_mem=$DEEP_mem."M";
				print DV "service $container start && \\\n";
				print DV "sudo $container pull google/deepvariant:\"$deepvariant_version-gpu\" && \\\n";
				print DV "sudo $container run -m $docker_mem -c $deepvariant_core  --privileged=true -v \"\${INPUT_DIR}\":\"/input\":rw -v \"\${OUTPUT_DIR}\":\"/output\":rw -v \"\${REF_DIR}\":\"/reference\":rw  google/deepvariant:\"$deepvariant_version-gpu\" $deepvariant --model_type=$deepvariant_mode --ref=/reference/$ref_fa --reads=/input/$dp_bam $BED --output_vcf=/output/$samp.vcf.gz --output_gvcf=/output/$samp.g.vcf.gz --num_shards=$deepvariant_core && \\\n";
			}
		}
		print DV echostring("$outdir/$samp/shell/deepVariant.$samp.sh");
		close DV;
	}
}
print "Please execute: \"sh $outdir/shell_run/run.$FileName.sh\" to start the pipeline.\n";

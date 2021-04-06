#!/usr/bin/perl

=Usage
  
  perl SLURM.pl --lines 1 --maxjob 100 --interval 600 your_work_shell.sh 
  --interval <num>  set interval time of checking by squeue, default 300 seconds
  --lines <num>     set number of lines to form a job, default 1
  --maxjob <num>    set the maximum number of jobs to run in parallel, default 50
  --help            output help information to screen  
=Usage
  
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

##get options from command line into variables and set default values
my ($Interval,$Lines, $Maxjob, $Help);
GetOptions(
	"lines:i"=>\$Lines,
	"maxjob:i"=>\$Maxjob,
	"interval:i"=>\$Interval,
	"help"=>\$Help
);
$Interval ||= 300;
$Lines ||= 1;
$Maxjob ||= 50;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $work_shell_file = shift;

##global variables
my $work_shell_file_globle = $work_shell_file.".$$.globle";
my $work_shell_file_error = $work_shell_file.".$$.log";
my $Work_dir = $work_shell_file.".$$.qsub";
my $current_dir = `pwd`; chomp $current_dir;

$work_shell_file_globle = $work_shell_file;

## read from input file, make the slurm shell files
my $line_mark = 0;
my $Job_mark="00001";
mkdir($Work_dir);
my @Shell;  ## store the file names of sbatch shell
open IN, $work_shell_file_globle || die "fail to open $work_shell_file_globle";
my $line_num='0000';
while(<IN>){
	chomp;
	s/&/;/g;
	next unless($_);
	if ($line_mark % $Lines == 0) {
		open OUT,">$Work_dir/$Job_prefix\_$Job_mark.sh" || die "failed creat $Job_prefix\_$Job_mark.sh";
		push @Shell,"$Job_prefix\_$Job_mark.sh";
		$Job_mark++;
		$line_num='0000';
	}
	s/;\s*$//;  ##delete the last character ";", because two ";;" characters will cause error in qsub
	s/;\s*;/;/g;
	$line_num++;
	print OUT $_."; echo $line_num is OK!\n";
#	print OUT $_."; echo This-Work-is-Completed!\n";

	if ($line_mark % $Lines == $Lines - 1) {
		print OUT "echo This-Work-is-Completed at `date`\n";
		close OUT;
	}
	
	$line_mark++;
}
close IN;
print OUT "echo This-Work-is-Completed!\n";
close OUT;


## run jobs by sbatch, until all the jobs are really finished
my $qsub_cycle = 1;
while (@Shell) {
	my %Alljob; ## store all the job IDs of this cycle
	my %Runjob; ## store the real running job IDs of this cycle
	my %Error;  ## store the unfinished jobs of this cycle
	chdir($Work_dir); ##enter into the qsub working directoy
	my $job_cmd = "sbatch "; 
	
	for  (my $i=0; $i<@Shell; $i++) {
		 while (1) {
			my $run_num = run_count(\%Alljob,\%Runjob);
			if ($i < $Maxjob || ($run_num != -1 && $run_num < $Maxjob) ) {
				print "$job_cmd $Shell[$i]\n"; 
				my $jod_return = `$job_cmd $Shell[$i]`;
				my $job_id = $1 if($jod_return =~ /Your job (\d+)/);
				$Alljob{$job_id} = $Shell[$i];  ## job id => shell file name
				last;
			}else{
				sleep $Interval;
			}
		}
	}

	chdir($current_dir); ##return into original directory 


	###waiting for all jobs fininshed
	while (1) {
		my $run_num = run_count(\%Alljob,\%Runjob);	
		last if($run_num == 0);
		sleep $Interval;
	}

	##run the secure mechanism to make sure all the jobs are really completed
	open OUT, ">>$work_shell_file_error" || die "fail create $$work_shell_file_error";
	chdir($Work_dir); ##enter into the qsub working directoy
	foreach my $job_id (sort keys %Alljob) {
		my $shell_file = $Alljob{$job_id};
		
		##read the .out file
		my $content;
		if (-f "$shell_file.$job_id.out") {
			open IN,"$shell_file.$job_id.out" || warn "fail to read $shell_file.$job_id.out";
			$content = join("",<IN>);
			close IN;
		}
		##check whether the job has been killed during running time
		if ($content !~ /This-Work-is-Completed/) {
			$Error{$job_id} = $shell_file;
			print OUT "In cycle $qsub_cycle, In $shell_file.$job_id.out,  \"This-Work-is-Completed\" is not found, so this work may be unfinished\n";
		}
		

	}

	##make @shell for next cycle, which contains unfinished tasks
	@Shell = ();
	foreach my $job_id (sort keys %Error) {
		my $shell_file = $Error{$job_id};
		push @Shell,$shell_file;
	}
	
	$qsub_cycle++;
	if($qsub_cycle > 10000){
		print OUT "\n\nProgram stopped because the reqsub cycle number has reached 10000, the following jobs unfinished:\n";
		foreach my $job_id (sort keys %Error) {
			my $shell_file = $Error{$job_id};
			print OUT $shell_file."\n";
		}
		print OUT "Please check carefully for what errors happen, and redo the work, good luck!";
		die "\nProgram stopped because the reqsub cycle number has reached 10000\n";
	}
	
	print OUT "All jobs finished!\n" unless(@Shell);

	chdir($current_dir); ##return into original directory 
	close OUT;
}


####################################################
################### Sub Routines ###################
####################################################


##get the IDs and count the number of running jobs
##the All job list and user id are used to make sure that the job id belongs to this program 
##add a function to detect jobs on the died computing nodes.
sub run_count {
	my $all_p = shift;
	my $run_p = shift;
	my $run_num = 0;

	%$run_p = ();
	my $user = `whoami`; chomp $user;
	my $qstat_result = `squeue`;
	if ($qstat_result =~ /failed receiving gdi request/) {
		$run_num = -1;
		return $run_num; 
	}
	my @jobs = split /\n/,$qstat_result; 
	foreach my $job_line (@jobs) {
		$job_line =~s/^\s+//;
		my @job_field = split /\s+/,$job_line;
		next if($job_field[3] ne $user);
	}

	return $run_num; 
}


#!/usr/bin/perl

### Version 1.0 -- 14 April 2015 -- Alex Buerkle

### wrap_slurm_moran.pl : a perl script to submit many serial slurm
### jobs to the queue using the sbatch command.

use warnings;
use strict;

### CHANGES -----------------------------------------------

### v1.0 - Alex - update of wrap_qsub to work with the SLURM system

## possible future addition: specify memory requirement per node, or other hardware requirements


### ------------------ JOB CONFIGURATION --------------------------------------------------------------- 
### all variables that typically need to be modified are in the
### following block

my $arccproject='microbiome';
my $runtime = '1-00'; ## 7-00 is 7 days (days-hours notation), 1:00:00 is 1 hour

# specify a reasonable value here.  If the job does not finish by the
# time reached, the job is terminated.  Your job should get greater
# priority if you specify a minimal walltime
my $sendmail = 'FALSE'; ##  if TRUE, sends an email when each job completes

## build array of jobs to be run individually (serially) by pbs
my @jobarray = ();

#Example job arrays
# my $bindir = '/project/$arccproject/bin/';
# my $filedir = '/project/$arccproject/tparchma/loxia_afreq_all220/';
#foreach my $j (0..5){
#	push @jobarray, "$bindir"."entropy-intel -i $filedir"."$infile -o out_lox_entropy"."$j".".hdf5 -l 1001 -b 1 -t 1 -k 3 -m 1 -w 0\n";
#}
#foreach my $infile (@ARGV) {	
#    push @jobarray, "/project/$arccproject/src/popmod-29july14/popmod -i $filedir"."$infile  -n 20000 -t 2 -b 6000 -o out_"."$infile\n";
#}
##### simple job array using gsl-randist, good for debugging
# my $outfile;
# foreach my $j (0..2){
#     $outfile = "test-gsl-out-"."$j".".txt";
#     push @jobarray, "gsl-randist 19187 10000 binomial 0.5 25 > $outfile\n";
# }
foreach my $j (1..289){
	push @jobarray, "./scripts/randomForest_eachTaxon.py $j 100\n";
}


#my $thin = 10;
#my $mcmcsteps = 60000;
#my $burnin = 50000;
#my $filedir = '/project/EvolutionaryGenetics/buerkle/barnswallows/entropy_june2015/input/';

#my $rand = int(rand(1000000));  ### entropy can take a random number seed

#foreach my $k(1..6){
 #   foreach my $chain(1..4){
#	foreach my $set ('global', 'east', 'west'){
#	    ##	my $set = 'global';
#	    $rand = int(rand(1000000));  ### entropy can take a random number seed
#	    my $job = "/project/$arccproject/bin/entropy -m 1 -r $rand ";
#	    $job .="-i $filedir"."variants_0.8_common_subset.entropyinput_".$set."_focallociGT6.mpgl ";
#	    $job .= "-q $filedir"."ldak".$k."_".$set.".txt ";
#	    $job .="-t $thin -b $burnin -l $mcmcsteps -k $k -w 0 -s 10 ";
#	    $job .="-o $set"."_GT6_k".$k."_chain".$chain.".hdf5\n";
#	    push @jobarray, $job, 
#	}
 #   }	
#}


# foreach(@jobarray){
#     print $_;
# }
# exit;

### -------------------------END JOB CONFIGURATION---------------------------------------------------------
### UNLIKELY you need to change anything beyond here
### automatically configured or infrequently changed variables
my $emailaddress = $ENV{USER}.'\@uwyo.edu'; ### mail address for notification that a job has completed 
my $jobname = 'slurm.'.$ENV{USER};

### modules to load:
my @modules =();
#push @modules, 'module load intel/14.0.0';
#push @modules, 'module load gsl/1.16';

#these have to be loaded first on teton for most other modules
#push @modules, 'module load swset/2018.05';
#push @modules, 'module load gcc/7.3.0';
#push @modules, 'module load gobject-introspection/1.49.2-py27';
#push @modules, 'module load harfbuzz';
#push @modules, 'module load pango/1.41.0-py27';
#push @modules, 'module load cairo/1.14.12-py27'; 
#push @modules, 'module load r/3.5.0-py27';
#push @modules, 'module load r/3.4.4';
#push @modules, 'module load r-rjags/4-6';
#push @modules, 'module load r-coda/0.19-1';
#push @modules, 'module load r-lattice/0.20-35';
push @modules, 'module load  python/3.8.7';

my $workdir = '/lscratch'; 
## Set to /lscratch on mtmoran.  This is local disk with high
## performance.  This is where results will be written temporarily on
## a node. 
my $basestoredir = "/project/$arccproject/users/$ENV{USER}";
## use a directory in /project/$arccproject/ as this has more
## allocated disk space than your home directory.  Will write to
## slurm_log and slurm_results inside this directory
my $logdir = "$basestoredir/slurm_log";
my $resultdir = "$basestoredir/slurm_results";
## -----------------------------------------------------------------------------------

print "First job:\n";
print "$jobarray[0]";

printf "Ready to submit %d jobs to SLURM (y/n): ", scalar @jobarray;
my $response = <STDIN>;
chomp $response;
if($response eq 'n'){
    print "Exiting without any SLURM submissions.\n";
    exit;
}
else{
    print "Proceeding with SLURM submission\n";
}
unless(-e $logdir){
    mkdir  $logdir or die "Failed to make log file directory";
}
unless(-e $resultdir ){
    mkdir  $resultdir or die "Failed to make result directory";
}

###### ------------------ BUILD SLURM SCRIPT---------------
my @slurmdirectives = "#!/bin/bash";

push @slurmdirectives, "#SBATCH --account=$arccproject";
push @slurmdirectives, "#SBATCH --job-name=$jobname";
push @slurmdirectives, "#SBATCH --time=$runtime"; 
push @slurmdirectives, "#SBATCH --nodes=1";
push @slurmdirectives, "#SBATCH --ntasks-per-node=16"; # one core per node
push @slurmdirectives, "#SBATCH --chdir=$logdir";
##push @slurmdirectives, "#SBATCH --partition=teton-hugemem";
#          SLURM can send informative email messages to you about the
#          status of your job.  
if($sendmail eq 'TRUE'){
    push @slurmdirectives, '#SBATCH --mail-type=END';
    push @slurmdirectives, "#SBATCH --mail-user $emailaddress";
}
my $pbsconf = join "\n", @slurmdirectives;

#  By default standard output and error streams are to be merged,
#  intermixed, as standard output.

##########################################
# 
#   Output some useful job information.  #
#
##########################################

my @slurmjob = ();
push @slurmjob, '
    echo ------------------------------------------------------
    echo wrapSLURM: Job is running on node $HOSTNAME
    echo ------------------------------------------------------
    echo wrapSLURM: submitting host was $SLURM_SUBMIT_HOST
    echo wrapSLURM: job identifier is $SLURM_JOB_ID
    echo wrapSLURM: job name is $SLURM_JOB_NAME
    echo ------------------------------------------------------
';

###############################################################                                                    
#   The prolog script automatically makes a directory on the local
#   disks for you.  The name of this directory depends on the job id,
#   but you need only refer to it using ${WORKDIR}.
##############################################################

push @slurmjob, "WORKDIR=$workdir/SLURM_\$SLURM_JOB_ID";
push @slurmjob, "SCP='/usr/bin/scp -o StrictHostKeyChecking=no'";

######################################################################
#   To minimize communications traffic, it is best for your job to
#   work with files on the local disk of the compute node.  Hence, one
#   needs to transfer files from your permanent home directory tree to
#   the directory ${WORKDIR} automatically created by wrapSLURM on the
#   local disk before program execution, and to transfer any important
#   output files from the local disk back to the permanent home
#   directory tree after program execution is completed.  We use
#   secure copy (scp) to do the file transfers to avoid distributed
#   filesystem bottlenecks.
######################################################################

#####################################################
#    Specify the permanent directory(ies) on the server host.  Note
#    that when the job begins execution, the current working directory
#    at the time the qsub command was issued becomes the current
#    working directory of the job.
#####################################################

push @slurmjob, "BASESTOREDIR=$basestoredir";
push @slurmjob, "RESULTDIR=$resultdir";

push @slurmjob, '
    echo workdir is $WORKDIR
    echo basestoredir is $BASESTOREDIR
    echo resultdir is $RESULTDIR
    echo ------------------------------------------------------
    echo \' \'
';

###############################################################
#                                                             #
#    setup WORKDIR (typically local disk)
#                                                             #
###############################################################

my $stagein = '
stagein()
{
    echo \' \'
	echo Setting up local working directory ${WORKDIR}
	echo Optionally transferring files from server to compute node
    mkdir ${WORKDIR}
    cd ${WORKDIR}
    ';
foreach my $module (@modules){
$stagein .= "$module\n";
}

$stagein .= '
    ###    ${SCP} ${BASESTOREDIR}/input_file .
}
';
push @slurmjob, $stagein;

############################################################
#                                                          #
#    Execute the run.  Do not run in the background.       #
#                                                          #
############################################################

push @slurmjob, "runprogram()";
push @slurmjob, "{\n";
## gather everything up to this point, use join with "\n"
my $prolog = join "\n", @slurmjob;
@slurmjob = ();

## this is where we would put the executable if this were not in a perl wrapper:
## program_executable < input_file > output_file

### --- beginning of epilog
push @slurmjob, "}";

###########################################################
#                                                         
#   Copy necessary files back to permanent directory and remove results on node     
#                                                         
###########################################################

push @slurmjob, '
stageout()
{
 echo \' \'
 echo Transferring files from compute nodes to server
 echo Writing files in permanent directory  ${RESULTDIR}
 
 cd ..  ## to parent directory of WORKDIR, for writing
 tar czf SLURM_${SLURM_JOBID}_results.tgz SLURM_${SLURM_JOBID}/
 cp SLURM_${SLURM_JOBID}_results.tgz  ${RESULTDIR}/
 ## clean up
 if [ -e "${RESULTDIR}/SLURM_${SLURM_JOBID}_results.tgz" ]
 then
    rm -rf ${WORKDIR}
    rm  SLURM_${SLURM_JOBID}_results.tgz
 fi
}
';

#####################################################################
#  A slurm command can be used to kill a running job.  It first sends
#  a SIGTERM signal, then after a delay (specified by the "kill_delay"
#  queue attribute (set to 30 seconds), it sends a SIGKILL signal
#  which eradicates the job.  During the time between the SIGTERM and
#  SIGKILL signals, the "cleanup" function below is run. You should
#  include in this function commands to copy files from the local disk
#  back to your home directory.  Note: if you need to transfer very
#  large files which make take longer than 30 seconds, change the
#  KillWait variable
#####################################################################

push @slurmjob, '
early()
{
 echo \' \'
 echo \' ############ WARNING:  EARLY TERMINATION ############# \'
 echo \' \'
}

trap \'early; stageout\' 2 9 15

';

##################################################                                               
#   Staging in, running the job, and staging out were specified above
#   as functions.  Now call these functions to perform the actual file
#   transfers and program execution.
#  #################################################

push @slurmjob, "stagein\n";
push @slurmjob, "runprogram\n";
push @slurmjob, "stageout\n";
my $epilog = join '', @slurmjob;

###### use loop to submit whole @jobarray ########------------
foreach my $job (0..($#jobarray - 1)){
  runserialjob($job);
}
## final job
runserialjob($#jobarray);


#### -------------------------------------------------------------------
sub runserialjob{
    my $j = $_[0];
    my $slurmjob = '';
    $slurmjob .= $pbsconf;
    $slurmjob .= $prolog;
    $slurmjob .= $jobarray[$j];
    $slurmjob .= $epilog;
    $slurmjob .= "exit\n";
    open SBATCH, "| sbatch 1>/dev/null" or die "Failed to fork for sbatch; $!";
    print SBATCH "$slurmjob";
    close(SBATCH) or die "Couldn't close SBATCH";
}

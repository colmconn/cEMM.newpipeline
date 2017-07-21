#!/bin/bash

## set -x

studyName=cEMM.machlearn

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/$studyName}
DATA=$ROOT/data
LOG_DIR=$ROOT/log
SCRIPTS_DIR=${ROOT}/scripts

. ${SCRIPTS_DIR}/logger_functions.sh

if [[ $# -gt 0 ]] ; then
    subjects="$*"
else

    # subjects="369_A 371_A 373_A 377_A"
    subjects=$( cd ../data; ls -1d *_[ACD] *_[ACD]2 )    
    subjectCount=$( cd ../data; ls -1d *_[ACD] *_[ACD]2 | wc -l )
fi

subjectCount=$( echo $subjects | wc -w )

taskName=alignTest
taskFile=$SCRIPTS_DIR/run/${taskName}-TaskFile.$BASHPID
info_message_ln "List of tasks to be executed is stored in $taskFile"

cat /dev/null > $taskFile


(( i=1 ))
for subject in ${subjects} ; do
    info_message_ln "$( printf "Adding script(s) for subject %s (%03d of %03d) to task file\n" $subject $i $subjectCount )"
    
    echo "$SCRIPTS_DIR/alignment_test.sh $subject" >> ${taskFile}
    
    (( i=i+1 ))
done

## cat ${taskFile}
## exit

## jobname
#$ -N $taskName

## queue
#$ -q all.q

## binary? 
#$ -b y

## rerunnable?
#$ -r y

## merge stdout and stderr?
#$ -j y

## send no mail
#$ -m n

## execute from the current working directory
#$ -cwd

## use a shell to run the command
#$ -shell yes 

## set the shell
#$ -S /bin/bash

## preserve environment
#$ -V 

[[ ! -d $LOG_DIR ]] && mkdir $LOG_DIR

nTasks=$( cat $taskFile | wc -l )
sge_command="qsub -N $taskName -q all.q -j y -m n -V -wd $( pwd ) -o $LOG_DIR -t 1-$nTasks" 
echo $sge_command
( exec $sge_command <<EOF
#!/bin/sh

#$ -S /bin/sh

command=\`sed -n -e "\${SGE_TASK_ID}p" $taskFile\`

exec /bin/sh -c "\$command"
EOF
)

echo "Running qstat"
qstat

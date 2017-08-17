#!/bin/bash

## set -x 

programName=`basename $0`

trap exit SIGHUP SIGINT SIGTERM

studyName=cEMM.newpipeline

GETOPT=$( which getopt )
ROOT=/data/sanDiego/$studyName
DATA=$ROOT/data
## this variable is used to search for a preexisting non-linear warped
## anatomy to avoid having to compute the warp again since it's
## so time consuming to do so
MACHLEARN_DATA=/data/sanDiego/cEMM.machlearn/data
LOG_DIR=$ROOT/log
SCRIPTS_DIR=${ROOT}/scripts

. ${SCRIPTS_DIR}/logger_functions.sh

GETOPT_OPTIONS=$( $GETOPT \
		      -o "fe:m:o:h:b:t:np:q" \
		      --longoptions "force,excessiveMotionThresholdFraction:,motionThreshold:,outlierThreshold:,threads:,blur:,tcat:,nonlinear,polort:,enqueue" \
		      -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

## 1 = force creation of zero padded files
force=0

## enqueue the job for execution
enqueue=0

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-f|--force)
	    force=1; shift 1;;
	-e|--excessiveMotionThresholdFraction)
	    excessiveMotionThresholdFraction=$2; shift 2 ;;	
	-m|--motionThreshold)
	    motionThreshold=$2; shift 2 ;;	
	-o|--outlierThreshold)
	    outlierThreshold=$2; shift 2 ;;	
	-h|--threads)
	    threads=$2; shift 2 ;;	
	-b|--blur)
	    blur=$2; shift 2 ;;	
	-t|--tcat)
	    tcat=$2; shift 2 ;;	
	-n|--nonlinear)
	    nonlinear=1; shift 1 ;;
	-p|--polort)
	    polort="-regress_polort $2"; shift 2;;	
	-q|--enqueue)
	    enqueue=1; shift 1 ;;	
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

# if [[ $force -eq 1 ]] ; then
#     info_message_ln "Forcing recreation of ZEROPADed files"
# fi

####################################################################################################
## Check that appropriate values are used to initialize arguments that
## control analysis if no values were provided on the command line

## The following values are used to exclude subjects based on the
## number of volumes censored during analysis
if [[ "x$excessiveMotionThresholdFraction" == "x" ]] ; then
    excessiveMotionThresholdFraction=0.2
    excessiveMotionThresholdPercentage=20
    warn_message_ln "No excessiveMotionThresholdFraction threshold was provided. Defaulting to $excessiveMotionThresholdFraction => ${excessiveMotionThresholdPercentage}%"
else
    excessiveMotionThresholdPercentage=$( echo "(($excessiveMotionThresholdFraction*100)+0.5)/1" | bc ) 

    info_message_ln "Using ${excessiveMotionThresholdFraction} as the subject exclusion motion cutoff fraction"
    info_message_ln "Using ${excessiveMotionThresholdPercentage}% as subject exclusion motion cutoff percentage"
    info_message_ln "Note that these values are used to exclude subjects based on the number of volumes censored during analysis"
fi


## motionThreshold and outlierThreshold are the values passed to
## afni_proc.py and are used when deciding to censor a volume or not
if [[ "x${motionThreshold}" == "x" ]] ; then
    motionThreshold=0.2
    warn_message_ln "No motionThreshold value was provided. Defaulting to $motionThreshold"
else
    info_message_ln "Using motionThreshold of ${motionThreshold}"
fi

if [[ "x${outlierThreshold}" == "x" ]] ; then
     outlierThreshold=0.1
     warn_message_ln "No outlierThreshold value was provided. Defaulting to $outlierThreshold"
else
    info_message_ln "Using outlierThreshold of ${outlierThreshold}"
fi

if [[ "x${threads}" == "x" ]] ; then
     threads=1
     warn_message_ln "No value for the number of parallel threads to use was provided. Defaulting to $threads"
else
    info_message_ln "Using threads value of ${threads}"
fi

if [[ "x${blur}" == "x" ]] ; then
     blur="6"
     warn_message_ln "No value for blur filter value to use was provided. Defaulting to $blur"
else
    info_message_ln "Using blur filter value of ${blur}"
fi

if [[ "x${tcat}" == "x" ]] ; then
     tcat="3"
     warn_message_ln "No value for tcat, the number of TRs to censor from the start of each volume, was provided. Defaulting to $tcat"
else
    info_message_ln "Using tcat filter value of ${tcat}"
fi

if [[ $nonlinear -eq 1 ]] ; then 
    info_message_ln "Using nonlinear alignment"
    scriptExt="NL"
else 
    info_message_ln "Using affine alignment only"
    scriptExt="aff"    
fi

####################################################################################################
if [[ "$#" -gt 0 ]] ; then
    subjects="$@"
else
    ## subjects=$( cd $DATA ; find ./ -maxdepth 1 \( -name '[0-9][0-9][0-9]_[ACD]' -o -name '[0-9][0-9][0-9]_[ACD]?' \) -printf "%f " )
    ## subjects=$( cd $DATA ; find ./ -maxdepth 1 \( -name '[0-9][0-9][0-9]_[A]' -o -name '[0-9][0-9][0-9]_[A]?' \) -printf "%f " )
    subjects=$( cd $DATA  ; ls -d [0-9][0-9][0-9]_{C,D} [0-9][0-9][0-9]_{C,D}2 2> /dev/null | grep -v 999 )
fi

[[ -d run ]] || mkdir run

for subject in $subjects ; do
    info_message_ln "#################################################################################################"
    info_message_ln "Generating script for subject $subject"

    if  [[ ! -f ${DATA}/$subject/${subject}EMM+orig.HEAD ]] ; then
	warn_message_ln "Can not find EMM EPI file for ${subject}. Skipping."
	continue
    else
	epiFile=${DATA}/$subject/${subject}EMM+orig.HEAD
    fi

    if  [[ ! -f ${DATA}/$subject/${subject}+orig.HEAD ]] ; then
	warn_message_ln "Can not find anatomy file for subject ${subject}. Skipping."
	continue
    else
	if [[ -f ${DATA}/$subject/${subject}_clp+orig.HEAD ]] ; then
	    anatFile=${DATA}/$subject/${subject}_clp+orig.HEAD
	else
	    anatFile=${DATA}/$subject/${subject}+orig.HEAD
	fi
    fi

    if [[ $nonlinear -eq 1 ]] ; then 
	outputScriptName=run/run-afniEmmPreproc-${subject}.${scriptExt}.sh
    else
	outputScriptName=run/run-afniEmmPreproc-${subject}.${scriptExt}.sh	
    fi

    if [[ -f EMM_alignment_parameters.sh ]] ; then
	info_message_ln "Loading alignment parameters from EMM_alignment_parameters.sh"
	source EMM_alignment_parameters.sh
    else
	extraAlignmentArgs="-align_opts_aea -partial_axial"	
	info_message_ln "Setting extra alignment options to default of $extraAlignmentArgs"
    fi  
    
    ## do non-linear warping? If so add the flag to the extra
    ## alignment args variable
    if [[ $nonlinear -eq 1 ]] ; then 
	extraAlignmentArgs="${extraAlignmentArgs} -tlrc_NL_warp"

	anat_base=$( basename $anatFile )
	anat_base=${anat_base%%+*}
	info_message_ln "Looking for pre-existing non-linear warped T1w images in ${MACHLEARN_DATA}/${subject}/afniEmmPreprocessed.NL"
	if [[ -f ${MACHLEARN_DATA}/${subject}/afniEmmPreprocessed.NL/${anat_base}_al_keep+tlrc.HEAD ]] && \
	   [[ -f ${MACHLEARN_DATA}/${subject}/afniEmmPreprocessed.NL/anat.un.aff.Xat.1D ]] && \
	   [[ -f ${MACHLEARN_DATA}/${subject}/afniEmmPreprocessed.NL/anat.un.aff.qw_WARP.nii.gz ]] ; then
	    info_message_ln "Yay! Found them!!!"
	    info_message_ln "Supplying pre-existing nonlinear warped anatomy to afni_proc.py"
	    extraAlignmentArgs="${extraAlignmentArgs} \\
	     -tlrc_NL_warped_dsets ${MACHLEARN_DATA}/${subject}/afniEmmPreprocessed.NL/${anat_base}_al_keep+tlrc.HEAD \\
                                   ${MACHLEARN_DATA}/${subject}/afniEmmPreprocessed.NL/anat.un.aff.Xat.1D \\
                                   ${MACHLEARN_DATA}/${subject}/afniEmmPreprocessed.NL/anat.un.aff.qw_WARP.nii.gz"
	fi
    fi


    cd ${DATA}/${subject} 
    regressorFile=$DATA/regressors/${subject}_EMM.wav.1D
    if  [[ ! -f ${regressorFile} ]] ; then
    	warn_message_ln "Can not find regressor file (${regressorFile}) for ${subject}. Skipping."
    	continue
    else
    	info_message_ln "Splitting regressors into individual files"
	
    	regressorLabels=("Oval" "Fear" "Happy" "Sad")
    	nRegressorLabels=${#regressorLabels[@]}
    	nRegressorColumns=$( head -1 ${regressorFile} | wc -w )

    	if [[ $nRegressorLabels -gt $nRegressorColumns ]] ; then
    	    warn_message_ln "Number of regressor labels ($nRegressorLabels) is greater than the number of regressor columns ($nRegressorColumns) for subject $subject. Skipping subject"
    	    continue
    	fi
    	splitRegressorFiles=""
    	for (( i=0; i<${nRegressorLabels}; i++ )) ; do
    	    info_message_ln "Mapping column: $i -> ${regressorLabels[$i]}"
    	    splitRegressorFile="${regressorFile%%.wav.1D}__${regressorLabels[$i]}.wav.1D"
    	    ## chop off any leading directory component so that the
    	    ## split files are put in the current working directory
    	    splitRegressorFile=${splitRegressorFile##*/}
    	    (( cc=i+1 ))
    	    awk -v col=$cc '{ print $col }' < $regressorFile > $splitRegressorFile
    	    splitRegressorFiles="${splitRegressorFiles} $splitRegressorFile"
    	done
    fi

    info_message_ln "Making regressor timeseries graph"
    1dplot -plabel $( echo ${subject} | sed -e "s/_/-/" ) \
    	   -sepscl \
    	   -png ${subject}.EMM.timeseries.png \
    	   -ynames $( echo ${regressorLabels[@]} |  sed -e "s/_/-/g" ) \
    	   - ${splitRegressorFiles}
    
    cd $SCRIPTS_DIR
    
    info_message_ln "Writing script: $outputScriptName"
    
    cat <<EOF > $outputScriptName
#!/bin/bash

set -x 

#$ -S /bin/bash

## disable compression of BRIKs/nii files
unset AFNI_COMPRESSOR

export PYTHONPATH=$AFNI_R_DIR

## use the newer faster despiking method. comment this out to get the
## old one back
export AFNI_3dDespike_NEW=YES

# turn off anoying colorization of info/warn/error messages since they
# only result in gobbledygook
export AFNI_MESSAGE_COLORIZE=NO

## only use a single thread since we're going to run so many subjects
## in parallel
export OMP_NUM_THREADS=${threads}

excessiveMotionThresholdFraction=$excessiveMotionThresholdFraction
excessiveMotionThresholdPercentage=$excessiveMotionThresholdPercentage

cd $DATA/$subject

preprocessingScript=${subject}.afniEmmPreprocess.$scriptExt.csh
rm -f \${preprocessingScript}

outputDir=afniEmmPreprocessed.$scriptExt
rm -fr \${outputDir}

motionThreshold=${motionThreshold}
outlierThreshold=${outlierThreshold}

##	     -tcat_remove_first_trs ${tcat}					\\
## -tlrc_opts_at -init_xform AUTO_CENTER \\
## 	     -regress_censor_outliers \$outlierThreshold                 	\\

afni_proc.py -subj_id ${subject}						\\
             -script \${preprocessingScript}					\\
	     -out_dir \${outputDir}						\\
	     -blocks tshift align tlrc volreg mask blur scale regress		\\
	     -copy_anat $anatFile						\\
	     -dsets $epiFile							\\
	     -tlrc_base MNI_caez_N27+tlrc					\\
	     -volreg_align_to first						\\
	     -volreg_tlrc_warp ${extraAlignmentArgs}				\\
	     -blur_size ${blur}							\\
	     -blur_to_fwhm							\\
	     -blur_opts_B2FW "-ACF -rate 0.2"					\\
	     -mask_apply group							\\
	     -regress_censor_first_trs ${tcat}					\\
	     -regress_apply_mot_types demean					\\
             -regress_censor_motion \$motionThreshold ${polort}			\\
	     -regress_censor_outliers \$outlierThreshold			\\
	     -regress_3dD_stop							\\
	     -regress_reml_exec							\\
	     -regress_use_stim_files						\\
	     -regress_stim_files  ${splitRegressorFiles}			\\
	     -regress_stim_labels ${regressorLabels[@]}				\\
	     -regress_run_clustsim no						\\
	     -regress_est_blur_epits						\\
	     -regress_est_blur_errts						\\
	     -regress_opts_3dD							\\
	     	 -num_glt 13							\\
   	     	 -gltsym 'SYM: +1*Fear -1*Oval'					\\
		 -glt_label 1 "Fear-Oval"					\\
		 -gltsym 'SYM: +1*Happy -1*Oval'				\\
		 -glt_label 2 "Happy-Oval"					\\
		 -gltsym 'SYM: +1*Sad -1*Oval'					\\
		 -glt_label 3 "Sad-Oval"					\\
		 -gltsym 'SYM: +0.5*Sad +0.5*Fear -1*Oval'			\\
		 -glt_label 4 "Negative-Oval"					\\
		 -gltsym 'SYM: +1*Fear -1*Happy'				\\
		 -glt_label 5 "Fear-Happy"					\\
		 -gltsym 'SYM: +1*Sad -1*Happy'					\\
		 -glt_label 6 "Sad-Happy"					\\
		 -gltsym 'SYM: +0.5*Sad +0.5*Fear -1*Happy'			\\
		 -glt_label 7 "Negative-Happy"					\\
		 -gltsym 'SYM: +1*Sad -1*Fear'					\\
		 -glt_label 8 "Sad-Fear"					\\
		 -gltsym 'SYM: +0.5*Fear +0.5*Happy -1*Sad'			\\
		 -glt_label 9 "FearHappy-Sad"					\\
		 -gltsym 'SYM: +0.5*Sad +0.5*Happy -1*Fear'			\\
		 -glt_label 10 "SadHappy-Fear"					\\
		 -gltsym 'SYM: +0.5*Fear +0.5*Happy -1*Oval'			\\
		 -glt_label 11 "FearHappy-Oval"					\\
		 -gltsym 'SYM: +0.5*Sad +0.5*Happy -1*Oval'			\\
		 -glt_label 12 "SadHappy-Oval"					\\
		 -gltsym 'SYM: +0.33*Sad +0.33*Fear +0.34*Happy -1*Oval'	\\
		 -glt_label 13 "Faces-Oval"

if [[ -f \${preprocessingScript} ]] ; then 
   tcsh -xef \${preprocessingScript}

    cd \${outputDir}
    xmat_regress=X.xmat.1D 

    if [[ -f \$xmat_regress ]] ; then 

        fractionOfCensoredVolumes=\$( 1d_tool.py -infile \$xmat_regress -show_tr_run_counts frac_cen )
        numberOfCensoredVolumes=\$( 1d_tool.py -infile \$xmat_regress -show_tr_run_counts trs_cen )
        totalNumberOfVolumes=\$( 1d_tool.py -infile \$xmat_regress -show_tr_run_counts trs_no_cen )

        ## rounding method from http://www.alecjacobson.com/weblog/?p=256
        cutoff=\$( echo "((\$excessiveMotionThresholdFraction*\$totalNumberOfVolumes)+0.5)/1" | bc )
	if [[ \$numberOfCensoredVolumes -gt \$cutoff ]] ; then 

	    echo "*** A total of \$numberOfCensoredVolumes of
	    \$totalNumberOfVolumes volumes were censored which is
	    greater than \$excessiveMotionThresholdFraction
	    (n=\$cutoff) of all total volumes of this subject" > \\
		00_DO_NOT_ANALYSE_${subject}_\${excessiveMotionThresholdPercentage}percent.txt

	    echo "*** WARNING: $subject will not be analysed due to having more than \${excessiveMotionThresholdPercentage}% of their volumes censored."
	fi
	
	# make an image to check alignment
	$SCRIPTS_DIR/snapshot_volreg.sh anat_final.${subject}+tlrc pb03.${subject}.r01.volreg+tlrc ${subject}.alignment
    else
	touch 00_DO_NOT_ANALYSE_${subject}_\${excessiveMotionThresholdPercentage}percent.txt
    fi
    echo "Compressing BRIKs and nii files"
    find ./ \( -name "*.BRIK" -o -name "*.nii" \) -print0 | xargs -0 gzip
else
    echo "*** No such file \${preprocessingScript}"
    echo "*** Cannot continue"
    exit 1
fi	

EOF

    chmod +x $outputScriptName
    if [[ $enqueue -eq 1 ]] ; then
	info_message_ln "Submitting job for execution to queuing system"
	LOG_FILE=$DATA/$subject/$subject-emm-afniPreproc.${scriptExt}.log
	info_message_ln "To see progress run: tail -f $LOG_FILE"
	rm -f ${LOG_FILE}
	qsub -N emm-$subject -q all.q -j y -m n -V -wd $( pwd )  -o ${LOG_FILE} $outputScriptName
    else
	info_message_ln "Job *NOT* submitted for execution to queuing system"
	info_message_ln "Pass -q or --enqueue options to this script to do so"	
    fi

done

if [[ $enqueue -eq 1 ]] ; then 
    qstat
fi

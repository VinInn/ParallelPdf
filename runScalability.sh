#!/bin/sh

if [ "$1" != "nhm8" -a "$1" != "wes12" -a "$1" != "amd48" -a "$1" != "nex64" -a "$1" != "snb8" -a "$1" != "snb32" -a "$1" != "snb64" -a "$1" != "wex80" -a "$1" != "mic" ]; then
    echo "ERROR: Unknown switch '$1'. Accepted values: nhm8, wes12, amd48, nex64, snb8, snb32, wex80, mic"
    exit 1;
fi

ARCH=$1
APPLICATION=./main_XEON
#Use MINOS
#OEVENTS="-m -e -n"
OEVENTS="-e -n"
#OEVENTS="-n"
BLOCKS=-b
TIMES=3
#if MPIRUN is uncommented, then use MPI
MPIRUN="mpirun -np"

#FILENAME=scalability_MIC_NoVec_NoFMA_BlockSplttingExternalLoop_oplaboum2_mic0
#FILENAME=scalability_MIC_Vec_BlockSplttingExternalLoopTest_MINOS_oplaboum3_mic0
#FILENAME=scalability_MIC_Vec_NoFMA_BlockSplttingExternalLoopTestSMT_3Core_oplaboum2_mic0_polblocks
#FILENAME=scalability_MIC_Vec_NoFMA_BlockSplttingExternalLoop_oplaboum2_mic0_polblocks_FPreduce_EDM2
#FILENAME=scalability_IVY_AVX_BlockSplttingExternalLoop_opladev37
#FILENAME=scalability_MIC_Vec_WithFMA_BlockSplttingExternalLoop_oplaboum2_mic0
#FILENAME=scalability_MIC_Vec_BlockSplttingExternalLoopTestSMT_1Core_oplaboum3_mic0
#FILENAME=scalability_MIC_Vec_BlockSplttingExternalLoop_oplaboum2_mic0
#FILENAME=scalability_SNB_SSE3_BlockSplttingExternalLoopNoTurbo_opladev36
#FILENAME=scalability_SNB_AVX_BlockSplttingExternalLoopNoTurboWithMPI_TAU_mcs
#FILENAME=scalability_SNB_AVX_BlockSplttingExternalLoop_opladev36
#FILENAME=scalability_SNB64_AVX_BlockSplttingExternalLoop_devrac5
FILENAME=scalability_WES_SSE_BlockSplttingExternalLoop_olslc6x

#TAUDIR=profiles_olwork09_olwork10_2MPIperHOST

LOGDIR=logs
LOG=${LOGDIR}/${FILENAME}

FLAG="OpenMP Real Time"
FLAGNC="OpenMP # FCN Calls"
FLAGFCN="OpenMP Result"

#EVENTS="500000"
EVENTS="1000000"

BLOCKSIZES="1000 2000 3000 5000 10000"


if [ "$ARCH" = "nhm8" ]; then
# Across sockets
    export KMP_AFFINITY="verbose,granularity=fine,proclist=[1,4,2,5,3,6,7,0,9,12,10,13,11,14,8,15],explicit"
# Inside sockets (not good for caches) 
#    export KMP_AFFINITY="verbose,granularity=fine,proclist=[4-7,1-3,0,12-15,9-11,8],explicit"
    THREADS="1 2 4 6 8 10 12 16"

elif [ "$ARCH" = "wes12" ]; then
# opladev31 & olslc6x topology 
# Across sockets
    export KMP_AFFINITY="verbose,granularity=fine,proclist=[1,6,2,7,3,8,4,9,5,10,0,11,13,18,14,19,15,20,16,21,17,22,12,23],explicit"

# Inside sockets (not good for caches) 
#    export KMP_AFFINITY="verbose,granularity=fine,proclist=[1-11,0,13-23,12],explicit"
    
# opladev32 topology
# Across sockets
#    export KMP_AFFINITY="verbose,granularity=fine,proclist=[2,1,4,3,6,5,8,7,10,9,0,11,14,13,16,15,18,17,20,19,22,21,12,23],explicit"

#    MPIPROCS="1 2"
#    THREADS="1 2 4 6 8 10 12 24"
#    NODESLIST="-hosts olwork09,olwork10"
##    NODESLIST="-hosts olwork09,olwork10,olwork12,olwork13,olwork14"
#    MPIPROCS="1 2 4"
    MPIPROCS="1 2"
#    THREADS="1 2 4 6 8 10 12"
    THREADS="1 12 24"
#    MPIPERHOST="2"
#    NODESLIST="-hosts olwork09,olwork10"

elif [ "$ARCH" = "amd48" ]; then
    TASKSET="taskset -c 1,6,12,18,24,30,36,42,2,7,13,19,25,31,37,43,3,8,14,20,26,32,38,44,4,9,15,21,27,33,39,45,5,10,16,22,28,34,40,46,11,17,23,29,35,41,47,0 "
#    export KMP_AFFINITY="verbose,granularity=fine,proclist=[1-47,0],explicit"
    THREADS="1 2 4 6 8 10 12 16 20 24 32 40 48"

elif [ "$ARCH" = "nex64" ]; then
    export KMP_AFFINITY="verbose,granularity=fine,proclist=[4-31,0-3,36-63,32-35],explicit"
    THREADS="1 2 4 6 8 10 12 16 20 24 32 40 48 56 64"

elif [ "$ARCH" = "snb8" ]; then
    export KMP_AFFINITY="verbose,granularity=fine,proclist=[1-3,0,5-7,4],explicit"
    BLOCKSIZES="1000 5000 10000"
#    MPIPROCS="1 2 3 4 5"
    THREADS="1 2 4 8"
#    NODESLIST="-hosts oplamcs01,oplamcs02,oplamcs03,oplamcs04,oplamcs05"

elif [ "$ARCH" = "snb32" ]; then
    export KMP_AFFINITY="verbose,granularity=fine,proclist=[1,8,2,9,3,10,4,11,5,12,6,13,7,14,0,15,17,24,18,25,19,26,20,27,21,28,22,29,23,30,16,31],explicit"
#    MPIPROCS="1 2 4 8 16 32"
    MPIPROCS="1 2 4"
#    THREADS="1 2 4 8 16 32"
    THREADS="1 2"

elif [ "$ARCH" = "snb64" ]; then
    export KMP_AFFINITY="verbose,granularity=fine,proclist=[4-31,0-3,36-63,32-35],explicit"
    MPIPROCS="1 2 4 8"
    THREADS="1 2 3 4 5 6 8 10 12 16 20 24 32 64"

elif [ "$ARCH" = "wex80" ]; then
    export KMP_AFFINITY="verbose,granularity=fine,proclist=[16,2,1,3,64,18,17,19,8,66,65,67,24,10,9,11,4,26,25,27,20,6,5,7,68,22,21,23,12,70,69,71,28,14,13,15,0,30,29,31,48,34,33,35,72,50,49,51,40,74,73,75,56,42,41,43,36,58,57,59,52,38,37,39,76,54,53,55,44,78,77,79,60,46,45,47,32,62,61,63],explicit"
    MPIPROCS="1 2 4 8"
    THREADS="1 2 3 4 5 6 8 10 12 16 20 24 32 40 80"

elif [ "$ARCH" = "mic" ]; then

    APPLICATION=./main_MIC
    source ${PWD}/sourceme_affinity_mic.sh

#    SMT studies
#    KMP_AFFINITY="verbose,granularity=fine,proclist=[5,6,7,8],explicit"
#    export KMP_AFFINITY

#    THREADS="1 2 3 4"
#    THREADS="1 2 4 6 8 10 12 16 24 31 32 63 64 95 96 127 128"
    THREADS="1 2 4 6 8 10 12 16 24 30 45 59 60 119 120 179 180 239 240"
#    THREADS="3"
#    THREADS="32 64 96 128"
#    THREADS="60 120 180 240"

#    BLOCKSIZES="100 500 1000 2000 2700 4000 8000 10000"
#    BLOCKSIZES="0 20000"
    BLOCKSIZES="0"

fi

check_threads_blocks() {
#    $1 is the NBLOCK
#    $2 is the NTHREADS
    if [ "$ARCH" = "mic" ]; then
	return 1
    fi

#:<<COMMENT
    if [[ \
	("$ARCH" = "snb8" && \
        (("$1" = "1000" && (("$2" = "4") || ("$2" = "8"))) || \
        ("$1" = "5000" && "$2" = "2") || \
        ("$1" = "10000" && "$2" = "1"))) || \
	("$ARCH" = "wes12" && \
        (("$1" = "1000" && (("$2" = "12") || ("$2" = "24"))) || \
        ("$1" = "2000" && ("$2" -ge "8" && "$2" -le "10")) || \
        ("$1" = "3000" && "$2" = "6") || \
        ("$1" = "5000" && "$2" = "4") || \
        ("$1" = "10000" && "$2" -le "2"))) || \
	("$ARCH" = "snb32" && \
        (("$1" = "1000" && (("$2" = "16") || ("$2" = "32"))) || \
        ("$1" = "2000" && ("$2" -ge "8" && "$2" -le "12")) || \
        ("$1" = "3000" && "$2" = "6") || \
        ("$1" = "5000" && "$2" = "4") || \
        ("$1" = "10000" && "$2" -le "2"))) || \
	("$ARCH" = "snb64" && \
        (("$1" = "1000" && (("$2" = "32") || ("$2" = "64"))) || \
        ("$1" = "2000" && ("$2" -ge "16" && "$2" -le "24")) || \
        ("$1" = "3000" && "$2" = "12") || \
        ("$1" = "5000" && (("$2" = "6") || ("$2" = "8") || ("$2" = "10"))) || \
        ("$1" = "10000" && (("$2" = "1") || ("$2" = "2") || ("$2" = "4"))))) || \
	("$ARCH" = "wex80" && \
        (("$1" = "1000" && (("$2" = "32") || ("$2" = "40") || ("$2" = "80"))) || \
        ("$1" = "2000" && ("$2" -ge "16" && "$2" -le "24")) || \
        ("$1" = "3000" && "$2" = "12") || \
        ("$1" = "5000" && (("$2" = "6") || ("$2" = "8") || ("$2" = "10"))) || \
        ("$1" = "10000" && (("$2" = "1") || ("$2" = "2") || ("$2" = "4"))))) \
        ]]; then
        return 1
    fi
#COMMENT

    return 0
}

echo "Affinity = "${KMP_AFFINITY}
if [ -n "${MPIRUN}" ]; then
    # use MPI
    # if MPIPROCS is not set, then use $THREADS
    MPIPROCS=${MPIPROCS:-$THREADS}
    echo "# Processes = "${MPIPROCS}

    if [ -n "${MPIPERHOST}" ]; then
	NODESLIST="${NODESLIST} -perhost ${MPIPERHOST}"
    fi

else
    # do not use MPI
    MPIPROCS="1"
fi
echo "# Threads = "${THREADS}

mkdir -p ${LOGDIR}

rm -f ${FILENAME}*_${ARCH}.txt
rm -f ${LOG}*.log

for NEVENTS in ${EVENTS}
do
    for NBLOCK in ${BLOCKSIZES}
    do 
	
	if [ "$NEVENTS" -gt "$NBLOCK" ]; then 

	    for NMPIPROCS in ${MPIPROCS}
	    do

		if [ -n "${MPIRUN}" ]; then
		    MPICMD="${MPIRUN} ${NMPIPROCS}"
		fi

		for NTHREADS in ${THREADS}
		do

		    if [ -n "${NODESLIST}" ]; then

			if [ -n "${MPIPERHOST}" ]; then
			    TOTNTHREADS=`expr \( \( $NMPIPROCS \- 1 \) \% $MPIPERHOST \+ 1 \) \* $NTHREADS`
			else
			    TOTNTHREADS="$NTHREADS"
			fi

			check_threads_blocks $NBLOCK $TOTNTHREADS
			if [ "$?" -eq "0" ]; then
			    continue
			fi

		    else

			TOTNTHREADS=`expr $NMPIPROCS \* $NTHREADS`

			check_threads_blocks $NBLOCK $TOTNTHREADS
			if [ "$?" -eq "0" ]; then
			    continue
			fi

			if [ -n "${NTHREADSMAX}" ]; then
			    if [ "$TOTNTHREADS" -gt "$NTHREADSMAX" ]; then
				break
			    fi
			fi
		    fi
		
		    SUMCPU=0.0
		    SUMCPU2=0.0

		    export OMP_NUM_THREADS=$NTHREADS

		    for NTIMES in $(seq 1 ${TIMES})
		    do
			echo "**** # of events = ${NEVENTS} / Blocksize = ${NBLOCK} / # Processes = ${NMPIPROCS} / # Threads = ${NTHREADS} / Execution ${NTIMES} ****"
			LOGNAME=${LOG}_${NEVENTS}_${NMPIPROCS}_${NBLOCK}_${NTHREADS}_${NTIMES}.log
			echo ${MPICMD} ${NODESLIST} ${TASKSET} ${APPLICATION} ${OEVENTS} ${NEVENTS} ${BLOCKS} ${NBLOCK}
			
			# TAU FILES
			if [ -n "${TAUDIR}" ]; then
			    export PROFILEDIR="${TAUDIR}/NEVENTS${NEVENTS}_BLOCK${NBLOCK}_MPI${NMPIPROCS}_OMP${NTHREADS}_${NTIMES}"
			    mkdir -p ${PROFILEDIR}
			fi


			${MPICMD} ${NODESLIST} ${TASKSET} ${APPLICATION} ${OEVENTS} ${NEVENTS} ${BLOCKS} ${NBLOCK} > ${LOGNAME} 2>&1
			NCALLS=`grep "${FLAGNC}" ${LOGNAME} | awk 'BEGIN { FS = "=" } ; { print $2 }'`
			TIMECPU=`grep "${FLAG}" ${LOGNAME} | awk 'BEGIN { FS = "=" } ; { print $2 }'`
			MINFCN=`grep "${FLAGFCN}" ${LOGNAME} | awk 'BEGIN { FS = "=" } ; { print $2 }'`
		    
			SUMCPU=`echo "sum" | awk -v var1=${SUMCPU} -v var2="${TIMECPU}" '{ print var1 + var2 }'`
			SUMCPU2=`echo "sum2" | awk -v var1=${SUMCPU2} -v var2="${TIMECPU}" '{ print var1 + var2 * var2 }'`
		    
		    done
		
		    LINE=`echo "div" | awk -v var1=$SUMCPU -v var2=$TIMES '{ print var1 / var2 }'`
		    LINE=`echo $LINE ; echo "div" | awk -v var1=$SUMCPU2 -v var2=$TIMES -v var3=$SUMCPU '{ print sqrt((var1/var2)-(var3*var3)/(var2*var2)) }'`
		
		    echo ${NMPIPROCS} $NTHREADS $LINE $NCALLS ${NBLOCK} ${NEVENTS} $MINFCN >> ${FILENAME}_${NEVENTS}_${NBLOCK}_${ARCH}.txt

		done

		if [ -z ${NTHREADSMAX} ]; then
		    NTHREADSMAX=${NTHREADS}
                fi

	    done
	fi
    done

done

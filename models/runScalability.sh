#!/bin/sh

if [ "$1" != "nhm8" -a "$1" != "wes12" -a "$1" != "amd48" -a "$1" != "nex32" -a "$1" != "snb4" -a "$1" != "snb16" -a "$1" != "wex40" -a "$1" != "mic" ]; then
    echo "ERROR: Unknown switch '$1'. Accepted values: nhm8, wes12, amd48, nex64, snb4, snb32, wex80, mic"
    exit 1;
fi

APPLICATION=./modelRooFit
OEVENTS=-n
OBLOCKS=-b
OTHREADS=-t
TIMES=3
FILENAME=scalability
LOGDIR=logs
LOG=${LOGDIR}/${FILENAME}
FLAG="Real Time"
FLAGNC="# FCN calls"
FLAGFCN="Result"
EVENTS="50000"
#EVENTS="10000 50000 100000 250000 500000 1000000"
#EVENTS="500000"

BLOCKSIZES="0"

if [ "$1" = "nhm8" ]; then
# Across sockets
    export KMP_AFFINITY="verbose,granularity=fine,proclist=[1,4,2,5,3,6,7,0,9,12,10,13,11,14,8,15],explicit"
# Inside sockets (not good for caches) 
#    export KMP_AFFINITY="verbose,granularity=fine,proclist=[4-7,1-3,0,12-15,9-11,8],explicit"
    THREADS="1 2 4 6 8 10 12 16"
#    BLOCKSIZES="0 30000 20000 15000 10000 5000 2000"

elif [ "$1" = "wes12" ]; then
# opladev31 topology
# Across sockets
#    export KMP_AFFINITY="verbose,granularity=fine,proclist=[1,6,2,7,3,8,4,9,5,10,0,11,13,18,14,19,15,20,16,21,17,22,12,23],explicit"
# Inside sockets (not good for caches) 
#    export KMP_AFFINITY="verbose,granularity=fine,proclist=[6-11,1-5,0,12-23],explicit"
    
# opladev32 topology
# Across sockets
    export KMP_AFFINITY="verbose,granularity=fine,proclist=[2,1,4,3,6,5,8,7,10,9,0,11,14,13,16,15,18,17,20,19,22,21,12,23],explicit"
    THREADS="1 2 4 6 8 10 12 16 20 24"

elif [ "$1" = "amd48" ]; then
    TASKSET="taskset -c 1,6,12,18,24,30,36,42,2,7,13,19,25,31,37,43,3,8,14,20,26,32,38,44,4,9,15,21,27,33,39,45,5,10,16,22,28,34,40,46,11,17,23,29,35,41,47,0 "
#    export KMP_AFFINITY="verbose,granularity=fine,proclist=[1-47,0],explicit"
    THREADS="1 2 4 6 8 10 12 16 20 24 32 40 48"

elif [ "$1" = "nex64" ]; then
    export KMP_AFFINITY="verbose,granularity=fine,proclist=[4-31,0-3,36-63,32-35],explicit"
    THREADS="1 2 4 6 8 10 12 16 20 24 32 40 48 56 64"

elif [ "$1" = "snb4" ]; then
    export KMP_AFFINITY="verbose,granularity=fine,proclist=[1-3,0,5-7,4],explicit"
    THREADS="1 2 4 6 8"

elif [ "$1" = "snb16" ]; then
    export KMP_AFFINITY="verbose,granularity=fine,proclist=[1,8,2,9,3,10,4,11,5,12,6,13,7,14,0,15],explicit"
    THREADS="1 2 4 6 8 10 12 16 20 24 32"

elif [ "$1" = "wex80" ]; then
    export KMP_AFFINITY="verbose,granularity=fine,proclist=[16,2,1,3,64,18,17,19,8,66,65,67,24,10,9,11,4,26,25,27,20,6,5,7,68,22,21,23,12,70,69,71,28,14,13,15,0,30,29,31,48,34,33,35,72,50,49,51,40,74,73,75,56,42,41,43,36,58,57,59,52,38,37,39,76,54,53,55,44,78,77,79,60,46,45,47,32,62,61,63],explicit"
    THREADS="1 2 4 6 8 10 12 16 20 24 32 40 48 56 64 72 80"

elif [ "$1" = "mic" ]; then

    APPLICATION=${APPLICATION}_$1
    THREADS="1 2 4 6 8 10 12 16 24 30 60 90 120"
    TIMES=1

fi

echo "Affinity = "${KMP_AFFINITY}
echo "# Threads = "${THREADS}

mkdir -p ${LOGDIR}

log_file_name()
{
    echo ${LOG}_${NEVENTS}_${NTHREADS}_${NTIMES}_$1$2.log
}

output_file_name()
{
    echo ${FILENAME}_${NEVENTS}_${NBLOCK}_$1$2.txt
}

echo "Delete old logs and ouput files"
for NEVENTS in ${EVENTS}
  do
  for NBLOCK in ${BLOCKSIZES}
    do
    if [ "$NEVENTS" -gt "$NBLOCK" ]; then 
   
	for NTHREADS in ${THREADS}
	  do
	  for NTIMES in $(seq 1 ${TIMES})
	    do
	    rm -f $(log_file_name $1 $2)
	  done
	done
	rm -f $(output_file_name $1 $2)
    fi
  done
done

for NEVENTS in ${EVENTS}
do
  for NBLOCK in ${BLOCKSIZES}
    do 
     
    if [ "$NEVENTS" -gt "$NBLOCK" ]; then 
   
    for NTHREADS in ${THREADS}
    do
      
	SUMCPU=0.0
	SUMCPU2=0.0

	export OMP_NUM_THREADS=${NTHREADS}

	for NTIMES in $(seq 1 ${TIMES})
	do
	    echo "**** # of events = ${NEVENTS} / # Blocksize = ${NBLOCK} / # Threads = ${NTHREADS} / Execution ${NTIMES} ****"
	    LOGNAME=$(log_file_name $1 $2)
	    ${TASKSET} ${APPLICATION} ${OEVENTS} ${NEVENTS} ${OBLOCKS} ${NBLOCKS} ${OTHREADS} ${NTHREADS} > ${LOGNAME} 2>&1
	    NCALLS=`grep "${FLAGNC}" ${LOGNAME} | awk 'BEGIN { FS = "=" } ; { print $2 }'`
	    TIMECPU=`grep "${FLAG}" ${LOGNAME} | awk 'BEGIN { FS = "=" } ; { print $2 }'`
	    MINFCN=`grep "${FLAGFCN}" ${LOGNAME} | awk 'BEGIN { FS = "=" } ; { print $2 }'`
	    
	    SUMCPU=`echo "sum" | awk -v var1=${SUMCPU} -v var2="${TIMECPU}" '{ print var1 + var2 }'`
	    SUMCPU2=`echo "sum2" | awk -v var1=${SUMCPU2} -v var2="${TIMECPU}" '{ print var1 + var2 * var2 }'`
	    
	done

	LINE=`echo "div" | awk -v var1=${SUMCPU} -v var2=${TIMES} '{ print var1 / var2 }'`
	LINE=`echo $LINE ; echo "div" | awk -v var1=${SUMCPU2} -v var2=${TIMES} -v var3=${SUMCPU} '{ print sqrt((var1/var2)-(var3*var3)/(var2*var2)) }'`

	echo ${NTHREADS} $LINE ${NCALLS} ${NBLOCK} ${NEVENTS} ${MINFCN} >> $(output_file_name $1 $2)

    done
   fi
done

done


# Set the affinity on MIC

KMP_AFFINITY="verbose,granularity=fine,proclist=["

# Maximum number of threads per core
THREADS=4

# Note that CORES starts from 0 and they are minus 1, i.e. for 32 cores you should have CORES=30
CORES=`expr \( \` tail -26 /proc/cpuinfo | grep processor | cut -f 2 -d : \` + 1 \) / 4 - 2`
echo "# Cores = " `expr ${CORES} + 2`

for THREAD in $(seq 1 ${THREADS})
do
    for CORE in $(seq 1 ${CORES})
    do
	ELEMENT=`echo "oper" | awk -v var1=${CORE} -v var2=${THREADS} -v var3=${THREAD} '{ print var1 * var2 + var3 }'`
	KMP_AFFINITY=${KMP_AFFINITY}${ELEMENT}","
    done
    KMP_AFFINITY=${KMP_AFFINITY}${THREAD}","

    if [ ${THREAD} -lt ${THREADS} ]; then
        ELEMENT=`echo "oper" | awk -v var1=${ELEMENT} -v var2=${THREADS} '{ print var1 + var2 }'`
        KMP_AFFINITY=${KMP_AFFINITY}${ELEMENT}","
    fi

done

KMP_AFFINITY=${KMP_AFFINITY}"0],explicit"

#echo "export KMP_AFFINITY="$KMP_AFFINITY

export KMP_AFFINITY


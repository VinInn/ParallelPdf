MY_ARGS="-e -n 500000"
NTHREADS=$1
if [ $# -lt 1 ]; then
    NTHREADS=32
fi

echo "Number of OpenMP threads: "$NTHREADS
COMMAND="cd /micfs/alazzaro; source sourcemic.sh; cd roofitMP; export OMP_NUM_THREADS=$NTHREADS; source sourceme_affinity_mic.sh; ./main $MY_ARGS"

echo $COMMAND
ssh mic0 $COMMAND

Run(){
  echo "$*"
  eval "$*"
}

prefix=$PWD
FLAGS="-no-vec"
# -fp-model precise -fp-model source -fimf-precision=high"

# add MIC options for default case (i.e. this script without a flag)
if [ ! "$*" ]; then
  MMIC="-mmic"
#  Run ./configure --prefix=$PWD/MIC --disable-openmp CXXFLAGS=\"-O2 -m64 -no-fma -Wall -DMKL -DWARNINGMSG $MMIC\" CXX=icpc CC=icc
#  Run ./configure --prefix=$PWD/MIC --disable-openmp CXXFLAGS=\"-O2 -m64 -no-fma -Wall -DWARNINGMSG $MMIC\" CXX=icpc CC=icc
  Run ./configure --prefix=$PWD/MIC --disable-openmp CXXFLAGS=\"-O2 -m64 -no-fma $FLAGS -Wall -DWARNINGMSG $MMIC\" CXX=icpc CC=icc

  # MIC case: correct libtool
  cp libtool libtool.MICerr
  sed -e 's+^postdeps=\"-L.*+postdeps=""+' libtool.MICerr > libtool
  Run diff libtool.MICerr libtool
else
#  Run ./configure --prefix=$PWD     --disable-openmp CXXFLAGS=\"-O2 -m64 -Wall -no-vec -DMKL -DWARNINGMSG\"       CXX=icpc CC=icc
#  Run ./configure --prefix=$PWD     --disable-openmp CXXFLAGS=\"-O2 -m64 -Wall -DMKL -DWARNINGMSG\"       CXX=icpc CC=icc
#  Run ./configure --prefix=$PWD     --disable-openmp CXXFLAGS=\"-O2 -m64 -Wall -mavx -DMKL -mkl=sequential -DWARNINGMSG\"       CXX=icpc CC=icc
#  Run ./configure --prefix=$PWD     --disable-openmp CXXFLAGS=\"-O2 -m64 -Wall -mavx -DMKL -DWARNINGMSG\"       CXX=icpc CC=icc
#  Run ./configure --prefix=$PWD     --disable-openmp CXXFLAGS=\"-O2 -m64 -Wall -mavx -DWARNINGMSG\"       CXX=icpc CC=icc
  Run ./configure --prefix=$PWD     --disable-openmp CXXFLAGS=\"-O2 -m64 -Wall $FLAGS -DWARNINGMSG\"       CXX=icpc CC=icc

fi


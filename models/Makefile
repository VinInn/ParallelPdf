
NAME     = modelRooFit
MICLABEL = "_mic"

MIC = $(shell root-config --cflags | grep "mmic" > /dev/null 2>&1 && echo $(MICLABEL))

target: icc

gcc: 
	g++ -O2 -m64 `root-config --libs --cflags` -lRooFitCore -lRooFit -lFoam -lMinuit -fopenmp $(NAME).cxx -o $(NAME)

icc: 
	icpc -O2 -m64 $(NAME).cxx -o $(NAME)$(MIC) `root-config --libs --cflags` -lRooFitCore -lRooFit -lFoam -lMinuit -openmp

turbotool:
	gcc turbotool.c -o turbotool

clean:
	rm -rf *~ modelRooFit modelRooFit$(MICLABEL) *.dot *.tar.bz2 turbotool

cleanlogs:
	rm -rf logs *.txt

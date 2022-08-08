EXE=analysis

CXX=/software/openmpi-intel/3.1.6/bin/mpicxx

CXXFLAGS=-Wall -O3 -DNDEBUG
CXXFLAGS+=`/home/apaul/beegfs/Zweig-Rule/HEPfit/build/HEPfit/bin/hepfit-config --cflags` 
CXXFLAGS+=-I/home/apaul/beegfs/Zweig-Rule/HEPfit/build/BAT/include -D_MPI -I/usr/include
CXXFLAGS+=`/home/apaul/opt/local/bin/root-config --cflags` 
CXXFLAGS+=`/home/apaul/opt/local/bin/gsl-config --cflags` 
CXXFLAGS+=-I ./src

LIBS=`/home/apaul/beegfs/Zweig-Rule/HEPfit/build/HEPfit/bin/hepfit-config --libs`
LIBS+=-L/home/apaul/beegfs/Zweig-Rule/HEPfit/build/BAT/lib -lBATmodels -lBAT  
LIBS+=`/home/apaul/opt/local/bin/root-config --libs` 
LIBS+=`/home/apaul/opt/local/bin/root-config --ldflags` 
LIBS+=`/home/apaul/opt/local/bin/gsl-config --libs` 
LIBS+=-Wl,-rpath,/home/apaul/beegfs/Zweig-Rule/HEPfit/build/BAT/lib -Wl,-rpath,/software/openmpi-intel/3.1.6/lib -Wl,-rpath,/home/apaul/opt/local/lib -Wl,-rpath,/home/apaul/opt/local/lib


all: $(EXE)

$(EXE): myModel_MCMC.o src/myModel.o src/myObservables.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f src/*.o *.o

distclean:
	rm -f $(EXE) src/*.o *.o


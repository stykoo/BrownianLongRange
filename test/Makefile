CXX=g++

MKLROOT=/home/lptl/intel/system_studio_2018/compilers_and_libraries_2018.2.199/linux/mkl/
MKL_CFLAGS=-DMKL_ILP64 -m64 -I${MKLROOT}/include
MKL_LDFLAGS=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group

CXXFLAGS=-g -pg -W -Wall -Ofast $(MKL_CFLAGS) -DUSE_MKL
LDFLAGS=-g -pg -lpthread -lm -ldl $(MKL_LDFLAGS)

all: test_change_alpha timing

ewald.o: ../src/ewald.cpp ../src/ewald.h
	$(CXX) -c ../src/ewald.cpp -o $@ $(CXXFLAGS) 

test_change_alpha: test_change_alpha.cpp ewald.o
	$(CXX) test_change_alpha.cpp ewald.o -o $@ $(CXXFLAGS) $(LDFLAGS)

timing: timing.cpp ewald.o
	$(CXX) timing.cpp ewald.o -o $@ $(CXXFLAGS) $(LDFLAGS)

clean:
	$(RM) *.o

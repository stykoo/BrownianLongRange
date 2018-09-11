CXX=g++
MKLROOT=/home/lptl/intel/system_studio_2018/mkl/
CXXFLAGS=-Og -g -Wall -Wshadow -Wextra -I${MKLROOT}/include -DUSE_MKL
MKLFLAGS=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -ldl

all: test_change_alpha timing

ewald.o: ../src/ewald.cpp ../src/ewald.h
	$(CXX) $(CXXFLAGS) -c ../src/ewald.cpp -o $@ $(MKLFLAGS)

test_change_alpha: test_change_alpha.cpp ewald.o
	$(CXX) $(CXXFLAGS) test_change_alpha.cpp ewald.o -o $@ $(MKLFLAGS)

timing: timing.cpp ewald.o
	$(CXX) $(CXXFLAGS) timing.cpp ewald.o -o $@ $(MKLFLAGS)

clean:
	$(RM) *.o test_change_alpha

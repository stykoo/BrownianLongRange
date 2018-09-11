#include <iostream>
#include <unistd.h>
#include <ctime>
#include <cmath>
#include "../src/ewald.h"

using namespace std;

void test(const Dim D, const int N, const int use_mkl) {
	const int n_tests = 10000; 
	const int n_change = 100; 
	std::cout << "D=" << D << ", N=" << N << ", use_mkl=" << use_mkl
		<< ", n_tests=" << n_tests << ", n_change=" << n_change << std::endl;

	double bias0 = 1.0;
	double bias1 = 1.5;

	double *positions = new double[D*N];
	double *charges = new double[N];

	for (int i = 0 ; i < D*N ; ++i) {
		positions[i] = drand48();
	}
	for (int i = 0 ; i < N ; ++i) {
		charges[i] = drand48();
	}

	double v0, v1;
	double *f0, *f1;

	Ewald *ew0 = new Ewald(D, N, positions, charges, bias0, use_mkl);
	Ewald *ew1 = new Ewald(D, N, positions, charges, bias1, use_mkl);

	f0 = ew0->update(v0);
	f1 = ew1->update(v1);

	//ew0->dump();
	//ew1->dump();
	
	double diff_max = 0., mx2_max = 0.;

	for (int k = 0 ; k < n_tests ; ++k) {
		if(k % n_change == 0){
			delete ew0;
			delete ew1;
			bias1 =  (0.25+2.5*drand48() )*bias0;
			//std::cout << "bias0=" << bias0 << ", bias1=" << bias1 << std::endl;
			//for (int i = 0 ; i < N ; ++i) {
			//	charges[i] = drand48();
			//}
			for (int i = 0 ; i < N ; ++i) {
				charges[i] = drand48();
			}
			ew0 = new Ewald(D, N, positions, charges, bias0, use_mkl);
			ew1 = new Ewald(D, N, positions, charges, bias1, use_mkl);
		}

		for (int i = 0 ; i < D*N ; ++i) {
			positions[i] = drand48();
		}
		f0 = ew0->update(v0);
		f1 = ew1->update(v1);

		double diff = fabs(v0 - v1);
		if (diff > diff_max) {
			diff_max = diff;
		}

		double mx2 = 0.0;
		for (int i = 0 ; i < D*N ; ++i) {
			mx2 += (f0[i] - f1[i]) * (f0[i] - f1[i]);
		}
		if (mx2 > mx2_max) {
			mx2_max = mx2;
		}
		//std::cout << "diff=" << diff << ", mx2=" << mx2 << std::endl;
	}

	std::cout << "==> diff_max=" << diff_max << ", mx2_max="
		<< mx2_max << std::endl;

	delete [] positions;
	delete [] charges;
}

int main() {
	test(DIM2, 2, 0);
	test(DIM2, 2, 1);
	test(DIM2, 2, 2);
	test(DIM3, 2, 0);
	test(DIM3, 2, 1);
	test(DIM3, 2, 2);

	return 0;
}

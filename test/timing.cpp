#include <iostream>
#include <unistd.h>
#include <ctime>
#include "../src/ewald.h"

using namespace std;

void test(const Dim D, const int N, const int n_iters, const double bias) {
	double *positions = new double[D * N];
	double *charges = new double[D * N];
	double v;

	for (int i = 0 ; i < D*N ; ++i) {
		charges[i] = 1.0;
		positions[i] = drand48();
	}

	Ewald ew1(D, N, positions, charges, bias, 0);
	clock_t beg = clock();
	for (int i = 0 ; i < n_iters ; ++i) {
		ew1.update(v);
	}
	double time1 = (clock() - beg) / (n_iters * (double) CLOCKS_PER_SEC);

	Ewald ew2(D, N, positions, charges, bias, 1);
	beg = clock();
	for (int i = 0 ; i < n_iters ; ++i) {
		ew2.update(v);
	}
	double time2 = (clock() - beg) / (n_iters * (double) CLOCKS_PER_SEC);

	Ewald ew3(D, N, positions, charges, bias, 2);
	beg = clock();
	for (int i = 0 ; i < n_iters ; ++i) {
		ew3.update(v);
	}
	double time3 = (clock() - beg) / (n_iters * (double) CLOCKS_PER_SEC);

	cout << N << " " << time1 << " " << time2 << " " << time3 << endl;

	delete[] positions;
	delete[] charges;
}

int main() {
	const int B = 4;
	double biases[B] = {1, 2, 3, 4};
	const int S = 7;
	int Ns[S] = {20, 50, 100, 200, 500, 1000, 2000};
	int n_iters = 5;

	for (int D = 2 ; D <= 3 ; ++D) {
		for (int i = 0 ; i < B ; ++i) {
			cout << "# Dim: " << D << ", Bias: " << biases[i]
				<< "\n# N no_mkl mkl mkl_full" << endl;

			for (int k = 0 ; k < S ; ++k) {
				test((Dim) D, Ns[k], n_iters, biases[i]);
			}

			cout << endl;
		}
	}

	return 0;
}

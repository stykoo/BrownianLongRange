/*
 * This file is an adaption of a code written by
 * Tony Maggs, PCT, Gulliver, ESPCI Paris.
 *
 * It implements the Evald method for the computation
 * of the Coulomb interaction between particles.
 */

//This is a conventional Ewald code with force with complexity N^{3/2}
//http://www.tcm.phy.cam.ac.uk/~mdt26/downloads/fraser.pdf
//Finite-size effects and Coulomb interactions in quantum Monte Carlo calculations
//for homogeneous systems with periodic boundary conditions
//Ewald with zero spatial average
//see also http://onlinelibrary.wiley.com/doi/10.1002/9780470164112.app6/pdf

#include <iostream>
#include <math.h>
#include <fstream>
#include <cassert>
#include "ewald.h"

#ifdef USE_MKL
	#include "mkl.h"
	#include "mkl_vsl.h"
#endif

using namespace std;

Ewald::Ewald(const int _N, double *_p, const double bias, const double *_Q,
	         const Dim _D) : D(_D), N(_N), p(_p), Q(_Q) {
	cout << "Ewald Error bound\t" << Ewald::error << endl;

	// Bias large to push to fourier space
	double twid = 1.7 * bias;
	alpha = sqrt(M_PI) * twid; // Equal work in real and fourier space
	alpha2 = alpha * alpha;
	double safety = sqrt(-log(Ewald::error));
	rRange = safety / alpha;
	fRange = safety * alpha / M_PI;

	cout << "\e[1;33mreal space Range " << rRange << " fourier Frange "
		<< fRange << " alpha " << alpha << "\e[0;03m" << endl;

	fr = new double[D*N]; // Realspace force
	ff = new double [D*N]; // Fourier force
	force = new double [D*N]; // Total force
	for(int i = 0 ; i < D*N ; i++){
		fr[i] = ff[i] = force[i] = 0;
	}

	initReal();
	initFourier();
}

void Ewald::initReal() {
	// Zero body potential
	vr0 = 0.0;
	for (int i = 0 ; i < N ; i++) {
		for (int j = 0 ; j < i ; j++) {
			vr0 += (Q[i] * Q[j]);
		}
	}
	vr0 *= -M_PI / alpha2;

#ifdef USE_MKL
	rs_hi = int(rRange + .5); // How many images in real space
	rs_lo = -rs_hi;
	rs_size = rs_hi - rs_lo + 1;
	if (D == DIM2) {
		rs_sizeTot = rs_size * rs_size;
	} else if (D == DIM3) {
		rs_sizeTot = rs_size * rs_size * rs_size;
	}
	std::cout << "rs_sizeTot: " << rs_sizeTot << std::endl;
	norms2 = new double[N * rs_sizeTot];
	norms = new double[N * rs_sizeTot];
	efs = new double[N * rs_sizeTot];
	cms = new double[N * rs_sizeTot];
	rrs = new double[D * N];
	in_range = new bool[N * rs_sizeTot];
#endif
}

void Ewald::initFourier() {
	lq = int(fRange + 0); //range in fourier space
	size = 2 * lq + 1;
	if (D == DIM2) {
		sizeTot = size * size;
	} else if (D == DIM3) {
		sizeTot = size * size * size;
	}
	cout << "Dim " << size << endl;

	a = new double[sizeTot]; // Fourier space sum
	c = new double[sizeTot]; // cos factor
	s = new double[sizeTot]; // sin factor
	qx = new double[sizeTot]; 
	qy = new double[sizeTot];
	qz = new double[sizeTot];
	for(int i = 0 ; i < sizeTot ; i++){
		qx[i] = qy[i] = qz[i] = a[i] = c[i] = s[i] =0; //fourier space arrays
	}

	cx = new double[size];
	cy = new double[size];
	cz = new double[size];

	sx = new double[size];
	sy = new double[size];
	sz = new double[size];

	// Build Fourier space weights
	for(int i = -lq ; i <= lq ; i++){
		int ii = lq + i;
		for(int j = -lq ; j <= lq ; j++){
			int jj = lq + j;
			if (D == DIM2) {
				double norm2 = i*i + j*j;
				int inx = ii + size*jj;
				if(norm2 != 0.) {
					a[inx] = 2 * exp(-M_PI * M_PI * norm2 / alpha2)
						/ (2 * M_PI * norm2);
				} else {
					a[inx] = 0;
				}

				// Wavevectors
				qx[inx] = 2 * M_PI * i;
				qy[inx] = 2 * M_PI * j;
			} else if (D == DIM3) {
				for(int k = -lq ; k <= lq ; k++){
					int kk = lq + k;
					double norm2 = i*i + j*j + k*k;
					int inx = ii + size*jj + size*size*kk;
					if(norm2 != 0.) {
						a[inx] = 2 * exp(-M_PI * M_PI * norm2 / alpha2)
							/ (2 * M_PI * norm2);
					} else {
						a[inx] = 0;
					}

					// Wavevectors
					qx[inx] = 2 * M_PI * i;
					qy[inx] = 2 * M_PI * j;
					qz[inx] = 2 * M_PI * k;
				}
			}
		}
	}

	suma = 0;
	for (int i = 0 ; i < sizeTot ; i++) {
		suma += a[i];
	}
}

Ewald::~Ewald() {
	delete [] fr; delete [] ff; delete [] force;

	delete [] a; delete [] c; delete [] s;

	delete [] qx; delete [] qy; delete [] qz;
	delete [] sx; delete [] sy; delete [] sz;
	delete [] cx; delete [] cy; delete [] cz;

#ifdef USE_MKL
	delete [] norms2;
	delete [] norms;
	delete [] efs;
	delete [] cms;
	delete [] rrs;
	delete [] in_range;
#endif
}

double* Ewald::fullforce(double * Ftot) {
	realSpace();
	fourierSpace();
	for(int i = 0 ; i < D*N ; i++){
		force[i] = fr[i] + ff[i];
	}
	*Ftot = vr + vf;
	return force;
}

void Ewald::dump(){
	cout << "\e[1;34mEnergy " << vr + vf <<"\e[0;39m"<<endl;
	cout << "Force" << " " << this << endl;
	for(int i=0 ; i < N ; i++){
		cout << "\e[1;31m"<<force[i] << "\t" << force[i+N] << "\t"
			<< force[i+2*N] << "\e[0;39m"<< endl;
	}
	cout << endl << endl;
}

void Ewald::realSpace() {
	for(int i = 0; i < D * N ; i++) {
		fr[i] = 0; //realspace force
	}

	vr = vr0; // Zero body potential

#ifdef USE_MKL
	for(int i = 0 ; i < N ; i++) {
		realSpaceAux2MKL(i);
	}
#else
	double vv;
	double f[3];
	// Double particle loop
	for(int i = 0 ; i < N ; i++) {
		for(int j = 0 ; j < i ; j++) {
			  realSpaceAux(i, j, &vv, f);

			  vr += vv;
			  
			  fr[i] += f[0];
			  fr[j] -= f[0];
			  
			  fr[i+N] += f[1];
			  fr[j+N] -= f[1];
			   
			  if (D == DIM3) {
				  fr[i+2*N] += f[2];
				  fr[j+2*N] -= f[2];
			  }
		}
	}
#endif
}

void Ewald::realSpaceAux(int ii, int jj, double *vv, double *f) {
	*vv = 0;
	double r[3];
	for (int l = 0 ; l < D ; l++) {
		f[l] = 0;
		r[l] = p[ii + N*l] - p[jj + N*l];
		if (r[l] > .5) {
			r[l] -= 1.0;
		} else if (r[l] < -.5) {
			r[l] += 1.0;
		}
	}

	int hi = int(rRange + .5); // How many images in real space
	int lo = -hi; 
	double qprod = Q[ii] * Q[jj];
	double pref = 2 * alpha / sqrt(M_PI);

	// Sum over periodic images
	for(int i=lo ; i <= hi ; i++){
		double d0 = r[0]-i;
		for(int j=lo; j<= hi; j++){
			double d1=r[1]-j;
			if (D == DIM2) {
				double norm2 = d0*d0 + d1*d1;
				if (norm2 !=0 && norm2 < rRange*rRange){
					double norm = sqrt(norm2);
					// Ewald energy
					double ef = qprod * erfc(alpha * norm) / norm;
					*vv += ef;
					// Forces
					double cm = pref * qprod * exp(-alpha2 * norm2) + ef;
					cm /= norm2; // Check formula
					f[0] += cm * d0;
					f[1] += cm * d1;
				}
			} else if (D == DIM3) {
				for(int k = lo; k <= hi; k++) {
					double d2 = r[2] - k;
					double norm2 = d0*d0 + d1*d1 + d2*d2;
					if (norm2 !=0 && norm2 < rRange*rRange){
						double norm = sqrt(norm2);
						// Ewald energy
						double ef = qprod * erfc(alpha * norm) / norm;
						*vv += ef;
						// Forces
						double cm = pref * qprod * exp(-alpha2 * norm2) + ef;
						cm /= norm2; // Check formula
						f[0] += cm * d0;
						f[1] += cm * d1;
						f[2] += cm * d2;
					}
				}
			}
		}
	}
}

#ifdef USE_MKL
void Ewald::realSpaceAux2MKL(int ii) {
	// Compute all the squared norms
	for (int jj = 0 ; jj < ii ; ++jj) {
		for (int l = 0 ; l < D ; l++) {
			double r = p[ii + N*l] - p[jj + N*l];
			if (r > .5) {
				r -= 1.0;
			} else if (r < -.5) {
				r += 1.0;
			}
			rrs[jj + N*l] = r;
		}
	}

	double d0i, d1i, d2i=0, d0, d1, d2;
	double *norms2_cur = norms2;
	for (int jj = 0 ; jj < ii ; ++jj) {
		d0i = rrs[jj] - rs_hi;
		d1i = rrs[jj+N] - rs_hi;
		if (D == DIM3) {
			d2i = rrs[jj+2*N] - rs_hi;
		}
		d0 = d0i;
		for(int i = 0 ; i < rs_size ; i++, d0 += 1.0){
			d1 = d1i;
			for(int j = 0 ; j < rs_size ; j++, d1 += 1.0){
				if (D == DIM2) {
					*norms2_cur++ = d0 * d0 + d1 * d1;
				} else if (D == DIM3) {
					d2 = d2i;
					for(int k = 0 ; k < rs_size ; k++, d2 += 1.0) {
						*norms2_cur++ = d0 * d0 + d1 * d1 + d2 * d2;
					}
				}
			}
		}
	}

	double rRange2 = rRange * rRange;
	double fac = 2 * alpha / sqrt(M_PI);

	norms2_cur = norms2;
	for (int p = 0 ; p < ii * rs_sizeTot ; ++p) {
		//if ((in_range[p] = std::signbit(norms2[p] - rRange2))) {
		if ((in_range[p] = (norms2[p] < rRange2))) {
			*norms2_cur++ = norms2[p];
		}
	}
	int inx = norms2_cur - norms2;

	// Use MKL functions on vectors
	vdSqrt(inx, norms2, norms);
	cblas_daxpby(inx, alpha, norms, 1, 0.0, efs, 1);
	vdErfc(inx, efs, efs);
	vdDiv(inx, efs, norms, efs);

	cblas_daxpby(inx, -alpha2, norms, 1, 0.0, cms, 1);
	vdExp(inx, cms, cms);
	cblas_daxpby(inx, 1.0, efs, 1, fac, cms, 1);
	vdDiv(inx, cms, norms2, cms);

	// Store the forces
	double f[3];
	bool *in_range_cur = in_range;
	double *efs_cur = efs;
	double *cms_cur = cms;
	for (int jj = 0 ; jj < ii ; ++jj) {
		double qprod = Q[ii] * Q[jj];
		for (int i = 0 ; i < D ; i++) {
			f[i] = 0;
		}

		d0i = rrs[jj] - rs_hi;
		d1i = rrs[jj+N] - rs_hi;
		if (D == DIM3) {
			d2i = rrs[jj+2*N] - rs_hi;
		}
		d0 = d0i;
		for(int i = 0 ; i < rs_size ; i++, d0 += 1.0){
			d1 = d1i;
			for(int j = 0 ; j < rs_size ; j++, d1 += 1.0){
				if (D == DIM2) {
					if (*in_range_cur++) {
						vr += qprod * (*efs_cur++);
						double cadd = qprod * (*cms_cur++);
						f[0] += cadd * d0;
						f[1] += cadd * d1;
					}
				} else if (D == DIM3) {
					d2 = d2i;
					for(int k = 0 ; k < rs_size ; k++, d2 += 1.0) {
						if (*in_range_cur++) {
							vr += qprod * (*efs_cur++);
							double cadd = qprod * (*cms_cur++);
							f[0] += cadd * d0;
							f[1] += cadd * d1;
							f[2] += cadd * d2;
						}
					}
				}
			}
		}

		fr[ii] += f[0];
		fr[jj] -= f[0];

		fr[ii+N] += f[1];
		fr[jj+N] -= f[1];

		if (D == DIM3) {
			fr[ii+2*N] += f[2];
			fr[jj+2*N] -= f[2];
		}
	}
}
#endif

double Ewald::fourierSpace() {
	for(int i = 0 ; i < D*N ; i++) {
		ff[i]=0; //fourier force
	}

	getStruct(); // Structure factor
	getPoten(); // Potential

	for (int i = 0 ; i < N ; i++){ //Force calcuation //particles loop
		arrays(i); // Fill cx[] and sx[] for particle i

		for(int ii = 0 ; ii < size ; ii++){
			for(int j = 0 ; j < size ; j++){
				if (D == DIM2) {
					int inx = ii + size*j;

					/* Taking k_z = 0 instead of summing over all k_z
					 * amounts to averaging the contributions on z.
					 * By symmetry, it does not change anything for the forces.
					 */
					// cz = 1, sz = 0
					double cc = cx[ii] * cy[j] - sx[ii] * sy[j]; 
					cc *= Q[i]; //charge
					// cz = 1, sz = 0
					double ss = sx[ii] * cy[j] + sy[j] * cx[ii];
					ss *= Q[i]; //charge
					double fact = 2 * a[inx] * (c[inx] * ss - s[inx] * cc);
					ff[i  ] += qx[inx] * fact;
					ff[i+N] += qy[inx] * fact;
				} else if (D == DIM3) {
					for(int k = 0 ; k < size ; k++){
						int inx = ii + size*j + size*size*k;

						double cc = cx[ii] * cy[j] * cz[k] 
							- cx[ii] * sy[j] * sz[k] 
							- sx[ii] * cy[j] * sz[k] 
							- sx[ii] * sy[j] * cz[k]; 
						cc *= Q[i]; //charge
						double ss = sx[ii] * cy[j] * cz[k]
							+ sy[j] * cx[ii] * cz[k]
							+ sz[k] * cx[ii] * cy[j]
							- sx[ii] * sy[j] * sz[k];
						ss *= Q[i]; //charge
						double fact = 2 * a[inx] * (c[inx] * ss - s[inx] * cc);
						ff[i    ] += qx[inx] * fact;
						ff[i+  N] += qy[inx] * fact;
						ff[i+2*N] += qz[inx] * fact;
					}
				}
			}
		}
	}

	return vf;
}

void Ewald::getPoten() {
	vf = -2 * N * suma; // Zero body
	for(int i = 0 ; i < sizeTot ; i++) { // Total energy in fourier space
		vf += (c[i] * c[i] + s[i] * s[i]) * a[i];
	}
}

void Ewald::getStruct(){
	for(int i = 0 ; i < sizeTot ; i++) { //zero arrays for structure factors
		c[i] = 0;
		s[i] = 0;
	}

	// Loop over particles
	for(int l = 0 ; l < N; l++){
		arrays(l); // cx[] sx[] for particle l
		for(int i = 0 ; i < size ; i++){
			for(int j = 0 ; j < size ; j++){
				if (D == DIM2) {
					int inx = i + size*j;
					// cz = 1, sz = 0
					double cc = cx[i] * cy[j] - sx[i] * sy[j]; 
					c[inx] += cc * Q[l];//charge 

					double ss = sx[i] * cy[j] + sy[j] * cx[i];
					s[inx] += ss * Q[l]; //charge
				} else if (D == DIM3) {
					for(int k = 0 ; k < size ; k++){
						int inx = i + size*j + size*size*k;
			  
						double cc = cx[i]* cy[j] * cz[k] 
							- cx[i] * sy[j] * sz[k] 
							- sx[i] * cy[j] * sz[k] 
							- sx[i] * sy[j] * cz[k]; 
						c[inx] += cc *Q[l]; //charge 

						double ss = sx[i] * cy[j] * cz[k]
							+ sy[j] * cx[i] * cz[k]
							+ sz[k] * cx[i] * cy[j]
							- sx[i] * sy[j] * sz[k];
						s[inx] += ss * Q[l]; //charge
					}
				}
			}
		}
	}
}

void Ewald::arrays(int l){
	for(int ii = 0 ; ii < size ; ii++) {
		double q = (-lq + ii) * 2 * M_PI; 

		cx[ii] = cos(q * p[l]);
		sx[ii] = sin(q * p[l]);
		cy[ii] = cos(q * p[l+N]);
		sy[ii] = sin(q * p[l+N]);
		if (D == DIM3) {
			cz[ii] = cos(q * p[l+2*N]);
			sz[ii] = sin(q * p[l+2*N]);
		}
	}
}

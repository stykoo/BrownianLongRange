/*
Copyright (C) ESPCI Paris / Sorbonne Université (2018)
Contributors: Tony Maggs, Alexis Poncet <aponcet@lptmc.jussieu.fr>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ActiveBrownian is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file ewald.cpp
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Ewald method for Coulomb interactions
 *
 * This is a conventional Ewald code with force with complexity N^{3/2}
 * http://www.tcm.phy.cam.ac.uk/~mdt26/downloads/fraser.pdf
 * Finite-size effects and Coulomb interactions in quantum Monte Carlo
 * calculations for homogeneous systems with periodic boundary conditions
 * Ewald with zero spatial average
 * see also http://onlinelibrary.wiley.com/doi/10.1002/9780470164112.app6/pdf
 */

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include "ewald.h"

#ifdef USE_MKL
#include "mkl.h"
#include "mkl_vsl.h"
#endif

using namespace std;

Ewald::Ewald(const Dim D_, const int N_, double *p_, const double *Q_,
			 const double bias_, int use_mkl_) : D(D_), N(N_) {
#ifndef USE_MKL
	assert(use_mkl_ == 0);
#endif
	//cout << "Ewald Error bound\t" << Ewald::error << endl;

	// Bias large to push to fourier space
	double twid = 1.7 * bias_;
	double alpha = sqrt(M_PI) * twid; // Equal work in real and fourier space
	double safety = sqrt(-log(Ewald::error));
	double rRange = safety / alpha;
	double fRange = safety * alpha / M_PI;

	//cout << "\e[1;33mreal space Range " << rRange << " fourier Frange "
		//<< fRange << " alpha " << alpha << "\e[0;03m" << endl;

	// We copy the charges because it would be bad if they changed from outside
	local_q = new double[N];
	for (int i = 0 ; i < N ; ++i) {
		local_q[i] = Q_[i];
	}

	energy = 0.0;
	forces = new double[D*N]; // Total force
	for (int i = 0 ; i < D*N ; ++i) {
		forces[i] = 0.0;
	}

	ewReal = new EwaldRS(D, N, p_, local_q, alpha, rRange, use_mkl_);
	ewFourier = new EwaldFS(D, N, p_, local_q, alpha, fRange);
}

Ewald::~Ewald() {
	delete [] local_q;
	delete [] forces;
	delete ewReal;
	delete ewFourier;
}

double* Ewald::update(double &vv) {
	ewReal->update();
	ewFourier->update();

	energy = ewReal->getEnergy() + ewFourier->getEnergy();

	double *fr = ewReal->getForces();
	double *ff = ewFourier->getForces();
	for(int i = 0 ; i < D*N ; i++){
		forces[i] = fr[i] + ff[i];
	}
	vv = energy;
	return forces;
}

void Ewald::dump() const {
	cout << "\e[1;34mEnergy " << ewReal->getEnergy() << " + "
	    << ewFourier->getEnergy() << " = " << energy <<"\e[0;39m"<<endl;
	cout << "Forces" << endl;
	for(int i=0 ; i < N ; i++){
		cout << "\e[1;31m";
		for (int l = 0 ; l < D ; ++l) {
			cout << forces[i + N*l] << "\t";
		}
		cout << "\e[0;39m"<< endl;
	}
	cout << endl << endl;
}


EwaldRS::EwaldRS(const Dim D_, const int N_, double *p_, const double *Q_, 
		         const double alpha_, const double range_, int use_mkl_) :
		D(D_), N(N_), p(p_), Q(Q_), alpha(alpha_), range(range_),
		use_mkl(use_mkl_) {
	energy = 0.0;
	forces = new double[D*N];
	for (int i = 0 ; i < D*N ; ++i) {
		forces[i] = 0.0;
	}

	// Zero body potential
	zero_body = 0.0;
	for (int i = 0 ; i < N ; i++) {
		for (int j = 0 ; j < i ; j++) {
			zero_body += (Q[i] * Q[j]);
		}
	}
	if (D == DIM2) {
		zero_body *= -2 * sqrt(M_PI) / alpha;
	} else if (D == DIM3) {
		zero_body *= -M_PI / (alpha * alpha);
	}

	shift_max = int(range + .5); // How many images in real space
	n_images_x = 2 * shift_max + 1;
	n_images_tot = 1;
	for (int l = 0 ; l < D ; ++l) {
		n_images_tot *= n_images_x;
	}
	//std::cout << "Real space images: " << n_images_tot << std::endl;

#ifdef USE_MKL
	if (use_mkl > 0) {
		int nn = 1;
		if (use_mkl == 1) {
			nn = N;
		} else if (use_mkl == 2) {
			nn = N * (N - 1) / 2;
		}
		int S = n_images_tot * nn;
		norms2 = new double[S];
		norms = new double[S];
		efs = new double[S];
		cms = new double[S];
		rrs = new double[D * nn];
		in_range = new bool[S];
	}
#endif
}

EwaldRS::~EwaldRS() {
	delete [] forces;
#ifdef USE_MKL
	if (use_mkl > 0) {
		delete [] norms2; delete [] norms;
		delete [] efs; delete [] cms; delete [] rrs;
		delete [] in_range;
	}
#endif
}

void EwaldRS::update() {
	for(int i = 0; i < D * N ; i++) {
		forces[i] = 0;
	}

	energy = zero_body; // Zero body energy

	if (use_mkl == 0) {
		for(int i = 0 ; i < N ; i++) {
			for(int j = 0 ; j < i ; j++) {
				updateAux(i, j);
			}
		}
	}
#ifdef USE_MKL
	else if (use_mkl == 1) {
		for(int i = 0 ; i < N ; i++) {
			updateAuxMKL1(i);
		}
	} else if (use_mkl == 2) {
		updateAuxMKL2();
	}
#endif
}

void EwaldRS::updateAux(const int ii, const int jj) {
	double vv = 0;
	double f[3], r[3];
	for (int l = 0 ; l < D ; l++) {
		f[l] = 0;
		r[l] = p[ii + N*l] - p[jj + N*l];
		r[l] -= std::round(r[l]);
	}

	double qprod = Q[ii] * Q[jj];
	double range2 = range * range;
	double alpha2 = alpha * alpha;
	double pref = 2 * alpha / sqrt(M_PI);

	// Sum over periodic images
	for(int i=-shift_max ; i <= shift_max ; i++){
		double d0 = r[0]-i;
		for(int j=-shift_max; j<= shift_max; j++){
			double d1=r[1]-j;
			if (D == DIM2) {
				double norm2 = d0*d0 + d1*d1;
				if (norm2 !=0 && norm2 < range2){
					double norm = sqrt(norm2);
					// Ewald energy
					double ef = qprod * erfc(alpha * norm) / norm;
					vv += ef;
					// Forces
					double cm = pref * qprod * exp(-alpha2 * norm2) + ef;
					cm /= norm2; // Check formula
					f[0] += cm * d0;
					f[1] += cm * d1;
				}
			} else if (D == DIM3) {
				for(int k = -shift_max; k <= shift_max; k++) {
					double d2 = r[2] - k;
					double norm2 = d0*d0 + d1*d1 + d2*d2;
					if (norm2 !=0 && norm2 < range2){
						double norm = sqrt(norm2);
						// Ewald energy
						double ef = qprod * erfc(alpha * norm) / norm;
						vv += ef;
						// Forces
						double cm = pref * qprod * exp(-alpha2*norm2) + ef;
						cm /= norm2; // Check formula
						f[0] += cm * d0;
						f[1] += cm * d1;
						f[2] += cm * d2;
					}
				}
			}
		}
	}

	energy += vv;
	for (int l = 0 ; l < D ; ++l) {
		forces[ii+l*N] += f[l];
		forces[jj+l*N] -= f[l];
	}
}

#ifdef USE_MKL
void EwaldRS::updateAuxMKL1(int ii) {
	// Compute all the squared norms
	for (int jj = 0 ; jj < ii ; ++jj) {
		for (int l = 0 ; l < D ; l++) {
			double r = p[ii + N*l] - p[jj + N*l];
			r -= std::round(r);
			rrs[jj + N*l] = r;
		}
	}

	double dd[3], dd0[3];
	double *norms2_cur = norms2;
	for (int jj = 0 ; jj < ii ; ++jj) {
		for (int l = 0 ; l < D ; ++l) {
			dd0[l] = rrs[jj+N*l] - shift_max;
		}
		dd[0] = dd0[0];
		for(int i = 0 ; i < n_images_x ; i++, dd[0] += 1.0){
			dd[1] = dd0[1];
			for(int j = 0 ; j < n_images_x ; j++, dd[1] += 1.0){
				if (D == DIM2) {
					*norms2_cur++ = dd[0] * dd[0] + dd[1] * dd[1];
				} else if (D == DIM3) {
					dd[2] = dd0[2];
					for(int k = 0 ; k < n_images_x ; k++, dd[2] += 1.0) {
						*norms2_cur++ = (dd[0] * dd[0] + dd[1] * dd[1]
						                 + dd[2] * dd[2]);
					}
				}
			}
		}
	}

	double range2 = range * range;

	norms2_cur = norms2;
	for (int k = 0 ; k < ii * n_images_tot ; ++k) {
		if ((in_range[k] = (norms2[k] < range2))) {
			*norms2_cur++ = norms2[k];
		}
	}

	// Use MKL functions on vectors
	int nn = norms2_cur - norms2;
	opsMKL(nn);

	// Store the forces
	double f[3];
	bool *in_range_cur = in_range;
	double *efs_cur = efs;
	double *cms_cur = cms;
	for (int jj = 0 ; jj < ii ; ++jj) {
		double qprod = Q[ii] * Q[jj];
		for (int l = 0 ; l < D ; ++l) {
			f[l] = 0;
		}

		for (int l = 0 ; l < D ; ++l) {
			dd0[l] = rrs[jj+N*l] - shift_max;
		}
		dd[0] = dd0[0];
		for(int i = 0 ; i < n_images_x ; i++, dd[0] += 1.0){
			dd[1] = dd0[1];
			for(int j = 0 ; j < n_images_x ; j++, dd[1] += 1.0){
				if (D == DIM2) {
					if (*in_range_cur++) {
						energy += qprod * (*efs_cur++);
						double cadd = qprod * (*cms_cur++);
						f[0] += cadd * dd[0];
						f[1] += cadd * dd[1];
					}
				} else if (D == DIM3) {
					dd[2] = dd0[2];
					for(int k = 0 ; k < n_images_x ; k++, dd[2] += 1.0) {
						if (*in_range_cur++) {
							energy += qprod * (*efs_cur++);
							double cadd = qprod * (*cms_cur++);
							f[0] += cadd * dd[0];
							f[1] += cadd * dd[1];
							f[2] += cadd * dd[2];
						}
					}
				}
			}
		}

		for (int l = 0 ; l < D ; ++l) {
			forces[ii+N*l] += f[l];
			forces[jj+N*l] -= f[l];
		}
	}
}

void EwaldRS::updateAuxMKL2() {
	// Compute all the squared norms
	double *rrs_curr = rrs;
	for(int ii = 0 ; ii < N ; ii++) {
		for (int jj = 0 ; jj < ii ; ++jj) {
			for (int l = 0 ; l < D ; l++) {
				double r = p[ii + N*l] - p[jj + N*l];
				r -= std::round(r);
				*rrs_curr++ = r;
			}
		}
	}

	double dd[3], dd0[3];
	double *norms2_cur = norms2;
	rrs_curr = rrs;
	for(int ii = 0 ; ii < N ; ii++) {
		for (int jj = 0 ; jj < ii ; ++jj) {
			for (int l = 0 ; l < D ; ++l) {
				dd0[l] = *rrs_curr++ - shift_max;
			}
			dd[0] = dd0[0];
			for(int i = 0 ; i < n_images_x ; i++, dd[0] += 1.0){
				dd[1] = dd0[1];
				for(int j = 0 ; j < n_images_x ; j++, dd[1] += 1.0){
					if (D == DIM2) {
						*norms2_cur++ = dd[0] * dd[0] + dd[1] * dd[1];
					} else if (D == DIM3) {
						dd[2] = dd0[2];
						for(int k = 0 ; k < n_images_x ; k++, dd[2] += 1.0) {
							*norms2_cur++ = (dd[0] * dd[0] + dd[1] * dd[1]
							                 + dd[2] * dd[2]);
						}
					}
				}
			}
		}
	}

	double range2 = range * range;
	int Sall = n_images_tot * N * (N-1) / 2;

	norms2_cur = norms2;
	for (int k = 0 ; k < Sall ; ++k) {
		if ((in_range[k] = (norms2[k] < range2))) {
			*norms2_cur++ = norms2[k];
		}
	}

	// Use MKL functions on vectors
	int nn = norms2_cur - norms2;
	opsMKL(nn);

	// Store the forces
	double f[3];
	bool *in_range_cur = in_range;
	double *efs_cur = efs;
	double *cms_cur = cms;
	rrs_curr = rrs;
	for(int ii = 0 ; ii < N ; ii++) {
		for (int jj = 0 ; jj < ii ; ++jj) {
			double qprod = Q[ii] * Q[jj];
			for (int i = 0 ; i < D ; i++) {
				f[i] = 0;
			}

			for (int l = 0 ; l < D ; ++l) {
				dd0[l] = *rrs_curr++ - shift_max;
			}
			dd[0] = dd0[0];
			for(int i = 0 ; i < n_images_x ; i++, dd[0] += 1.0) {
				dd[1] = dd0[1];
				for(int j = 0 ; j < n_images_x ; j++, dd[1] += 1.0) {
					if (D == DIM2) {
						if (*in_range_cur++) {
							energy += qprod * (*efs_cur++);
							double cadd = qprod * (*cms_cur++);
							f[0] += cadd * dd[0];
							f[1] += cadd * dd[1];
						}
					} else if (D == DIM3) {
						dd[2] = dd0[2];
						for(int k = 0 ; k < n_images_x ; k++, dd[2] += 1.0) {
							if (*in_range_cur++) {
								energy += qprod * (*efs_cur++);
								double cadd = qprod * (*cms_cur++);
								f[0] += cadd * dd[0];
								f[1] += cadd * dd[1];
								f[2] += cadd * dd[2];
							}
						}
					}
				}
			}

			for (int l = 0 ; l < D ; ++l) {
				forces[ii+N*l] += f[l];
				forces[jj+N*l] -= f[l];
			}
		}
	}
}

void EwaldRS::opsMKL(const int n) {
	double fac = 2 * alpha / sqrt(M_PI);
	double alpha2 = alpha * alpha;

	vdSqrt(n, norms2, norms);
	cblas_daxpby(n, alpha, norms, 1, 0.0, efs, 1);
	vdErfc(n, efs, efs);
	vdDiv(n, efs, norms, efs);

	cblas_daxpby(n, -alpha2, norms2, 1, 0.0, cms, 1);
	vdExp(n, cms, cms);
	cblas_daxpby(n, 1.0, efs, 1, fac, cms, 1);
	vdDiv(n, cms, norms2, cms);
}
#endif


EwaldFS::EwaldFS(const Dim D_, const int N_, double *p_, const double *q_,
		         const double alpha_, const double range_) :
		D(D_), N(N_), p(p_), Q(q_), alpha(alpha_), range(range_) {
	lq = int(range + 0); // Range in fourier space
	size = 2 * lq + 1;
	sizeTot = 1;
	for (int l = 0 ; l < D ; ++l) {
		sizeTot *= size;
	}
	//cout << "Size " << size << endl;

	energy = 0.0;
	forces = new double[D*N];
	for (int i = 0 ; i < D*N ; ++i) {
		forces[i] = 0.0;
	}

	a = new double[sizeTot]; // Fourier space sum
	c = new double[sizeTot]; // cos factor
	s = new double[sizeTot]; // sin factor
	qx = new double[sizeTot]; 
	qy = new double[sizeTot];
	qz = new double[sizeTot];
	for (int i = 0 ; i < sizeTot ; ++i) {
		a[i] = c[i] = s[i] = qx[i] = qy[i] = qz[i] = 0.0;
	}


	cx = new double[size]; cy = new double[size]; cz = new double[size];
	sx = new double[size]; sy = new double[size]; sz = new double[size];

	double alpha2 = alpha * alpha;

	// Build Fourier space weights
	for(int i = -lq ; i <= lq ; i++){
		int ii = lq + i;
		for(int j = -lq ; j <= lq ; j++){
			int jj = lq + j;
			if (D == DIM2) {
				double norm2 = i*i + j*j;
				double norm = sqrt(norm2);
				int inx = ii + size*jj;
				if(norm != 0.) {
					a[inx] = erfc(M_PI * norm / alpha) / (2 * norm);
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
						a[inx] = exp(-M_PI * M_PI * norm2 / alpha2)
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

	double suma = 0;
	for (int i = 0 ; i < sizeTot ; i++) {
		suma += a[i];
	}
	zero_body = 0;
	for (int i = 0 ; i < N ; i++) {
		zero_body += (Q[i] * Q[i]);
	}
	zero_body *= -suma;
}

EwaldFS::~EwaldFS() {
	delete [] forces;
	delete [] a; delete [] c; delete [] s;
	delete [] qx; delete [] qy; delete [] qz;
	delete [] sx; delete [] sy; delete [] sz;
	delete [] cx; delete [] cy; delete [] cz;
}

void EwaldFS::update() {
	for(int i = 0 ; i < D*N ; i++) {
		forces[i]=0; //fourier force
	}

	computeStructureFactor(); // Structure factor
	computeEnergy(); // Potential

	for (int i = 0 ; i < N ; i++){ //Force calcuation //particles loop
		fillCosSin(i); // Fill cx[] and sx[] for particle i

		for(int ii = 0 ; ii < size ; ii++){
			for(int j = 0 ; j < size ; j++){
				if (D == DIM2) {
					int inx = ii + size*j;

					// cos(k0*ri0 + k1*rj1)
					double cc = cx[ii] * cy[j] - sx[ii] * sy[j]; 
					cc *= Q[i]; //charge
					// sin(k0*ri0 + k1*rj1)
					double ss = sx[ii] * cy[j] + sy[j] * cx[ii];
					ss *= Q[i]; //charge
					// ... sin()
					double fact = 2 * a[inx] * (c[inx] * ss - s[inx] * cc);
					forces[i  ] += qx[inx] * fact;
					forces[i+N] += qy[inx] * fact;
				} else if (D == DIM3) {
					for(int k = 0 ; k < size ; k++){
						int inx = ii + size*j + size*size*k;

						// cos(k0*ri0 + k1*rj1 + k2*rk2)
						double cc = cx[ii] * cy[j] * cz[k] 
							- cx[ii] * sy[j] * sz[k] 
							- sx[ii] * cy[j] * sz[k] 
							- sx[ii] * sy[j] * cz[k]; 
						cc *= Q[i]; //charge
						// sin(k0*ri0 + k1*rj1 + k2*rk2)
						double ss = sx[ii] * cy[j] * cz[k]
							+ sy[j] * cx[ii] * cz[k]
							+ sz[k] * cx[ii] * cy[j]
							- sx[ii] * sy[j] * sz[k];
						ss *= Q[i]; //charge
						double fact = 2 * a[inx]
									  * (c[inx] * ss - s[inx] * cc);
						forces[i    ] += qx[inx] * fact;
						forces[i+  N] += qy[inx] * fact;
						forces[i+2*N] += qz[inx] * fact;
					}
				}
			}
		}
	}
}

void EwaldFS::computeEnergy() {
	energy = zero_body; // Zero body
	for(int i = 0 ; i < sizeTot ; i++) { // Total energy in fourier space
		energy += (c[i] * c[i] + s[i] * s[i]) * a[i];
	}
}

void EwaldFS::computeStructureFactor(){
	for(int i = 0 ; i < sizeTot ; i++) { //zero arrays for structure factors
		c[i] = 0; s[i] = 0;
	}

	// Loop over particles
	for(int l = 0 ; l < N; l++){
		fillCosSin(l); // cx[] sx[] for particle l
		for(int i = 0 ; i < size ; i++){
			for(int j = 0 ; j < size ; j++){
				if (D == DIM2) {
					int inx = i + size*j;
					
					// cos(k0*r0 + k1*r1)
					double cc = cx[i] * cy[j] - sx[i] * sy[j]; 
					c[inx] += cc * Q[l];//charge 
					// sin(k0*r0 + k1*r1)
					double ss = sx[i] * cy[j] + sy[j] * cx[i];
					s[inx] += ss * Q[l]; //charge
				} else if (D == DIM3) {
					for(int k = 0 ; k < size ; k++){
						int inx = i + size*j + size*size*k;
			  
						// cos(k0*r0 + k1*r1 + k2*r2)
						double cc = cx[i] * cy[j] * cz[k] 
							- cx[i] * sy[j] * sz[k] 
							- sx[i] * cy[j] * sz[k] 
							- sx[i] * sy[j] * cz[k]; 
						c[inx] += cc * Q[l]; //charge 

						// sin(k0*r0 + k1*r1 + k2*r2)
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

void EwaldFS::fillCosSin(int l){
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

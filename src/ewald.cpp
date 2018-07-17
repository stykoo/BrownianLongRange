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

using namespace std;

#define PRINT 0

Ewald::Ewald(int _N, double *_p, double bias, double *_Q) :
		N(_N), p(_p), Q(_Q) {
	// Bias large to push to fourier space
	double twid = 1.7 * bias;
	double error = 1.e-8;

	cout << "Ewald Error bound\t" << error << endl;

	alpha = sqrt(M_PI) * twid; // Equal work in real and fourier space
	alpha2 = alpha * alpha;
	double safety = sqrt(-log(error));
	rRange = safety / alpha;
	fRange = safety * alpha / M_PI;

	cout << "\e[1;33mreal space Range " << rRange << " fourier Frange "
		<< fRange << " alpha " << alpha << "\e[0;03m" << endl;

	lq=int(fRange + 0); //range in fourier space
	dim = 2 * lq + 1;
	dim3 = dim * dim * dim;
	cout << "Dim " << dim << endl;

	fr = new double[3*N]; // Realspace force
	ff = new double [3*N]; // Fourier force
	force = new double [3*N]; // Fourier force (??)
	for(int i = 0 ; i < 3*N ; i++){
		fr[i] = ff[i] = force[i] = 0; // Real space arrays
	}

	a = new double[dim3]; // Fourier space sum
	c = new double[dim3]; // cos factor
	s = new double[dim3]; // sin factor
	qx = new double[dim3]; 
	qy = new double[dim3];
	qz = new double[dim3];
	for(int i = 0 ; i < dim3 ; i++){
		qx[i] = qy[i] = qz[i] = a[i] = c[i] = s[i] =0; //fourier space arrays
	}

	cx = new double[dim];
	cy = new double[dim];
	cz = new double[dim];

	sx = new double[dim];
	sy = new double[dim];
	sz = new double[dim];

	// Build Fourier space weights
	for(int i = -lq ; i <= lq ; i++){
		int ii = lq + i;
		for(int j = -lq ; j <= lq ; j++){
			int jj = lq + j;
	  		for(int k = -lq ; k <= lq ; k++){
				int kk = lq + k;
				double norm2 = i*i + j*j + k*k;
				int inx = ii + dim*jj + dim*dim*kk;
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

	suma = 0;
	for (int i = 0 ; i < dim3 ; i++) {
		suma += a[i];
	}

	// Zero body potential
	vr0 = 0.0;
	for (int i = 0 ; i < N ; i++) {
		for (int j = 0 ; j < i ; j++) {
			vr0 += (Q[i] * Q[j]);
		}
	}
	vr0 *= -M_PI / alpha2;

	if(PRINT){
		cout << "Qx\n";
		for(int i = 0 ; i < dim3 ; i++) {
			cout << qx[i] << " ";
		}
		cout << endl;
	}
}

Ewald::~Ewald() {
	delete [] fr;
	delete [] ff;
	delete [] force;

	delete [] a;
	delete [] c;
	delete [] s;

	delete [] qx;
	delete [] qy;
	delete [] qz;
	delete [] sx;
	delete [] sy;
	delete [] sz;
	delete [] cx;
	delete [] cy;
	delete [] cz;
}

double* Ewald::fullforce(double * Ftot) {
	realSpace();
	fourierSpace();
	for(int i = 0 ; i < 3*N ; i++){
		force[i] = fr[i] + ff[i];
	}
	*Ftot = vr + vf;
	return force;
}

/* void Ewald::samplefourier(int ns) {
	ofstream fs;
	fs.open("fourier.dat");
	fs.precision(15);
	assert(N==1);
	p[0]=0; p[1]=0;
	p[2] =0; p[3]=0;
	p[4]=sz/2.; p[5] = sz/2.;
	for (int i=0;i<=ns;i++){
		for(int j=0;j<=ns;j++){
			p[1] = double(i)/double(ns);
			p[3] = double(j)/double(ns);
			getStruct();
			getPoten();
			fs<<i<<" "<<j<<" "<<vf<<endl;
		}
	}
	fs.close();
} */

void Ewald::dump(){
	cout << "\e[1;34mEnergy " << vr + vf <<"\e[0;39m"<<endl;
	cout << "Force" << " " << this << endl;
	for(int i=0 ; i < 2*N ; i++){
		cout << "\e[1;31m"<<force[i] << "\t" << force[i+2*N] << "\t"
			<< force[i+4*N] << "\e[0;39m"<< endl;
	}
	cout << endl << endl;
}

void Ewald::realSpace() {
	double f[3];
	double vv;

	for(int i = 0; i < 3 * N ; i++) {
		fr[i]=0; //realspace force
	}

	vr = vr0; // Zero body potential

	// Double particle loop
	for(int i = 0 ; i < N ; i++) {
		for(int j = 0; j < i ; j++) {
			  realSpaceAux(i, j, &vv, f);
			  vr += vv;
			  
			  fr[i] += f[0];
			  fr[j] -= f[0];
			  
			  fr[i+N] += f[1];
			  fr[j+N] -= f[1];
			   
			  fr[i+2*N] += f[2];
			  fr[j+2*N] -= f[2];
		}
	}
}

void Ewald::realSpaceAux(int ii, int jj, double * vv, double *f) {
	*vv = 0;
	for (int i = 0 ; i < 3 ; i++) {
		f[i] = 0;
	}

	double r[3];
	for (int l = 0 ; l < 3 ; l++) {
		r[l] = p[ii + N*l] - p[jj + N*l];
		if (r[l] > .5) {
			r[l] -= 1.0;
		} else if (r[l] < -.5) {
			r[l] += 1.0;
		}
	}

	int hi = int(rRange + .5); // How many images in real space
	int lo = -hi;

	// Sum over periodic images
	for(int i=lo ; i <= hi ; i++ ){
		double d0 = r[0]-i;
		for(int j=lo; j<= hi; j++ ){
			double d1=r[1]-j;
			for(int k = lo; k <= hi; k++) {
				double d2 = r[2] - k;
				double norm2 = d0*d0 + d1*d1 + d2*d2;
				if (norm2 !=0 && norm2 < rRange*rRange){
					double norm = sqrt(norm2);
					// Ewald energy
					double ef = Q[ii] * Q[jj] * erfc(alpha * norm) / norm;
					*vv += ef;
					// Forces
					double cm = (2 * alpha / sqrt(M_PI)) * exp(-alpha2 * norm2)
						* Q[ii] * Q[jj] + ef;
					f[0] += cm * d0 / norm2; //check formula
					f[1] += cm * d1 / norm2;
					f[2] += cm * d2 / norm2;
				}
			}
		}
	}
}

double Ewald::fourierSpace() {
	for(int i = 0 ; i < 3*N ; i++) {
		ff[i]=0; //fourier force
	}

	getStruct(); // Structure factor
	getPoten(); // Potential

	for (int i = 0 ; i < 2*N ; i++){ //Force calcuation //particles loop
		arrays(i); // Fill cx[] and sx[] for particle i

		for(int ii = 0 ; ii < dim ; ii++){
			for(int j = 0 ; j < dim ; j++){
				for(int k = 0 ; k < dim ; k++){
					int inx = ii + dim*j + dim*dim*k;

					double cc = cx[ii]* cy[j] *cz[k] 
						- cx[ii] * sy[j] * sz[k] 
						- sx[ii] * cy[j] * sz[k] 
						- sx[ii] * sy[j] * cz[k]; 
					cc *= Q[i]; //charge
					double ss = sx[ii] * cy[j] * cz[k]
						+  sy[j] * cx[ii] * cz[k]
						+ sz[k] * cx[ii] * cy[j]
						- sx[ii]*sy[j]*sz[k];
					ss *= Q[i]; //charge
					double fact = 2* a[inx] * (c[inx] * ss - s[inx] * cc);
					ff[i    ] += qx[inx] * fact;
					ff[i+2*N] += qy[inx] * fact;
					ff[i+4*N] += qz[inx] * fact;
				}
			}
		}
	}

	return vf;
}

void Ewald::getPoten() {
	vf = -2 * N * suma; // Zero body
	for(int i = 0 ; i < dim3 ; i++){ // Total energy in fourier space
		double cc = c[i];
		double ss = s[i];
		vf += (cc*cc + ss*ss) * a[i];
	}
}

void Ewald::getStruct(){
	for(int i = 0; i < dim3; i++) { //zero arrays for structure factors
		c[i] = 0;
		s[i] = 0;
	}

	// Loop over particles
	for(int l = 0 ; l < N; l++){
		arrays(l); // cx[] sx[] for particle l
		for(int i = 0 ; i < dim ; i++){
			for(int j = 0 ; j < dim ; j++){
				for(int k = 0 ; k < dim ; k++){
					int inx = i + dim*j + dim*dim*k;
		  
					double cc = cx[i]* cy[j] *cz[k] 
						- cx[i] * sy[j] * sz[k] 
						- sx[i] * cy[j] * sz[k] 
						- sx[i] * sy[j] * cz[k]; 
						c[inx] += cc *Q[l];//charge 

					double ss = sx[i] * cy[j] * cz[k]
						+ sy[j] * cx[i] * cz[k]
						+ sz[k] * cx[i] * cy[j]
						- sx[i]*sy[j]*sz[k];
					s[inx] += ss*Q[l]; //charge
				}
			}
		}
	}
}

void Ewald::arrays(int l){
	for(int ii = 0 ; ii < dim ; ii++) {
		double q = (-lq + ii) * 2 * M_PI; 

		cx[ii] = cos(q * p[l]);
		sx[ii] = sin(q * p[l]);
		cy[ii] = cos(q * p[l+N]);
		sy[ii] = sin(q * p[l+N]);
		cz[ii] = cos(q * p[l+2*N]);
		sz[ii] = sin(q * p[l+2*N]);
	}
}

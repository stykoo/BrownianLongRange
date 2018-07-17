/*
 * This file is an adaption of a code written by
 * Tony Maggs, PCT, Gulliver, ESPCI Paris.
 *
 * It implements the Evald method for the computation
 * of the Coulomb interaction between particles.
 */


#ifndef BROWNIANLONGRANGE_EVALD_H_
#define BROWNIANLONGRANGE_EVALD_H_

#include <math.h>

class Ewald{
	public:
		Ewald(int _N, double *_p, double bias, double *_q);
		~Ewald();

		double* fullforce(double *);
		//void samplefourier(int ns);
		void dump();
		double* getf() { return force; }

	private:
		void setup(double, double, double);
		void realSpace();
		void realSpaceAux(int i, int j, double *v, double *f);
		double fourierSpace();
		void getPoten();
		void getStruct();
		void arrays(int i);

		int N; //!< Number of particles
		double *p; //!< Positions of the particles
		double *Q; //!< Particle charges
		double alpha; //!< Parameter of Evald method
		double alpha2; //!< Square of parameter of Evald method
		double rRange; //!< Range in real space
		double fRange; //!< Range in Fourier space
		int lq; //!< Integer range in Fourier space
		int dim; //!< Size of arrays in 1d Fourier space
		int dim3; //!< Size of arrays in 3d Fourier space
		double *ff; //!< Forces in Fourier space
		double *fr; //!< Forces in real space
	    double *force; //!< Total forces
		double *a; //!< Fourier space sum
		double *c; //!< Cos factors
		double *s; //!< Sin factors
		double *qx, *qy, *qz; //!< Wavevectors
		double suma; //!< Sum of the Fourier space contributions
		double vr0; //!< Zero body potential
		double vr; //!< Energy in real space
	    double vf; //!< Energy in Fourier space
		double *cx, *cy, *cz, *sx, *sy, *sz; //!< Arrays for structure factor

};

#endif

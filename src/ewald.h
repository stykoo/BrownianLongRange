/*
 * This file is an adaption of a code written by
 * Tony Maggs, PCT, Gulliver, ESPCI Paris.
 *
 * It implements the Evald method for the computation
 * of the Coulomb interaction between particles.
 */


#ifndef BROWNIANLONGRANGE_EVALD_H_
#define BROWNIANLONGRANGE_EVALD_H_

#define USE_MKL

#include <math.h>

enum Dim { DIM2 = 2, DIM3 = 3 };

class Ewald{
	public:
		Ewald(const int _N, double *_p, const double bias, const double *_q,
			  const Dim _D=DIM3);
		~Ewald();

		double* fullforce(double *);
		void dump();
		double* getForce() { return force; }

		static constexpr double error = 1e-8; //!< Error

	private:
		void initReal();
		void initFourier();
		void realSpace();
		void realSpaceAux(int i, int j, double *v, double *f);
#ifdef USE_MKL
		void realSpaceAux2MKL(int ii);
#endif
		double fourierSpace();
		void getPoten();
		void getStruct();
		void arrays(int i);

		const int D; //!< Dimension of the space in which the particles are
		const int N; //!< Number of particles
		double *p; //!< Positions of the particles
		const double *Q; //!< Particle charges
		double alpha; //!< Parameter of Evald method
		double alpha2; //!< Square of parameter of Evald method
		double rRange; //!< Range in real space
		double fRange; //!< Range in Fourier space
		int lq; //!< Integer range in Fourier space
		int size; //!< Size of arrays in 1d Fourier space
		int sizeTot; //!< Size of arrays in 2d/3d Fourier space
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

#ifdef USE_MKL
		int rs_lo, rs_hi, rs_size, rs_sizeTot;
		double *norms2, *norms, *efs, *cms;
		double *rrs;
		bool *in_range;
#endif
};

#endif

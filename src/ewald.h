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
 * \file ewald.h
 * \brief Ewald method for Coulomb interactions
 *
 * Definition of the main class for Ewald methods,
 * and of the two classes taking care of real space and Fourier space.
 */


#ifndef _EVALD_H_
#define _EVALD_H_

#define M_PI 3.14159265358979323846

//! Dimension of the system
enum Dim {
	DIM2 = 2, //!< Two dimensions
	DIM3 = 3  //!< Three dimensions
};

class EwaldRS;
class EwaldFS;

//!  Class for the Ewald method
class Ewald {
	public:
		Ewald(const Dim D_, const int N_, double *p_, const double *q_,
			  const double bias_, int use_mkl_ = 0); //!< Constructor
		~Ewald(); //!< Destructor

		double* update(double &vv); //!< Update the energy and the forces
		void dump() const; //!< Dump values to stdout
		//! Return the current total energy
		double getEnergy() const { return energy; }
		//! Return the current total forces
		double* getForces() { return forces; }

		static constexpr double error = 1e-8; //!< Error of Ewald method

	private:
		const Dim D; //!< Dimension of the space in which the particles are
		const int N; //!< Number of particles
		double *local_q; //!< Local copy of the charges
		double energy; //!< Total energy
	    double *forces; //!< Total forces

		EwaldRS *ewReal; //!< Ewald method in real space
		EwaldFS *ewFourier; //!< Ewald method in Fourier space
};

//!  Class for the real space contribution of the Ewald method
class EwaldRS {
	public:
		//! Constructor
		EwaldRS(const Dim D_, const int N_, double *p_, const double *q_, 
		        const double alpha_, const double range_, int use_mkl_ = 0);
		~EwaldRS(); //!< Destructor

		void update(); //!< Update the real space energy and forces
		//! Return the current energy in real space
		double getEnergy() const { return energy; }
		//! Return the current forces in real space
		double* getForces() { return forces; }

	private:
		void updateAux(const int ii, const int jj);
#ifdef USE_MKL
		void updateAuxMKL1(const int ii);
		void updateAuxMKL2();
		void opsMKL(const int n);
#endif

		const Dim D; //!< Dimension of the space in which the particles are
		const int N; //!< Number of particles
		double *p; //!< Positions of the particles
		const double *Q; //!< Charges of the particles
		const double alpha; //!< Parameter of Evald method
		const double range; //!< Range in real space
		int use_mkl; //!< Use MKL library

		double energy; //!< Energy in real space
		double *forces; //!< Forces in real space
		double zero_body; //!< Zero body potential

		int shift_max; //!< Number of images in real space
		int n_images_x; //!< Number of images in a given direction
		int n_images_tot; //!< Total number of images

		// Arrays needed when using MKL
#ifdef USE_MKL
		double *norms2, *norms, *efs, *cms;
		double *rrs;
		bool *in_range;
#endif
};

//!  Class for the Fourier space contribution of the Ewald method
class EwaldFS {
	public:
		EwaldFS(const Dim D_, const int N_, double *p_, const double *Q_,
		        const double alpha_, const double range_); //!< Constructor
		~EwaldFS(); //!< Destructor

		void update(); //!< Update the Fourier space energy and forces
		//! Return the current energy in Fourier space
		double getEnergy() const { return energy; }
		//! Return the current forces in Fourier space
		double* getForces() { return forces; }
		
	private:
		void computeEnergy(); //!< Compute the Fourier space energy
		//! Compute the arrays for the structure factor
		void computeStructureFactor();
		//! Compute the arrays of cosine and sin related to the positions
		void fillCosSin(int i);

		const Dim D; //!< Dimension of the space in which the particles are
		const int N; //!< Number of particles
		double *p; //!< Positions of the particles
		const double *Q; //!< Particle charges
		const double alpha; //!< Parameter of Evald method
		const double range; //!< Range in Fourier space

		int lq; //!< Integer range in Fourier space
		int size; //!< Size of arrays in 1d Fourier space
		int sizeTot; //!< Size of arrays in 2d/3d Fourier space

		double energy; //!< Energy
		double *forces; //!< Forces

		double *a; //!< Fourier space sum
		double *c; //!< Cos factors
		double *s; //!< Sin factors
		double *qx, *qy, *qz; //!< Wavevectors
		double zero_body; //!< Zero body potential in Fourier space
		double *cx, *cy, *cz, *sx, *sy, *sz; //!< Arrays for structure factor
};

#endif // _EVALD_H_

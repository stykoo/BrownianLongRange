/*
Copyright (C) Sorbonne Université (2018)
Contributor: Alexis Poncet <aponcet@lptmc.jussieu.fr>

This file is part of BrownianLongRange.

BrownianLongRange is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BrownianLongRange is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BrownianLongRange.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file state.h
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief State of the system
 *
 * Header file for state.h.
 * It defines the class State.
 */

#ifndef BROWNIANLONGRANGE_STATE_H_
#define BROWNIANLONGRANGE_STATE_H_

#include <vector>
#include <array>
#include <random>

#define DIM 3


/*!
 * \brief Class for the state of the system
 *
 * This class takes care of the initialization
 * and the evolution of the state of the system.
 */
class State {
	public:
		//! Constructor of State
		State(const double _len, const long _n_parts, const long _n_parts_1,
		      const double _charge1, const double _charge2,
		      const double _pot_strength, const double _temperature,
			  const double _dt, const double _mass);
		~State() { }
		void evolveInertia(); //!< Do one time step with inertia
		void evolveNoInertia(); //!< Do one time step without inertia

		//! Get the positions 
		const std::vector<double> & getPos() const {
			return positions;
		}


	private:
		const double len; //!< Length of the box
		const long n_parts; //!< Number of particles
		const long n_parts_1; //!< Number of particles of species 1
		double charge1; //!< Charge of species 1
		double charge2; //!< Charge of species 1
		const double pot_strength; //!< Strength of the interparticle potential
		const double dt; //!< Timestep
		double mass; //!< Mass of the particles if inertial dynamics

		std::mt19937 rng; //!< Random number generator
		//! Gaussian noise for temperature
		std::normal_distribution<double> noise;

		std::vector<double> positions; //!< Positions of the particles
		std::vector<double> forces;  //!< Internal forces
};

/*! 
 * \brief Periodic boundary conditions on a segment
 * 
 * Update x to be between 0 and L.
 *
 * \param x Value
 * \param L Length of the box
 */
template<typename T>
void pbc(T &x, const T L){
	x -= L * std::floor(x / L);
}

/*! 
 * \brief Periodic boundary conditions on a segment (symmetric)
 * 
 * Update x to be between -L/2 and L/2.
 *
 * \param x Value
 * \param L Length of the box
 */
template<typename T>
inline void pbcSym(T &x, const T L) {
	x -= L * std::round(x / L);
}

// Trick to avoid round ASSUMING LITTLE ENDIAN
union i_cast {double d; int i[2];};
#define double2int(i, d, t)  \
    {volatile union i_cast u; u.d = (d) + 6755399441055744.0; \
    (i) = (t)u.i[0];}

template<>
inline void pbcSym<double>(double &x, const double L) {
	double d = x / L;
	int i;
	double2int(i, d, int);
	x -= L * i;
}

#endif // BROWNIANLONGRANGE_STATE_H_

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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "ewald.h"

/*!
 * \brief Class for the state of the system
 *
 * This class takes care of the initialization
 * and the evolution of the state of the system.
 */
class State {
	public:
		//! Constructor of State
		State(const long _n_parts, const long _n_parts_1,
		      const double _charge1, const double _charge2,
		      const double _pot_strength, const double _temperature,
			  const double _field, const double _dt, const double _mass,
			  const double _bias, const int _D);
		~State();
		void evolveNoInertia(); //!< Do one time step without inertia
		void evolveInertia(); //!< Do one time step with inertia
		void resetPosIni(); //!< Reset the initial positions
		void getDisplacements(double &X1, double &X2) const;

		//! Get the positions 
		const std::vector<double> & getPos() const {
			return positions;
		}


	private:
		const int D; //!< Dimension (2 or 3)
		const long n_parts; //!< Number of particles
		const long n_parts_1; //!< Number of particles of species 1
		const double pot_strength; //!< Strength of the interparticle potential
		const double field; //!< External field
		const double dt; //!< Timestep
		double mass; //!< Mass of the particles if inertial dynamics

		double sigma;
		gsl_rng *rng;

		std::vector<double> positions; //!< Positions of the particles
		std::vector<double> positions_ini; //!< Initial positions
		std::vector<double> offsets; //!< Number of turns around the box
		std::vector<double> velocities; //!< Velocities of the particles
		std::vector<double> charges; //!< Charges of the particles
		Ewald *ew; //!< To perform Ewald summation
};

#endif // BROWNIANLONGRANGE_STATE_H_

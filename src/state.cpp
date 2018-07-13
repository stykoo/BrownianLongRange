/*
Copyright (C) Sorbonne Université (2018)
Contributor: Alexis Poncet <aponcet@lptmc.jussieu.fr>

This file is part of ActiveBrownian.

ActiveBrownian is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ActiveBrownian is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ActiveBrownian.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file simul.cpp
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief State of the system
 *
 * Implementation of the methods of the class State to simulate
 * interacting active Brownian particles in dimension 2.
*/

#include <cmath>
#include <chrono>
#include <algorithm>
#include "state.h"

/*!
 * \brief Constructor of State
 *
 * Initializes the state of the system: particles randomly placed in a 2d box.
 *
 * \param _len Length of the box
 * \param _n_parts Number of particles
 * \param _pot_strength Strength of the interparticle potential
 * \param _temperature Temperature
 * \param _rot_dif Rotational diffusivity
 * \param _activity Activity
 * \param _dt Timestep
 * \param _fac_boxes Factor for the boxes
 */
State::State(const double _len, const long _n_parts, const long _n_parts_1,
	         const double _charge1, const double _charge2,
	         const double _pot_strength, const double _temperature,
	         const double _dt, const double _mass) :
	len(_len), n_parts(_n_parts), n_parts_1(_n_parts_1),
	charge1(_charge1), charge2(_charge2), pot_strength(_pot_strength),
	dt(_dt), mass(_mass),
	// We seed the RNG with the current time
	rng(std::chrono::system_clock::now().time_since_epoch().count()),
	// Gaussian noise from the temperature
	noise(0.0, std::sqrt(2.0 * _temperature * dt))
{
	positions.resize(DIM * n_parts);
	forces.assign(DIM * n_parts, 0);

    std::uniform_real_distribution<double> rndPos(0, len);

	for (long i = 0 ; i < DIM * n_parts ; ++i) {
		positions[i] = rndPos(rng);
		forces[i] = 0;
	}
}

/*!
 * \brief Do one time step
 *
 * Evolve the system for one time step according to coupled Langevin equation.
 */
void State::evolveNoInertia() {
	//calcInternalForces();

	for (long i = 0 ; i < DIM * n_parts ; ++i) {
		// Internal forces + Gaussian noise
		positions[i] += dt * forces[i] + noise(rng);
	}

	// Enforce PBC
	for (long i = 0 ; i < DIM * n_parts ; ++i) {
		pbc(positions[i], len);
	}
}

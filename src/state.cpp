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
#include <random>
#include <algorithm>
#include "state.h"

/*!
 * \brief Constructor of State
 *
 * Initializes the state of the system: particles randomly placed in a 2d box.
 */
State::State(const long _n_parts, const long _n_parts_1,
	         const double _charge1, const double _charge2,
	         const double _pot_strength, const double _temperature,
	         const double _dt, const double _mass, const double _bias=1.0) :
	n_parts(_n_parts), n_parts_1(_n_parts_1), pot_strength(_pot_strength),
	dt(_dt), mass(_mass), sigma(std::sqrt(2.0 * _temperature * dt))
{
	std::random_device rn;
	unsigned int sd = rn() ;
	rng = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(rng, sd);

	positions.resize(DIM * n_parts);
	for (long i = 0 ; i < DIM * n_parts ; ++i) {
	//for (long i = 0 ; i < 2 * n_parts ; ++i) {
		positions[i] = gsl_ran_flat(rng, 0.0, 1.0);
	}
	/*for (long i = 2 * n_parts ; i < DIM * n_parts ; ++i) {
		positions[i] = 0;
	}*/

	charges.resize(n_parts);
	for (long i = 0 ; i < n_parts_1 ; ++i) {
		charges[i] = _charge1;
	}
	for (long i = n_parts_1 ; i < n_parts ; ++i) {
		charges[i] = _charge2;
	}

	ew = new Ewald(n_parts, positions.data(), _bias, charges.data());
}

State::~State() {
	gsl_rng_free(rng);
	delete ew;
}

/*!
 * \brief Do one time step
 *
 * Evolve the system for one time step according to coupled Langevin equation.
 */
void State::evolveNoInertia() {
	// Compute internal forces
	double pot;
	double *forces = ew->fullforce(&pot);

	for (long i = 0 ; i < DIM * n_parts ; ++i) {
		// Internal forces + Gaussian noise
		positions[i] += dt * forces[i];
		positions[i] += gsl_ran_gaussian_ziggurat(rng, sigma);
	}

	// Enforce PBC
	for (long i = 0 ; i < DIM * n_parts ; ++i) {
		positions[i] -= std::floor(positions[i]);
	}

	/*for (long i = 2 * n_parts ; i < DIM * n_parts ; ++i) {
		positions[i] = 0;
	}*/
}

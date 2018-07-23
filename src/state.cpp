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
#ifdef RESTRICT_2D
	for (long i = 0 ; i < 2 * n_parts ; ++i) {
		positions[i] = gsl_ran_flat(rng, 0.0, 1.0);
	}
	for (long i = 2 * n_parts ; i < DIM * n_parts ; ++i) {
		positions[i] = 0;
	}
#else
	for (long i = 0 ; i < DIM * n_parts ; ++i) {
		positions[i] = gsl_ran_flat(rng, 0.0, 1.0);
	}
#endif

	velocities.assign(DIM * n_parts, 0.0);

	charges.resize(n_parts);
	for (long i = 0 ; i < n_parts_1 ; ++i) {
		charges[i] = _charge1;
	}
	for (long i = n_parts_1 ; i < n_parts ; ++i) {
		charges[i] = _charge2;
	}

	ew = new Ewald(n_parts, positions.data(), _bias, charges.data());
	// Initialize the forces (needed for inertial dynamics)
	double pot;
	ew->fullforce(&pot);
}

State::~State() {
	gsl_rng_free(rng);
	delete ew;
}

/*!
 * \brief Do one time step
 *
 * Evolve the system for one time step according to coupled
 * overdamped Langevin equation.
 */
void State::evolveNoInertia() {
	// Compute internal forces
	double pot;
	double *forces = ew->fullforce(&pot);

#ifdef RESTRICT_2D
	long hi = 2 * n_parts;
#else
	long hi = DIM * n_parts;
#endif

	for (long i = 0 ; i < hi ; ++i) {
		// Internal forces + External force + Gaussian noise
		positions[i] += dt * forces[i];
		positions[i] += gsl_ran_gaussian_ziggurat(rng, sigma);
		// Enforce PBC
		positions[i] -= std::floor(positions[i]);
	}
}

/*!
 * \brief Do one time step
 *
 * Evolve the system for one time step according to coupled Langevin equation
 * with intertia.
 *
 * See doi/10.1080/00268976.2012.760055 for algorithm (Eq. 21-22)
 */
void State::evolveInertia() {
	double b = 1.0 / (1.0 + dt / (2.0 * mass));
	double a = b * (1.0 - dt / (2.0 * mass));
	double c_rv = b * dt;
	double c_rf = b * dt * dt / (2.0 * mass);
	double c_rn = b * dt / (2.0 * mass);
	double c_vf = dt / (2.0 * mass);
	double c_vn = b / mass;

#ifdef RESTRICT_2D
	long hi = 2 * n_parts;
#else
	long hi = DIM * n_parts;
#endif

	double *forces = ew->getForce();

	for (long i = 0 ; i < hi ; ++i) {
		double u = gsl_ran_gaussian_ziggurat(rng, sigma);

		// Update positions
		positions[i] += c_rv * velocities[i];
		positions[i] += c_rf * forces[i];
		positions[i] += c_rn * u;
		positions[i] -= std::floor(positions[i]); // PBC

		// Update velocities (old forces)
		velocities[i] *= a;
		velocities[i] += c_vf * a * forces[i];
		velocities[i] += c_vn * u;
	}

	// New forces
	double pot;
	forces = ew->fullforce(&pot);

	for (long i = 0 ; i < hi ; ++i) {
		// Update velocities (new forces)
		velocities[i] += c_vf * forces[i];
	}
}

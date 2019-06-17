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

#include <iostream>

/*!
 * \brief Constructor of State
 *
 * Initializes the state of the system: particles randomly placed in a 2d box.
 */
State::State(const long _n_parts, const long _n_parts_1,
	         const double _charge1, const double _charge2,
	         const double _pot_strength, const double _temperature,
	         const double _field, const double _dt, const double _mobility1,
		 const double _mobility2, const double _mass1, const double _mass2,
		 const double _bias=1.0, const int _D=3) :
	D(_D), n_parts(_n_parts), n_parts_1(_n_parts_1), pot_strength(_pot_strength),
	field(_field), dt(_dt)
{
	std::random_device rn;
	unsigned int sd = rn() ;
	rng = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(rng, sd);

	positions.resize(D * n_parts);
	positions_ini.resize(D * n_parts);
	offsets.resize(D * n_parts);

	for (long i = 0 ; i < D * n_parts ; ++i) {
		positions[i] = gsl_ran_flat(rng, 0.0, 1.0);
	}
	resetPosIni();

	velocities.assign(D * n_parts, 0.0);
	
	charges.resize(D * n_parts);
	mobilities.resize(D * n_parts);
	masses.resize(D * n_parts);
	sigma.resize(D * n_parts);
	for (long j = 0 ; j < D ; ++j) {
		for (long i = 0 ; i < n_parts_1 ; ++i) {
			charges[j * n_parts + i] = _charge1;
			mobilities[j * n_parts + i] = _mobility1;
			masses[j * n_parts + i] = _mass1;
			sigma[j * n_parts + i] = std::sqrt(2 * _temperature * _mobility1 * dt);
		}
		for (long i = n_parts_1 ; i < n_parts ; ++i) {
			charges[j * n_parts + i] = _charge2;
			mobilities[j * n_parts + i] = _mobility2;
			masses[j * n_parts + i] = _mass2;
			sigma[j * n_parts + i] = std::sqrt(2 * _temperature * _mobility2 * dt);
		}
	}
	
	Dim DD = (D == 3) ? DIM3 : DIM2;

#ifdef USE_MKL
	ew = new Ewald(DD, n_parts, positions.data(), charges.data(), _bias, 1);
#else
	ew = new Ewald(DD, n_parts, positions.data(), charges.data(), _bias);
#endif

	// Initialize the forces (needed for inertial dynamics)
	double pot;
	ew->update(pot);
}

State::~State() {
	gsl_rng_free(rng);
	delete ew;
}

/*!
 * \brief Do one time step without inertia
 *
 * Evolve the system for one time step according to coupled
 * overdamped Langevin equation.
 */
void State::evolveNoInertia() {
	// Compute internal forces
	double pot;
	double *forces = ew->update(pot);

	for (long i = 0 ; i < n_parts ; ++i) {
		// External forces
		positions[i] += mobilities[i] * dt * field * charges[i];
	}
	for (long i = 0 ; i < D * n_parts ; ++i) {
		// Internal forces + Gaussian noise
		positions[i] += mobilities[i] * dt * pot_strength * forces[i];
		positions[i] += gsl_ran_gaussian_ziggurat(rng, sigma[i]);
		// Enforce PBC (periodic boundary conditions)
		double o = std::floor(positions[i]);
		positions[i] -= o;
		offsets[i] += o;
	}
}

/*!
 * \brief Do one time step with inertia
 *
 * Evolve the system for one time step according to coupled Langevin equations
 * with intertia.
 *
 * See doi/10.1080/00268976.2012.760055 for algorithm (Eq. 21-22)
 */
void State::evolveInertia() {
	// Old forces
	double *forces = ew->getForces();

	for (long i = 0 ; i < D * n_parts ; ++i) {
		double b = 1.0 / (1.0 + dt / (2.0 * mobilities[i] * masses[i]));
		double a = b * (1.0 - dt / (2.0 * mobilities[i] * masses[i]));
		double c_rv = b * dt;
		double c_rf = b * dt * dt / (2.0 * masses[i]);
		double c_rn = b * dt / (2.0 * masses[i]);
		double c_vf = dt / (2.0 * masses[i]);
		double c_vn = b / masses[i];

		if (i / n_parts_1 == 0) {
			// External forces
			positions[i] += c_rf * field * charges[i];
		}

		double u = gsl_ran_gaussian_ziggurat(rng, sigma[i]);
		
		// Update positions
		positions[i] += c_rv * velocities[i];
		positions[i] += c_rf * pot_strength * forces[i];
		positions[i] += c_rn * u;
		// Enforce PBC
		double o = std::floor(positions[i]);
		positions[i] -= o;
		offsets[i] += o;

		// Update velocities (old internal forces)
		velocities[i] *= a;
		velocities[i] += c_vf * a * pot_strength * forces[i];
		velocities[i] += c_vn * u;

		if (i / n_parts_1 == 0) {
			// External forces (old + new)
			velocities[i] += c_vf * (a + 1.0) * field * charges[i];
		}
	}

	// New forces
	double pot;
	forces = ew->update(pot);

	for (long i = 0 ; i < D * n_parts ; ++i) {
		double c_vf = dt / (2.0 * masses[i]);

		// Update velocities (new internal forces)
		velocities[i] += c_vf * pot_strength * forces[i];
	}
}

//! Reset the initial positions
void State::resetPosIni() {
	for (long i = 0 ; i < D * n_parts ; ++i) {
		positions_ini[i] = positions[i];
		offsets[i] = 0.0;
	}
}

//! Get the average displacement of particles 1 (resp. 2) along 1st axis
void State::getDisplacements(double &X1, double &X2) const {
	X1 = 0.0;
	for (long i = 0 ; i < n_parts_1 ; ++i) {
		X1 += (positions[i] - positions_ini[i] + offsets[i]);
	}
	X1 /= n_parts_1;

	X2 = 0.0;
	for (long i = n_parts_1 ; i < n_parts ; ++i) {
		X2 += (positions[i] - positions_ini[i] + offsets[i]);
	}
	X2 /= (n_parts - n_parts_1);
}

//! Get the average internal force (ie not the electric field) applied on particles 1 (resp. 2) along 1st axis
void State::getInternalForces(double &F1, double &F2) const {
	double *forces = ew->getForces();
	
	F1 = 0.0;
	for (long i = 0 ; i < n_parts_1 ; ++i) {
		F1 += (pot_strength * forces[i]);
	}
	F1 /= n_parts_1;
	
	F2 = 0.0;
	for (long i = n_parts_1 ; i < n_parts ; ++i) {
		F2 += (pot_strength * forces[i]);
	}
	F2 /= (n_parts - n_parts_1);
}


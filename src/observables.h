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
 * \file observables.h
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Computation and export of the observables
 *
 * Header file for observables.cpp.
*/

#ifndef BROWNIANLONGRANGE_OBSERVABLES_H
#define BROWNIANLONGRANGE_OBSERVABLES_H

#include <vector>
#include "state.h"

class Observables {
	public:
		Observables(const long n_parts_, const long n_parts_1_,
				    const long n_div_x_, const long n_pts_, const int D_);
		//! Compute the observables for a given state
		void compute(const State *state, const long t);
		//! Export to hdf5
		void writeH5(const std::string fname, double charge1, double charge2,
	                 double pot_strength, double temperature, double field,
					 double dt, long n_iters, long n_iters_th, double bias,
					 long skip, int inertia) const;

	private:
		size_t boxOfPair(const long i, const long j,
				         const std::vector<double> & pos);

		const int D; //!< Dimension (2 or 3)
		const long n_parts; //!< Number of particles 
		const long n_parts_1; //!< Number of particles of species 1
		const long n_div_x; //!< Number of divisions in x
		const long n_div_r; //!< Number of divisions for radial coordinate
		const long n_div_tot; //!< Total number of divisions
		const long n_pts; //!< Number of points in time

		std::vector<double> displ1; //!< Displacement of particles 1
		std::vector<double> displ2; //!< Displacement of particles 1
		std::vector<long long> correls11; //!< Correlations
		std::vector<long long> correls22; //!< Correlations
		std::vector<long long> correls12; //!< Correlations
};

#endif // BROWNIANLONGRANGE_OBSERVABLES_H

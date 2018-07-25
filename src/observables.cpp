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
 * \file observables.cpp
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Computation and export of the observables
*/

#include <iostream>
//#include <cassert>
#include <cmath>
#include "H5Cpp.h"
#include "observables.h"

/*
 * \brief Constructor of Observables
 *
 * Initialize the vector for correlations.
 */
Observables::Observables(const long n_parts_, const long n_parts_1_,
			             const long n_div_x_) :
		n_parts(n_parts_), n_parts_1(n_parts_1_),
		n_div_x(n_div_x_),
#ifdef RESTRICT_2D
		n_div_r((n_div_x_ + 1) / 2),
#else
		n_div_r((long) std::ceil(std::sqrt(0.5) * n_div_x)),
#endif
		n_div_tot(n_div_x * n_div_r)
		//n_div_tot(n_div_x * n_div_x * n_div_x)
{
	correls11.assign(n_div_tot, 0);
	correls22.assign(n_div_tot, 0);
	correls12.assign(n_div_tot, 0);
}

/*
 * \brief Compute the observables for a given state
 */
void Observables::compute(const State *state) {
	const std::vector<double> & pos = state->getPos();

	// For each pair of particles
	for (long i = 0 ; i < n_parts_1 ; ++i) {
		for (long j = i + 1 ; j < n_parts_1 ; ++j) {
			correls11[boxOfPair(i, j, pos)]++; // Add 1 in the right box
		}
		for (long j = n_parts_1 ; j < n_parts ; ++j) {
			correls12[boxOfPair(i, j, pos)]++; // Add 1 in the right box
		}
	}

	for (long i = n_parts_1 ; i < n_parts ; ++i) {
		for (long j = i + 1 ; j < n_parts ; ++j) {
			correls22[boxOfPair(i, j, pos)]++; // Add 1 in the right box
		}
	}
}

size_t Observables::boxOfPair(const long i, const long j,
		                      const std::vector<double> & pos) {
#ifdef RESTRICT_2D
	double dr[2];
	for (int a = 0 ; a < 2 ; ++a) {
		dr[a] = pos[a * n_parts + j] - pos[a * n_parts + i];
		dr[a] -= std::round(dr[a]);
	}
	size_t bx = (size_t) (n_div_x * (dr[0] + 0.5));
	size_t br = (size_t) (n_div_x * std::abs(dr[1]));
#else
	//size_t b[3];
	double dr[3];
	for (int a = 0 ; a < 3 ; ++a) {
		dr[a] = pos[a * n_parts + j] - pos[a * n_parts + i];
		dr[a] -= std::round(dr[a]);
		//dr[a] -= std::floor(dr[a]);
		//b[a] = (size_t) (dr[a] * n_div_x);
	}
	size_t bx = (size_t) (n_div_x * (dr[0] + 0.5));
	size_t br = (size_t) (n_div_x * std::sqrt(dr[1]*dr[1] + dr[2]*dr[2]));
#endif
	return bx * n_div_r + br;
	//return b[0] * n_div_x * n_div_x + b[1] * n_div_x + b[2];	
}

/*
 * \brief Export the observables to a hdf5 file
 */
void Observables::writeH5(const std::string fname, double charge1,
		                  double charge2, double pot_strength,
						  double temperature, double field, double dt,
						  long n_iters, long n_iters_th, double bias, long skip,
						  int inertia) const {
	try {
		H5::H5File file(fname, H5F_ACC_TRUNC);

		// General attributes
		H5::DataSpace default_ds;
		H5::Attribute a_n_parts = file.createAttribute(
				"n_parts", H5::PredType::NATIVE_LONG, default_ds);
		a_n_parts.write(H5::PredType::NATIVE_LONG, &n_parts);
		H5::Attribute a_n_parts_1 = file.createAttribute(
				"n_parts_1", H5::PredType::NATIVE_LONG, default_ds);
		a_n_parts_1.write(H5::PredType::NATIVE_LONG, &n_parts_1);
		H5::Attribute a_charge1 = file.createAttribute(
				"charge1", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_charge1.write(H5::PredType::NATIVE_DOUBLE, &charge1);
		H5::Attribute a_charge2 = file.createAttribute(
				"charge2", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_charge2.write(H5::PredType::NATIVE_DOUBLE, &charge2);
		H5::Attribute a_pot_strength = file.createAttribute(
				"pot_strength", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_pot_strength.write(H5::PredType::NATIVE_DOUBLE, &pot_strength);
		H5::Attribute a_temperature = file.createAttribute(
				"temperature", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_temperature.write(H5::PredType::NATIVE_DOUBLE, &temperature);
		H5::Attribute a_field = file.createAttribute(
				"field", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_field.write(H5::PredType::NATIVE_DOUBLE, &field);
		H5::Attribute a_dt = file.createAttribute(
				"dt", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_dt.write(H5::PredType::NATIVE_DOUBLE, &dt);
		H5::Attribute a_n_iters = file.createAttribute(
				"n_iters", H5::PredType::NATIVE_LONG, default_ds);
		a_n_iters.write(H5::PredType::NATIVE_LONG, &n_iters);
		H5::Attribute a_n_iters_th = file.createAttribute(
				"n_iters_th", H5::PredType::NATIVE_LONG, default_ds);
		a_n_iters_th.write(H5::PredType::NATIVE_LONG, &n_iters_th);
		H5::Attribute a_bias = file.createAttribute(
				"bias", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_bias.write(H5::PredType::NATIVE_DOUBLE, &bias);
		H5::Attribute a_skip = file.createAttribute(
				"skip", H5::PredType::NATIVE_LONG, default_ds);
		a_skip.write(H5::PredType::NATIVE_LONG, &skip);
		H5::Attribute a_inertia = file.createAttribute(
				"inertia", H5::PredType::NATIVE_INT, default_ds);
		a_inertia.write(H5::PredType::NATIVE_INT, &inertia);
		
		// We chunk the data and compress it
		// Chunking should depend on how we intend to read the data
		H5::DataSet dataset11, dataset22, dataset12;
		H5::DSetCreatPropList plist;
		plist.setDeflate(6);
		/*long chunk_w = std::min(100l, n_div_x);
		hsize_t chunk_dims[3] = {1, (hsize_t) chunk_w, (hsize_t) chunk_w};
		plist.setChunk(3, chunk_dims);*/
		long chunk_w = std::min(1024l, n_div_r);
		hsize_t chunk_dims[2] = {1, (hsize_t) chunk_w};
		plist.setChunk(2, chunk_dims);

		// Dimensions of the data
		/*hsize_t dims[3] = {(hsize_t) n_div_x, (hsize_t) n_div_x,
			               (hsize_t) n_div_x};
		H5::DataSpace dataspace(3, dims);*/
		hsize_t dims[3] = {(hsize_t) n_div_x, (hsize_t) n_div_r};
		H5::DataSpace dataspace(2, dims);

		// Write datasets for correlations
		dataset11 = file.createDataSet("correlations11",
				   					   H5::PredType::NATIVE_LLONG,
									   dataspace, plist);
		dataset11.write(correls11.data(), H5::PredType::NATIVE_LLONG);
		dataset22 = file.createDataSet("correlations22",
				   					   H5::PredType::NATIVE_LLONG,
									   dataspace, plist);
		dataset22.write(correls22.data(), H5::PredType::NATIVE_LLONG);
		dataset12 = file.createDataSet("correlations12",
				   					   H5::PredType::NATIVE_LLONG,
									   dataspace, plist);
		dataset12.write(correls12.data(), H5::PredType::NATIVE_LLONG);

		// Attributes for correlations
		H5::Attribute a_n_div_x11 = dataset11.createAttribute(
				"n_div_x", H5::PredType::NATIVE_LONG, default_ds);
		a_n_div_x11.write(H5::PredType::NATIVE_LONG, &n_div_x);
		H5::Attribute a_n_div_x22 = dataset22.createAttribute(
				"n_div_x", H5::PredType::NATIVE_LONG, default_ds);
		a_n_div_x22.write(H5::PredType::NATIVE_LONG, &n_div_x);
		H5::Attribute a_n_div_x12 = dataset12.createAttribute(
				"n_div_x", H5::PredType::NATIVE_LONG, default_ds);
		a_n_div_x12.write(H5::PredType::NATIVE_LONG, &n_div_x);
	} catch (H5::Exception& err) {
        err.printError();
	}
}

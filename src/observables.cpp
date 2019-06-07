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
			             const long n_div_x_, const long n_pts_, const int D_) :
		D(D_), n_parts(n_parts_), n_parts_1(n_parts_1_), n_div_x(n_div_x_),
		n_div_r((D == 3)
				? (long) std::ceil(std::sqrt(0.5) * n_div_x)
				: (n_div_x_ + 1) / 2),
		n_div_tot(n_div_x * n_div_r), n_pts(n_pts_)
{
	displ1.assign(n_pts, 0.0);
	displ2.assign(n_pts, 0.0);
	intForces1.assign(n_pts, 0.0);
	intForces2.assign(n_pts, 0.0);
	correls12.assign(n_div_tot, 0);
	correls11.assign(n_div_tot, 0);
	correls22.assign(n_div_tot, 0);
	correls12.assign(n_div_tot, 0);
}

/*
 * \brief Compute the observables for a given state at 
 */
void Observables::compute(const State *state, const long t) {
	double X1, X2;
	state->getDisplacements(X1, X2);
	displ1[t] = X1;
	displ2[t] = X2;
	
	double F1, F2;
	state->getInternalForces(F1, F2);
	intForces1[t] = F1;
	intForces2[t] = F2;
	
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
	size_t bx = 0, br = 0;
	if (D == 2) {
		double dr[2];
		for (int a = 0 ; a < 2 ; ++a) {
			dr[a] = pos[a * n_parts + j] - pos[a * n_parts + i];
			dr[a] -= std::round(dr[a]);
		}
		bx = (size_t) (n_div_x * (dr[0] + 0.5));
		br = (size_t) (n_div_x * std::abs(dr[1]));
	} else if (D == 3) {
		double dr[3];
		for (int a = 0 ; a < 3 ; ++a) {
			dr[a] = pos[a * n_parts + j] - pos[a * n_parts + i];
			dr[a] -= std::round(dr[a]);
		}
		bx = (size_t) (n_div_x * (dr[0] + 0.5));
		br = (size_t) (n_div_x * std::sqrt(dr[1]*dr[1] + dr[2]*dr[2]));
	}
	return bx * n_div_r + br;
}

/*
 * \brief Export the observables to a hdf5 file
 */
void Observables::writeH5(const std::string fname, double charge1, double charge2,
				double mobility1, double mobility2, double pot_strength,
				double temperature, double field, double dt, long n_iters,
				long n_iters_th, double bias, long skip, int inertia,
				double mass1, double mass2) const {
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
		H5::Attribute a_mobility1 = file.createAttribute(
				"mobility1", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_mobility1.write(H5::PredType::NATIVE_DOUBLE, &mobility1);
		H5::Attribute a_mobility2 = file.createAttribute(
				"mobility2", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_mobility2.write(H5::PredType::NATIVE_DOUBLE, &mobility2);
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
		H5::Attribute a_mass1 = file.createAttribute(
				"mass1", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_mass1.write(H5::PredType::NATIVE_DOUBLE, &mass1);
		H5::Attribute a_mass2 = file.createAttribute(
				"mass2", H5::PredType::NATIVE_DOUBLE, default_ds);
		a_mass2.write(H5::PredType::NATIVE_DOUBLE, &mass2);
		H5::Attribute a_dimension = file.createAttribute(
				"dimension", H5::PredType::NATIVE_INT, default_ds);
		a_dimension.write(H5::PredType::NATIVE_INT, &D);

		/* DISPLACEMENTS */
		std::vector<double> data_displ(3 * n_pts);
		for (long i = 0 ; i < n_pts ; ++i) {
			data_displ[3 * i] = i * skip * dt;
			data_displ[3 * i + 1] = displ1[i];
			data_displ[3 * i + 2] = displ2[i];
		}

		H5::DataSet datasetDispl;
		H5::DSetCreatPropList plistDispl;
		plistDispl.setDeflate(6);
		hsize_t dimsDispl[2] = {(hsize_t) n_pts, 3};
		H5::DataSpace dataspaceDispl(2, dimsDispl);
		hsize_t chunk_displ[2] = {(hsize_t) std::min(1024l, n_pts), 1};
		plistDispl.setChunk(2, chunk_displ);

		datasetDispl = file.createDataSet("displacements",
				H5::PredType::NATIVE_DOUBLE, dataspaceDispl, plistDispl);
		datasetDispl.write(data_displ.data(), H5::PredType::NATIVE_DOUBLE);
		
		/* INTERNAL FORCES */
		std::vector<double> data_intForces(3 * n_pts);
		for (long i = 0 ; i < n_pts ; ++i) {
			data_intForces[3 * i] = i * skip * dt;
			data_intForces[3 * i + 1] = intForces1[i];
			data_intForces[3 * i + 2] = intForces2[i];
		}

		H5::DataSet datasetIntForces;
		H5::DSetCreatPropList plistIntForces;
		plistIntForces.setDeflate(6);
		hsize_t dimsIntForces[2] = {(hsize_t) n_pts, 3};
		H5::DataSpace dataspaceIntForces(2, dimsIntForces);
		hsize_t chunk_intForces[2] = {(hsize_t) std::min(1024l, n_pts), 1};
		plistIntForces.setChunk(2, chunk_intForces);

		datasetIntForces = file.createDataSet("internal_forces",
				H5::PredType::NATIVE_DOUBLE, dataspaceIntForces, plistIntForces);
		datasetIntForces.write(data_intForces.data(), H5::PredType::NATIVE_DOUBLE);
		
		/* CORRELATIONS */
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
		hsize_t dims[2] = {(hsize_t) n_div_x, (hsize_t) n_div_r};
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
        err.printErrorStack();
	}
}

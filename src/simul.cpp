/*
Copyright (C) Sorbonne Université (2018)
Contributor: Alexis Poncet <aponcet@lptmc.jussieu.fr>

This file is part of BrownianLongRange.

BrownianLongRange is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ActiveBrownian is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file simul.cpp
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Simulation of the system
 *
 * Implementation of the methods of the class Simul to simulate
 * Brownian particles interacting with electrostatic interactions.
*/

#include <exception>
#include <thread>
#include <boost/program_options.hpp>
#include "observables.h"
#include "simul.h"
#include "state.h"

#ifdef VISU
#ifdef RESTRICT_2D
#include "visu2d.h"
//#include "visu3d.h"
#endif
#endif

namespace po = boost::program_options;

/*!
 * \brief Constructor of Simul
 *
 * Initializes the parameters of structure Simul
 * from the command-line arguments using boost::program_options.
 *
 * \param argc Number of arguments
 * \param argv Arguments
 */
Simul::Simul(int argc, char **argv) {
	status = SIMUL_INIT_SUCCESS;

	po::options_description opts("Options");
	opts.add_options()
		("parts,n", po::value<long>(&n_parts)->required(),
		 "Number of particles")
		("parts1,m", po::value<long>(&n_parts_1)->required(),
		 "Number of particles of species 1")
		("q1", po::value<double>(&charge1)->required(),
		 "Charge of particles of species 1")
		("q2", po::value<double>(&charge2)->required(),
		 "Charge of particles of species 2")
		("eps,e", po::value<double>(&pot_strength)->default_value(1.0),
		 "Strength of interparticle potential")
		("temp,T", po::value<double>(&temperature)->required(), "Temperature")
		("dt,t", po::value<double>(&dt)->required(), "Timestep")
		("iters,I", po::value<long>(&n_iters)->required(),
		 "Number of time iterations")
		("itersTh,J", po::value<long>(&n_iters_th)->default_value(0),
		 "Number of time iterations of thermalization")
		("inertia,i", po::bool_switch(&inertia),
		 "Simulation with inertia")
		("mass,M", po::value<double>(&mass)->default_value(0.0),
		 "Mass of the particles")
		("bias", po::value<double>(&bias)->default_value(1.0),
		 "Bias toward Fourier space")
		("skip,S", po::value<long>(&skip)->default_value(100),
		 "Iterations between two computations of observables")
		("output,O",
		 po::value<std::string>(&output)->default_value("observables.h5"),
		 "Name of the output file")
		("div,d",
		 po::value<long>(&n_div_x)->default_value(100),
		 "Number of divisions for correlations")
		("sleep", po::value<int>(&sleep)->default_value(0),
		 "Number of milliseconds to sleep for between iterations")
		("help,h", "Print help message and exit")
		;

	try {
		po::variables_map vars;
		po::store(po::parse_command_line(argc, argv, opts), vars);

		// Display help and exit
		if (vars.count("help")) {
			std::cout << "Usage: " << argv[0] << " options\n";
			std::cout << opts << std::endl;
			status = SIMUL_INIT_HELP;
			return;
		}

        po::notify(vars);
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		status = SIMUL_INIT_FAILED;
		return;
	}

	// Check if the values of the parameters are allowed
	if (notStrPositive(n_parts, "n_parts")
		|| notPositive(n_parts_1, "n_parts_1")
		|| notPositive(pot_strength, "eps") || notPositive(temperature, "T")
		|| notStrPositive(dt, "dt") || notPositive(n_iters, "n_iters")
		|| notPositive(n_iters_th, "n_iters_th")
		|| notPositive(mass, "mass") || notStrPositive(bias, "bias")
		|| notStrPositive(n_div_x, "n_div_x")) {
		status = SIMUL_INIT_FAILED;
		return;
	}
}

/*!
 * \brief Run the simulation
 *
 * Construct the state of the system and update it for the number
 * of iterations wanted. Also take care of launching the thread for
 * visualization.
 */
void Simul::run() {
	if (status != SIMUL_INIT_SUCCESS) {
		std::cerr << "You should not be runing a failed simulation..."
		          << std::endl;
		return;
	}

	// Initialize the state of the system
	State state(n_parts, n_parts_1, charge1, charge2, pot_strength,
	            temperature, dt, mass, bias);
	Observables obs(n_parts, n_parts_1, n_div_x);

#ifdef VISU
#ifdef RESTRICT_2D
	Visu visu(&state, n_parts, n_parts_1);
	std::thread thVisu(&Visu::run, &visu); 
	//Visu3d visu(&state, n_parts, n_parts_1);
	//std::thread thVisu(&Visu3d::run, &visu); 
#endif
#endif

	if (inertia) {
		// Thermalization
		for (long t = 0 ; t < n_iters_th ; ++t) {
			state.evolveInertia();
		}
		// Time evolution
		for (long t = 0 ; t < n_iters ; ++t) {
			state.evolveInertia();
			if (t % skip == 0) {
				obs.compute(&state);
			}
			if (sleep > 0) {
				std::this_thread::sleep_for(std::chrono::milliseconds(sleep));
			}
		}
	} else {
		// Thermalization
		for (long t = 0 ; t < n_iters_th ; ++t) {
			state.evolveNoInertia();
		}
		// Time evolution
		for (long t = 0 ; t < n_iters ; ++t) {
			state.evolveNoInertia();
			if (t % skip == 0) {
				obs.compute(&state);
			}
			if (sleep > 0) {
				std::this_thread::sleep_for(std::chrono::milliseconds(sleep));
			}
		}
	}

	obs.writeH5(output, charge1, charge2, pot_strength, temperature,
				dt, n_iters, n_iters_th, bias, skip, (int) inertia);

#ifdef VISU
#ifdef RESTRICT_2D
	thVisu.join();
#endif
#endif
}

/*!
 * \brief Print the parameters of the simulation
 */
void Simul::print() const {
	if (inertia) {
		std::cout << "# [WithInertia] ";
	} else {
		std::cout << "# [WithoutInertia] ";
	}
	std::cout << "n_parts=" << n_parts
	          << ", n_parts_1=" << n_parts_1 << ", charge1="
			  << charge1 << ", charge2=" << charge2
	          << ", pot_strength=" << pot_strength << ", temperature="
			  << temperature << ", dt=" << dt << ", n_iters=" << n_iters
			  << ", n_iters_th=" << n_iters_th << ", mass=" << mass
			  << ", bias=" << bias << ", skip=" << skip << "\n";
	std::cout << std::endl;
}

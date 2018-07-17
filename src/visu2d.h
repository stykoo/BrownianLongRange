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
 * \file visu.h
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \version 1.0
 * \date 2017-11-27
 * \brief Visualization of the system
 *
 * Header file for visu.cpp.
*/

#ifndef BROWNIANLONGRANGE_VISU2D_H_
#define BROWNIANLONGRANGE_VISU2D_H_

#include <SFML/Graphics.hpp>
#include "state.h"

class Visu {
	public:
		Visu(const State *state, const long n_parts,
		     const long n_parts_1);
		void run();

	private:
		// Constant variables for visualization
		// const int windowSize = 600; //!< Size of the window
		const int FPS = 24; //!< Number of frames per second

		const State *state; //!< Pointer to the state of the system
		const long n_parts; //!< Number of particles
		const long n_parts_1; //!< Number of particles of species 1
};

#endif // BROWNIANLONGRANGE2D_VISU_H_

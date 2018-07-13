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
 * \file simul.cpp
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Visualization of the system
 *
 * The system is visualized using the SFML library.
*/

#include <iostream>
#include <cmath>
#include "visu2d.h"

/*!
 * \brief Constructor for visualization
 *
 * \param state Pointer to the state of the system
 * \param len Length of the box
 * \param n_parts Number of particles
 * \param n_parts Number of particles of species 1
 */
Visu::Visu(const State *state, const double len,
           const long n_parts, const long n_parts_1) :
	state(state), len(len), n_parts(n_parts), n_parts_1(n_parts_1) {
}

/*!
 * \brief Thread for visualization.
 *
 * Open a window, draw the particles and update their
 * positions at a certain number of FPS while the simulation is runing.
 */
void Visu::run() {
	sf::VideoMode mode = sf::VideoMode::getDesktopMode();
	const float windowSize = std::min(mode.width, mode.height) * 9 / 10;

    sf::RenderWindow window;
    window.create(sf::VideoMode(windowSize, windowSize),
	              "Brownian Particles with Coulomb interactions");

	float scale = windowSize / len;
	// We assume that the particles have diameter 1
    sf::CircleShape circle1(scale / 2.0);
    sf::CircleShape circle2(scale / 2.0);
    circle1.setFillColor(sf::Color::Transparent);
    circle2.setFillColor(sf::Color::Transparent);
    circle1.setOutlineThickness(-3);
    circle2.setOutlineThickness(-3);
    circle1.setOutlineColor(sf::Color::Red);
    circle2.setOutlineColor(sf::Color::Blue);

    window.setFramerateLimit(FPS);

	const std::vector<double> & pos = state->getPos();

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(sf::Color::White);

		for (long i = 0 ; i < n_parts ; ++i) {
			float x = pos[i] * scale; 
			pbc(x, (float) windowSize);
			float y = pos[n_parts + i] * scale; 
			pbc(y, (float) windowSize);

			// Draw multiple times if on the boundary
			int per_x = (x > windowSize - scale);
			int per_y = (y > windowSize - scale);

			for (int px = 0 ; px <= per_x ; ++px) {
				for (int py = 0 ; py <= per_y ; ++py) {
					if (i < n_parts_1) {
						circle1.setPosition(x - px * windowSize,
										    y - py * windowSize);
						window.draw(circle1);
					} else {
						circle2.setPosition(x - px * windowSize,
							   			    y - py * windowSize);
						window.draw(circle2);
					}
				}
			}
		}
        window.display();
    }
}

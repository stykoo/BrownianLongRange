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
 * \file visu3d.cpp
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \brief Visualization of the system in 3d.
 *
 * The system is visualized using the VTK library.
*/

#include <iostream>
#include <vector>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkCubeSource.h>
#include <vtkAxesActor.h>
#include <vtkTransform.h>
#include <vtkSmartPointer.h>
#include <vtkCommand.h>
#include <vtkProperty.h>

#include "visu3d.h"

template<typename T>
using vSP = vtkSmartPointer<T>; // Shortcut

struct Visu3dTimer : public vtkCommand {
	public:
		vtkTypeMacro(Visu3dTimer, vtkCommand);
		static Visu3dTimer *New() {
    		return new Visu3dTimer();
  		}

		void Execute(vtkObject *caller, unsigned long vtkNotUsed(eventId),
					 void *vtkNotUsed(callData)) {
			const std::vector<double> & pos = state->getPos();

			vtkRenderWindowInteractor *iren =
				static_cast<vtkRenderWindowInteractor*>(caller);

			// Update the positions and the colors
			for (size_t i = 0 ; i < sphereActors->size() ; ++i) {
				(*sphereActors)[i]->SetPosition(pos[i],
                                                pos[n_parts+i],
											    pos[2*n_parts+i]);
			}

			iren->Render();
		}

	// Public attributes (this is dirty...)
	long n_parts;
	std::vector< vSP<vtkActor> > *sphereActors;
	const State *state;
};

/*!
 * \brief Constructor for visualization
 *
 * \param state Pointer to the state of the system
 * \param n_parts Number of particles
 * \param n_parts_1 Number of particles of species 1
 */
Visu3d::Visu3d(const State *state, const long n_parts, const long n_parts_1) :
	state(state), n_parts(n_parts), n_parts_1(n_parts_1) {
}

/*!
 * \brief Thread for visualization.
 *
 * Open a window, draw the particles and update their
 * positions while the simulation is runing.
 */
void Visu3d::run() {
	const std::vector<double> & pos = state->getPos();

	std::vector< vSP<vtkSphereSource> > sphereSources;
	std::vector< vSP<vtkPolyDataMapper> > sphereMappers;
	std::vector< vSP<vtkActor> > sphereActors;

	for (long i = 0 ; i < n_parts ; ++i) {
		sphereSources.push_back(vSP<vtkSphereSource>::New());
		sphereSources[i]->SetRadius(0.5);
		sphereSources[i]->SetThetaResolution(sphere_res);
		sphereSources[i]->SetPhiResolution(sphere_res);

		//create a mapper
		sphereMappers.push_back(vSP<vtkPolyDataMapper>::New());
		sphereMappers[i]->SetInputConnection(sphereSources[i]->GetOutputPort());

		// create an actor
		sphereActors.push_back(vSP<vtkActor>::New());
		sphereActors[i]->SetMapper(sphereMappers[i]);
		sphereActors[i]->SetPosition(pos[i], pos[n_parts+i], pos[2*n_parts+i]);
		sphereActors[i]->GetProperty()->SetColor(1.0 * (i < n_parts_1),
				                                 0.0,
												 1.0 * (i >= n_parts_1));
		sphereActors[i]->GetProperty()->SetOpacity(sphere_opa);

	}

	// Box
	vSP<vtkCubeSource> cubeSource = vSP<vtkCubeSource>::New();
	cubeSource->SetXLength(1.0);
	cubeSource->SetYLength(1.0);
	cubeSource->SetZLength(1.0);
	cubeSource->SetCenter(0.5, 0.5, 0.5);
	cubeSource->Update();
	vSP<vtkPolyDataMapper> cubeMapper = vSP<vtkPolyDataMapper>::New();
	cubeMapper->SetInputConnection(cubeSource->GetOutputPort());
	vSP<vtkActor> cubeActor = vSP<vtkActor>::New();
	cubeActor->GetProperty()->SetRepresentationToWireframe();
	cubeActor->SetMapper(cubeMapper);

	// Axes
	vSP<vtkAxesActor> axes = vSP<vtkAxesActor>::New();
	axes->SetXAxisLabelText("");
	axes->SetYAxisLabelText("");
	axes->SetZAxisLabelText("");

	// Renderer and render window
	vSP<vtkRenderer> renderer = vSP<vtkRenderer>::New();
	vSP<vtkRenderWindow> renderWindow = vSP<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);

	// Interactor
	vSP<vtkRenderWindowInteractor> renderWindowInteractor =
		vSP<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Add the actors to the scene
	for (long i = 0 ; i < n_parts ; ++i) {
		renderer->AddActor(sphereActors[i]);
	}
	renderer->AddActor(cubeActor);
	renderer->AddActor(axes);
	renderer->SetBackground(.1,.2,.3); // Background dark blue

	// Take care of animation
	renderWindowInteractor->Initialize();
	renderWindowInteractor->CreateRepeatingTimer(delay);
	vtkSmartPointer<Visu3dTimer> timerCallback =
		vtkSmartPointer<Visu3dTimer>::New();
	timerCallback->n_parts = n_parts;
	timerCallback->sphereActors = &sphereActors;
	timerCallback->state = state;
	renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, timerCallback);

	// Render and interact
	renderWindow->Render();
	renderWindowInteractor->Start();
}

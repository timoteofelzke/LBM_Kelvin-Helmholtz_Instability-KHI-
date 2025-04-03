/**
 *  @file   initial.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Header of the functions to read and apply initial conditions
 */

#ifndef __INITIAL_H_INCLUDED__
#define __INITIAL_H_INCLUDED__

#include <string>

struct stGeometry; // Foward Declare of stGeometry (see geometry.h for more details)

/**
 * @brief Allocates memory and set the initial (homogenous) condition for the distribution function
 * @param geo reference to the stGeometry struct containing the information about the medium/domain.
 * @param rho value of the density to be set in all fluid nodes.
 * @return the position in memory of the allocated distriubution function
 * 
 * This function allocates memory for all (fluid) sites and it sets and null velocity with homosgenous density in the fluid.
 **/

double * setInitial(stGeometry& geo, double rho);

/**
 * @brief Allocates memory and set the initial condition (density and velocity) from a vtk file.
 *
 * @param geo reference to the stGeometry struct containing the information about the medium/domain.
 * @param rhoFilename String with the name of vtk file with density field to be set as initial condition
 * @param velFilename String with the name of vtk file with velocity field to be set as initial condition
 * @return the position in memory of the allocated distriubution function.
 *
 * Allocates memory for the distributions functions and set the initial condition  (density and velocity) from a vtk file.
 **/

double * setInitialFromVtk(stGeometry& geo,std::string rhoFilename, std::string velFilename);

#endif

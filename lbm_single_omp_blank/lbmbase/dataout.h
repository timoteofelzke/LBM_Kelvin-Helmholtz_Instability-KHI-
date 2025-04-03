/**
 *  @file   dataout.h
 *  @author Diogo Nardelli Siebert
 *  @brief Header of functions to write output files (VTK) of  the flow fields (density and velocity)
 */

#ifndef __DATAOUT_H_INCLUDED__   // if datain.h hasn't been included yet...
#define __DATAOUT_H_INCLUDED__

#include <string>
#include "vtkio/vtkio.h"

struct stGeometry;
struct stForce; // Foward Declared

/**
 *  Write the density field to a vtk file.
 *
 *  @param      geo                 Pointer to the array containg the information about the medium geometry
 *  @param      iniN                Pointer to the array containg the distribution function (population) of all fluid sites.
 *  @param      step                The current step (time iteration) of the simulation.
 *  @param      fileFormat          The format of vtk to used ASCII or BINARY.
 *  @param      prefix              The prefix using in the filename          
 **/

void vtkDensity( stGeometry& geo, double *ini_N, unsigned int step, int fileFormat, std::string prefix = "rho");

/**
 *  Write the velocity field to a vtk file.
 *
 *  @param      geo                 Pointer to the array containg the information about the medium geometry
 *  @param      iniN                Pointer to the array containg the distribution function (population) of all fluid sites.
 *  @param      step                The current step (time iteration) of the simulation.
 *  @param      fileFormat          The format of vtk to used ASCII (vtkASCII = 0) or BINARY (vtkBinary = 1)
 **/

void vtkVelocity(stGeometry& geo, double *ini_N, stForce& force , unsigned int step, int fileFormat = vtkASCII);

/**
 *  Convert a integer to string (with zeros to the left)
 *
 *  @param      number      Integer number to be converted.
 *  @return     String containing the integer.
 */

std::string intToString(int number);

/**
 *  Backup the distributions functions to a file.
 *
 *  @param      filename            Name of the file where the backup will be saved.
 *  @param      iniN                 Pointer to the array containg the distribution function (population) of all fluid sites.
 *  @param      numberOfPoints      The number of fluid points in the simulation.
 *  @param      step                The current step (time iteration) of the simulation.
 **/

void saveRecovery(const char* filename, double *iniN, int numberOfPoints, unsigned int step);

/**
 *  Load the distributions functions from a files.
 *
 *  @param      filename            Name of the file where the backup is saved.
 *  @param      iniN                Pointer to the array containg the distribution function (population) of all fluid sites.
 *  @param      numberOfPoints      The number of fluid points in the simulation.
 *  @param      step                The current step (time iteration) of the simulation.
 **/

void loadRecovery(const char* filename, double *pDistribution, int numberOfPoints, unsigned int& step);

void vtkDensity( stGeometry& geo, double *ini_R, double* ini_B, unsigned int step, int fileFormat);

/**
 *  Write the density field to a vtk file.
 *
 *  @param      geo                 Pointer to the array containg the information about the medium geometry
 *  @param      iniN                Pointer to the array containg the distribution function (population) of all fluid sites.
 *  @param      step                The current step (time iteration) of the simulation.
 *  @param      fileFormat          The format of vtk to used ASCII or BINARY.
 **/

void vtkConcentration( stGeometry& geo, double *ini_R,double* ini_B, unsigned int step, int fileFormat);

/**
 *  Convert a integer to string (with zeros to the left)
 *
 *  @param      number      Integer number to be converted.
 *  @return     String containing the integer.
 */


void vtkVelocity( stGeometry& geo, double *ini_R, double *ini_B, unsigned int step, int fileFormat);

#endif

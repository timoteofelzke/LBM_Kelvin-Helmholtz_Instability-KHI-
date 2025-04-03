/**
 *  @file   parameters.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Header of functions needed to read INI files with simulation parameters
 *
 *  Header of functions needed to read INI files with simulation parameters
 *  It is important to nota that command line options are parsed/handled using options.h functions
 */

#ifndef __DATAIN_H_INCLUDED__   // if datain.h hasn't been included yet...
#define __DATAIN_H_INCLUDED__

#include <string>

struct stParameters
{
    std::string geoFilename;      /*!< The Name of file with the geometry */
    std::string geoFiletype;      /*!< The type of geometry file (ex: vtk, raw) */
    std::string iniRhoFilename;   /*!< The name of the file with density field for the initial condition  */
    std::string iniVelFilename;   /*!< The name of the file with velocity field for the initial condition    */
    int geoSizeX;                 /*!< Length (in number of sites) of the geometry in the X axis (used only for reading raw files) */
    int geoSizeY;                 /*!< Length (in number of sites) of the geometry in the Y axis (used only for reading raw files) */
    int geoSizeZ;                 /*!< Length (in number of sites) of the geometry in the Z axis (used only for reading raw files) */
    double tau;                   /*!< Collision relaxation Time for BGK operator */
    double viscosity;             /*!< Value of the viscosity for the BGK operator */
    double initialRho;            /*!< Initial valuer for the density (if it is uniform r */
    int numberOfSteps;            /*!< The number of step in the simulation */
    int numberOfFiles;            /*!< The number of output (vtk) file recorder in the simulation */
    int recordInterval;           /*!< The interval of number of steps to record the output (vtk) files  */
    int fileFormat;               /*!< The format of the geometry file */
    bool fileDensity;             /*!< A boolean to indicate if a vtk output for the density must be generated */
    bool fileVelocity;            /*!< A boolean to indicate if a vtk output for the velocity must be generated */
    bool setInitial;              /*!< A boolean to indicate if initial conditions from a vtk file will be used */
    bool writeZero;               /*!< A boolean to indicate if the output for the initial set must be generated */
    bool writeFinal;              /*!< A boolean to indicate if the output for the final step set must be generated */
    double accX , accY, accZ;     /*!< The external acceleration field */
};

/**
 *  @bfief Read the parameters store in a INI file and store the information in a stParameters struct.
 *
 *  @param  filename        String with the name of the INI file.
 *  @return A stParameters struct with information read from the INI file.
 *
 *  Read the parameters store in a INI file and store the information in a stParameters struct.
 **/

stParameters readDataINI(std::string filename);

/**
 *  @bfief Print to the standard output (screen) the value of the parameters recorded in a stParameters struct
 *
 *  @param prm The stParameters struct with the parameters information (read from the INI file)
 **/

void reportParameters(stParameters prm) ;

#endif

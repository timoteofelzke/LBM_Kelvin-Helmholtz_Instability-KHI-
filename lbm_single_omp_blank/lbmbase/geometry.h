/**
 *  @file    geometry.h
 *  @author  Diogo Nardelli Siebert
 *  @brief   Header of the functions to read and manage geometry information.
 */

#ifndef __GEOMETRY_H_INCLUDED__
#define __GEOMETRY_H_INCLUDED__

#define SOLID 0
#define FLUID 1

#include <string>

struct stGeometry
{
    unsigned int *position;    /*!< Array with the position of each indexed (fluid) site */
    unsigned int *index;       /*!< Index of all sites (solid sites have index 0 and fluid sites are numbered starting from 1) */
    int numberOfPoints;        /*!< Total number of fluid sites in the geometry */
    int nx;                    /*!< Length (in number of sites) of the geometry in the X axis */
    int ny;                    /*!< Length (in number of sites) of the geometry in the Y axis  */
    int nz;                    /*!< Length (in number of sites) of the geometry in the Z axis  */
};

/**
 *  Computes the coordinate (x,y,z) for a given site index.  It is important to note that only fluid nodes
 *  are indexed and the enumeration starts with one. All solid nodes recieve the index zero.
 *
 *  @param  geo                 Reference to the stGeometry containing all the geometry information
 *  @param  index               Index of the site whose position must be determined.
 *  @param  x                   Reference to the variable which will store the x coordinate ot the position
 *  @param  y                   Reference to the variable which will store the y coordinate ot the position
 *  @param  z                   Reference to the variable which will store the z coordinate ot the position
 */

void getPosition(stGeometry &geo,int id, int& x,int &y,int &z);

unsigned int getPosition(stGeometry &geo,int id);

/**
 *  Computes the index of the site for a given position (x,y,z). It is important to note that only fluid nodes
 *  are indexed and the enumeration starts with one. All solid nodes recieve the index zero.
 *
 *  @param  geo                 Reference to the stGeometry containing all the geometry information
 *  @param  x                   Reference to x coordinate ot the site.
 *  @param  y                   Reference to y coordinate ot the site.
 *  @param  z                   Reference to z coordinate ot the site.
 *  @return                     Index of the position if it is a fluid node and zero if it is a solid node.
 */

unsigned int getIndex(stGeometry &geo,int x,int y, int z);

stGeometry createFullGeometry(int nx, int ny, int nz);

stGeometry readGeometry(std::string filename ,std::string filetype, int sizeX, int sizeY, int sizeZ); 

/**
 *  Read the geometry, i.e., the 3D matrix containing information whether a pixel is a solid (0) or a fluid (1) . This
 *  function can be used with both Binary or ASCII VTK. The data type must be either INT (integer of 32 bits) or CHAR
 *  (unsigned integer of 8 bits)
 *
 *  @param  filename        A string with name of the VTK file.
 *  @return                 An stGeometry struct with all the (organized) information about the geometry.
 **/

stGeometry readGeometryVTK(std::string filename);

/**
 *  Read the geometry, i.e., the 3D matrix containing information whether a pixel is a solid (0) or a fluid (1). The data
 *  must be stored as unsigned int of 8 bits. No offset in the begining of the file is allowed, i. e., the size of the file
 *  must be exactly Nx*Ny*Nz*8 bits.
 *
 *  @param  filename        A string with name of the RAW file.
 *  @param  nx              Length (in number of pixels) of the image in the X axis
 *  @param  ny              Length (in number of pixels) of the image in the X axis
 *  @param  nz              Length (in number of pixels) of the image in the X axi
 *  @return                 An stGeometry struct with all the (organized) information about the geometry.s
 **/

stGeometry readGeometryRAW(std::string geoFilename,int nx,int ny,int nz);

/**
 *  Print information about the geometyr in the standart output (screen).
 *
 *  @param  geo                 Reference to the stGeometry containing all the geometry information
 **/

void reportGeometry(stGeometry& geo);

/**
 *  Update the position of the sites in a stGeometry stuct  if the geometry has been reindexed.
 *
 *  @param  geo                 Reference to the stGeometry containing all the geometry information
 **/

void updateGeometry(stGeometry& geo);

#endif

/**
 * @author Diogo Nardellli Siebert
 * @date 27/07/15
 * @file vtkio.h
 * @brief Header file for functions and data structure used for reading and writing to VTK files.
 */

#ifndef __VTKIO_H_INCLUDED__
#define __VTKIO_H_INCLUDED__

#include <string>
#include <fstream>
#include <stdint.h>

enum { vtkASCII = 0 , vtkBinary = 1, xmlASCII = 2, xmlBinary = 3, xmlCompressed = 4};
enum vtkDataFormat { vtkScalar , vtkVector };

/**
 * @struct stVtkHeader
 * @brief Stores the information contained the header of a VTK legacy file
 * 
 * This struct contain variables for all the field present in the header of 
 * VTK Legacy file. It is used to pass information between different functions
 * which read data from a VTK.  
 *
 **/

bool isLittleEndian();

struct stVtkHeader
{
    std::string vtkVersion;               /*!< The Vtk File Format Version of the file */
    std::string fileTitle;                /*!< The title of the file (do not confuse with the name of the file) */
    std::string fileFormat;               /*!< If the format that data is stored (BINARY OR ASCII) */
    std::string dataSetType;              /*!< The type of geometry (grid) that data is associeted to (STRUCTURED GRID for LBM applications) */
    int   sizeX,  sizeY,  sizeZ;          /*!< Lenght in pixels of the image in each axis */
    int  ratioX, ratioY, ratioZ;          /*!< Ratio of the different axis */
    int originX,originY,originZ;          /*!< Position of the origin of the image */
    int64_t pointData;                        /*!< Number of points in the geometry  */
    int64_t dataSize;                         /*!< Number of data store in the file, i. e., number of points multiplied by the number of components, */
    std::string dataType;                 /*!< The type of the field store in the file (SCALARS or VECTORS) */
    std::string dataName;                 /*!< The name of the data field stored in the file */
    std::string variableType;             /*!< The type of the data (ex: float, int, char) */
    std::string lookupTable;              /*!< ? */
    std::string byteOrder; 
    bool xml;
    bool compress;
    int vtkType;
    unsigned char * dataBuffer;
    unsigned int dataOffset;
    unsigned int dataCount;
    unsigned int dataUnit;
    std::streampos headerSize;

};

/**
 *  @brief Exchange a float from big endian to little endian and vice-versa.
 *  @param  data        Pointer to the float to be inverted.
 *
 *  Invert a float number store in memory byte by byte (big endian <-> little endian) to be compatible with the VTK format.
 **/

void vtkFloatSwap( float* pData);

/**
 *  @brief  Exchange a short int from big endian to little endian and vice-versa.
 *  @param  data        Pointer to the float to be inverted.
 *
 *  Invert a short number store in memory byte by byte (big endian <-> little endian) to be compatible with the VTK format.
 **/

void vtkShortSwap( short* pData);

/**
 *  @brief Creates a stVtkHeader struct with the  given information.
 *
 *  @param  fileTitle       The title of the file (do not confuse with the name of the file)
 *  @param  fileFormat      If the format that data is stored (BINARY OR ASCII)
 *  @param  dataType        The type of the field store in the file (SCALARS or VECTORS)
 *  @param  variableType    The type of the data (ex: float, int, char)
 *  @param  sizeX           Lenght in pixels of the image in the X axis
 *  @param  sizeY           Lenght in pixels of the image in the Y axis
 *  @param  sizeZ           Lenght in pixels of the image in the Z axis
 *
 *  @return  A stVtkHeader struct with the  given information.
 **/

stVtkHeader setVtkHeader(std::string fileTitle,int fileFormat,std::string dataType
			,std::string dataName,std::string variableType,int sizeX,int sizeY
            ,int sizeZ);

/**
 *  @brief Generates the header in the ASCII format and writes the generated string to a streaming output.
 *
 *  @param  file        The output streaming
 *  @param  header      stVtkHeader with the header information.
 **/

void writeVtkHeader(std::ostream &file, stVtkHeader& header);
void writeVtkLegacyHeader(std::ostream &file, stVtkHeader& header);
void writeVtkXMLHeader(std::ostream &file, stVtkHeader& header);

/**
 * @brief Read the header of a VTK from a file 
 * @param vtkFile   A ifstream object created from a VTK legacy file. To avoid problems, always open the file in binary mode (std::ifstream::binary).
 * 
 * This function read all information in a header of a VTK legacy file and stores the data
 * in a stVtkHeader struct.  
 */

stVtkHeader readVtkHeader(std::ifstream& vtkFile);

/**
 * @brief Read a float from  a vtk file.
 *
 * @param vtkFile   A ifstream object created from a VTK legacy file. To avoid problems, always open the file in binary mode (std::ifstream::binary).
 * @param header    The header containing all the information about the data stored in the file. (See documention of the function getVtkHeader)

 * @return Return the next the double in the current vtkFile ifstream position.
 * 
 * This function reads a single float contained in a VTK file (the float in the current position of the vtkFile ifstream).
 * If the file format is ASCII it will work even for double data. If the file format is binary than will only work if data in the file 
 * is really saved in float.
 *
 * @return return the data read into as a double.
 */

double readVtkFloat(std::ifstream &vtkFile,stVtkHeader& header);       

/**
 * @brief Write a float to a vtk file.
 *
 * @param vtkFile   A ifstream object created from a VTK legacy file. To avoid problems, always open the file in binary mode (std::ifstream::binary).
 * @param header    The header containing all the information about the data that will be stored in the file. (See documention of the function getVtkHeader)
 * @param data      The float data that is going to be writting in the file.
 *
 * This function writes a single float to a VTK file.
 */

void writeVtkFloat(std::ostream &vtkFile,stVtkHeader& header, float data);

/**
 * @brief Read a char from a vtk file.
 *
 * @param vtkFile   A ifstream object created from a VTK legacy file. To avoid problems, always open the file in binary mode (std::ifstream::binary).
 * @param header    The header containing all the information about the data stored in the file. (See documention of the function getVtkHeader)

 * @return Return the next char in the current vtkFile ifstream position.
 *
 * This function reads a single char contained in a VTK file (the char in the current position of the vtkFile ifstream).
 *
 * @return return the read char.
 */

char readVtkChar(std::ifstream &vtkFile,stVtkHeader& header);

/**
 * @brief Read a int from a vtk file.
 *
 * @param vtkFile   A ifstream object created from a VTK legacy file. To avoid problems, always open the file in binary mode (std::ifstream::binary).
 * @param header    The header containing all the information about the data stored in the file. (See documention of the function getVtkHeader)

 * @return Return the next int in the current vtkFile ifstream position.
 *
 * This function reads a single int contained in a VTK file (the int in the current position of the vtkFile ifstream).
 *
 * @return return the variable read.
 */

int readVtkInt(std::ifstream &vtkFile,stVtkHeader& header);

void writeVtkShort(std::ostream &vtkFile,stVtkHeader& header, short data);
void writeVtkFooter(std::ostream &file, stVtkHeader& header);
void writeVtkXMLFooter(std::ostream &file, stVtkHeader& header);

unsigned char readVtkUnsignedChar(std::ifstream &vtkFile,stVtkHeader& header);

unsigned char* getVtkUnsignedChar(std::ifstream &vtkFile,stVtkHeader& header,unsigned long int size);

void vtkIntSwap( int* pData);


#endif

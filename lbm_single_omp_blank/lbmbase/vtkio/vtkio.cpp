#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <zlib.h>   
#include "vtkio.h"

/**
 * @author Diogo Nardellli Siebert
 * @date 27/07/15
 * @file vtkio.cpp
 * @brief Implementation of functions for reading and writing to VTK files.
 */

using namespace std;

static char b64table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

unsigned char *spc_base64_encode(unsigned char *input, size_t len, int wrap) 
{
	unsigned char *output, *p;
	size_t        i = 0, mod = len % 3, toalloc;
	
	toalloc = (len / 3) * 4 + (3 - mod) % 3 + 1;
	if (wrap) 
	{
		toalloc += len / 57;
		if (len % 57) toalloc++;    
	}
	
	p = output = (unsigned char *) malloc(((len / 3) + (mod ? 1 : 0)) * 4 + 1);
	if (!p) return 0;
	while (i < len - mod) 
	{
		*p++ = b64table[input[i++] >> 2];
		*p++ = b64table[((input[i - 1] << 4) | (input[i] >> 4)) & 0x3f];
		*p++ = b64table[((input[i] << 2) | (input[i + 1] >> 6)) & 0x3f];
		*p++ = b64table[input[i + 1] & 0x3f];
		i += 2;
		if (wrap && !(i % 57)) *p++ = '\n';
	}
	if (!mod) 
	{
		if (wrap && i % 57) *p++ = '\n';
		*p = 0;
		return output;
	} 
	else 
	{
		*p++ = b64table[input[i++] >> 2];
		*p++ = b64table[((input[i - 1] << 4) | (input[i] >> 4)) & 0x3f];
		if (mod == 1) 
		{
		*p++ = '=';
		*p++ = '=';
		if (wrap) *p++ = '\n';
		*p = 0;
		return output;
		} 
		else 
		{
			*p++ = b64table[(input[i] << 2) & 0x3f];
			*p++ = '=';
			if (wrap) *p++ = '\n';
			*p = 0;
			return output;
		}
	}
}

stVtkHeader readVtkHeader( ifstream& vtkFile)
{
	stVtkHeader header;

	string keyword;

	getline(vtkFile, header.vtkVersion);    // Read First Line (Comment about vtk file version)
	getline(vtkFile, header.fileTitle);

	vtkFile >> header.fileFormat;
	vtkFile >> keyword >> header.dataSetType;
	vtkFile >> keyword >> header.sizeX >> header.sizeY >> header.sizeZ;
	vtkFile >> keyword >> header.ratioX >> header.ratioY >> header.ratioZ;
	vtkFile >> keyword >> header.originX >> header.originY >> header.originZ;
	vtkFile >> keyword >> header.pointData;
	vtkFile >> header.dataType;
	vtkFile >> header.dataName >> header.variableType;

	//vtkFile >> keyword >> header.lookupTable;

	// Computes the amount of data stored in the file (multiple the pointdata by the number of components)
	if (header.dataType == "SCALARS") 
	{
			header.dataSize = 1 * header.pointData;
			vtkFile >> keyword >> header.lookupTable;
	}
	else if (header.dataType == "VECTORS") header.dataSize = 3 * header.pointData;

	// Read the new line char in the end of the line to avoid problem reading this instead of the float number
	if (header.fileFormat == "BINARY")
	{
		char c;
		vtkFile.read(&c, sizeof(char) );	
	}
	
	header.headerSize = vtkFile.tellg();
	
	return header;
}

char readVtkChar(ifstream &vtkFile,stVtkHeader& header)
{
	char data;
	if (header.fileFormat == "ASCII")
	{
		vtkFile >> data;
		data = static_cast<char>( atoi( &data)  );
	}
	else if(header.fileFormat == "BINARY")
		vtkFile.read( (char*)  &data, sizeof(char) );
	return data;
}

unsigned char readVtkUnsignedChar(ifstream &vtkFile,stVtkHeader& header)
{
        unsigned char data;
    if (header.fileFormat == "ASCII")
            vtkFile >> data;
    else if(header.fileFormat == "BINARY")
            vtkFile.read( (char*)  &data, sizeof(unsigned char) );
    return data;
}

unsigned char* getVtkUnsignedChar(ifstream &vtkFile,stVtkHeader& header,unsigned long int size)
{
    unsigned char* data = new unsigned char[size];
    if (header.fileFormat == "BINARY")
            vtkFile.read( (char*) data, size * sizeof(unsigned char) );
    return data;
}

int readVtkInt(ifstream &vtkFile,stVtkHeader& header)
{
	int data;
	if (header.fileFormat == "ASCII")
		vtkFile >> data;
	else if(header.fileFormat == "BINARY")
    {
		vtkFile.read( (char*) &data, sizeof(int) );
        vtkIntSwap( &data );
    }
	return data;
}

double readVtkFloat(ifstream &vtkFile,stVtkHeader& header)
{
	float data;
	
	// Verify if the format of the file is ASCII or BINARY and read the data accordingly
	if (header.fileFormat == "ASCII")
		vtkFile >> data;
	else if(header.fileFormat == "BINARY")
	{
		vtkFile.read( (char*) (&data), sizeof(float) );
		vtkFloatSwap( &data );
	}
	
	return data;
}

void vtkShortSwap( short* pData)
{
	char oneByte; //Temporary variable to store a char in the swap.
	char* swap = (char*) pData; // Convert the pointer to a char* to be possible to point to individual bytes in the float.

	// Invert the 4 bytes in the float number  0123 -> 3210
	oneByte = swap[0];
	swap[0] = swap[1];
	swap[1] = oneByte;
}

void vtkFloatSwap( float* pData)
{
	char oneByte; //Temporary variable to store a char in the swap.
	char* swap = (char*) pData; // Convert the pointer to a char* to be possible to point to individual bytes in the float.

	// Invert the 4 bytes in the float number  0123 -> 3210
	oneByte = swap[0];
	swap[0] = swap[3];
	swap[3] = oneByte;
	oneByte = swap[1];
	swap[1] = swap[2];
	swap[2] = oneByte;
}

void vtkIntSwap( int* pData)
{
    char oneByte; //Temporary variable to store a char in the swap.
    char* swap = (char*) pData; // Convert the pointer to a char* to be possible to point to individual bytes in the float.

    // Invert the 4 bytes in the float number  0123 -> 3210
    oneByte = swap[0];
    swap[0] = swap[3];
    swap[3] = oneByte;
    oneByte = swap[1];
    swap[1] = swap[2];
    swap[2] = oneByte;
}

void writeVtkHeader(ostream &file, stVtkHeader& header)
{
    if (header.xml) writeVtkXMLHeader(file,header);
    else writeVtkLegacyHeader(file,header);
}

void writeVtkLegacyHeader(ostream &file, stVtkHeader& header)
{
	file << header.vtkVersion << endl;
	file << header.fileTitle << endl;
	file << header.fileFormat << endl;
	file << "DATASET " << header.dataSetType << endl;
	file << "DIMENSIONS " << header.sizeX << " " << header.sizeY << " " << header.sizeZ << endl;
	file << "SPACING " << header.ratioX << " " << header.ratioY << " " << header.ratioZ << endl;
	file << "ORIGIN " << header.originX << " " << header.originY << " " << header.originZ << endl;
	file << "POINT_DATA " << header.pointData << endl;
	file << header.dataType << " " << header.dataName << " " << header.variableType << endl;
	if ( header.dataType == "SCALARS" ) file << "LOOKUP_TABLE " << header.lookupTable << endl;
}

void writeVtkXMLHeader(ostream &file, stVtkHeader& header)
{
 	file << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl;
	file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"" << header.byteOrder << "\"";
	if ( header.compress ) file << " compressor=\"vtkZLibDataCompressor\"";
	file << ">" << endl;

	file << "\t<ImageData WholeExtent=\" 0 " << header.sizeX - 1 << " 0 " << header.sizeY - 1 << " 0 " << header.sizeZ - 1 << "\"";
	file << " Origin=\""  << header.originX << " " << header.originY << " " << header.originZ << "\"";
	file << " Spacing=\"" << header.ratioX << " " << header.ratioY << " " << header.ratioZ << "\">" << endl;
        
	file << "\t\t<Piece Extent= \"0 " << header.sizeX - 1 << " 0 " << header.sizeY - 1 << " 0 " << header.sizeZ - 1 << "\">" << endl;
	
	file << "\t\t\t<PointData";

	if (header.dataType == "SCALARS") 
	{
		file << " Scalars=";		
	}
	else if (header.dataType == "VECTORS") file << " Vectors=";
	
	file << "\"" << header.dataName <<"\">" << endl;
	
	file << "\t\t\t\t<DataArray Name=\"" << header.dataName << "\"" << " type=";
	
    if (header.variableType == "float" ) { file << "\"Float32\""; header.dataUnit = sizeof(float); }
    if (header.variableType == "int"   ) { file << "\"Int32\""; header.dataUnit = sizeof(int); }
    if (header.variableType == "char"  ) { file << "\"UInt8\"";  header.dataUnit = sizeof(char); }

    if (header.dataType == "VECTORS") file << " NumberOfComponents=\"3\"";

    file << " format=";
	if (header.fileFormat == "ASCII" ) file << "\"ascii\">";
	if (header.fileFormat == "BINARY") 
	{
		header.dataOffset = (header.compress ? 0 : 1) * sizeof(unsigned int);
		header.dataBuffer = new unsigned char[header.dataOffset  + header.dataSize * header.dataUnit]; 
        file << "\"binary\">";
	}
}

void writeVtkFooter(ostream &file, stVtkHeader& header)
{
    if (header.xml) writeVtkXMLFooter(file, header);
}

void writeVtkXMLFooter(ostream &file, stVtkHeader& header)
{
	if (header.fileFormat == "BINARY") 
	{
		if (header.compress)
		{
			char* compressedData = new char[ header.dataCount * header.dataUnit] ;
 
			z_stream defstream;
			defstream.zalloc = Z_NULL; defstream.zfree = Z_NULL; defstream.opaque = Z_NULL;
			defstream.avail_in = (uInt)  header.dataSize * header.dataUnit; 
			defstream.next_in = (Bytef *) header.dataBuffer; 
			defstream.avail_out = (uInt) header.dataSize * header.dataUnit; // size of output
			defstream.next_out = (Bytef *) compressedData ;     
			// the actual compression work.

			deflateInit(&defstream, Z_BEST_COMPRESSION);
			deflate(&defstream, Z_FINISH);
			deflateEnd(&defstream);

			unsigned int compressedSize = (char*) defstream.next_out - compressedData;
            unsigned int compressedInfo[] = {1, header.dataCount * header.dataUnit, header.dataCount* header.dataUnit, compressedSize };
	
            unsigned char* base64Data =  spc_base64_encode( (unsigned char*) compressedInfo , 4*sizeof(unsigned int) ,  0);
            delete header.dataBuffer;
            file << base64Data;
            free(base64Data);
            base64Data =   spc_base64_encode( (unsigned char*) compressedData ,compressedSize , 0 );
            file << base64Data;
            free(base64Data);
			delete compressedData;
		}
		else
		{
			int* dataSize = (int*) header.dataBuffer;
			*dataSize = header.dataSize * header.dataUnit;
            unsigned char* base64Data = spc_base64_encode( header.dataBuffer , header.dataOffset + header.dataSize*header.dataUnit , 0 );
            delete header.dataBuffer;
            file << base64Data;
            free(base64Data);
		}	

	}
    file << "</DataArray>" << endl;
    file << "\t\t\t</PointData>" << endl;
	file << "\t\t</Piece>" << endl;
	file << "\t</ImageData>" << endl;
	file << "</VTKFile>" << endl;
}
	
stVtkHeader setVtkHeader(string fileTitle,int fileFormat,string dataType,string dataName,string variableType,int sizeX,int sizeY,int sizeZ)
{
	stVtkHeader header;
	
	header.fileTitle = fileTitle;   
    header.dataType = dataType;
	header.dataName = dataName;
	header.variableType = variableType;
	header.sizeX = sizeX;
	header.sizeY = sizeY;
	header.sizeZ = sizeZ;
	header.vtkVersion = "# vtk DataFile Version 2.0";
	header.ratioX = 1;
	header.ratioY = 1;
	header.ratioZ = 1;
	header.originX = 0;
	header.originY = 0;
	header.originZ = 0;
	header.pointData = static_cast<int64_t>(header.sizeX) *  static_cast<int64_t>(header.sizeY) *  static_cast<int64_t>(header.sizeZ);
	header.dataSetType = "STRUCTURED_POINTS";
	header.lookupTable = "default";    

    switch ( fileFormat )
    {
        case vtkASCII:
            header.fileFormat = "ASCII";
            header.xml = false;
            header.compress = false;
            break;
        case xmlASCII:
            header.fileFormat = "ASCII";
            header.xml = true;
            header.compress = false;
            break;
        case xmlBinary:
            header.fileFormat = "BINARY";
            header.xml = true;
            header.compress = false;
            break;
        case xmlCompressed:
            header.fileFormat = "BINARY";
            header.xml = true;
            header.compress = true;
            break;
        default:
            header.fileFormat = "BINARY";
            header.xml = false;
            header.compress = false;
            break;
    }

	if ( isLittleEndian() ) header.byteOrder = "LittleEndian";
	else header.byteOrder = "BigEndian";
    if (header.dataType == "SCALARS")
	{
			header.dataSize = 1 * header.pointData;			
	}
	else if (header.dataType == "VECTORS") header.dataSize = 3 * header.pointData;
	
    header.dataCount = 0;
	return header;	
}

void writeVtkFloat(ostream &vtkFile,stVtkHeader& header, float data)
{
	// Verify if the format of the file is ASCII or BINARY and read the data accordingly
	if (header.fileFormat == "ASCII")
		vtkFile << data << " ";
    else if ( (header.fileFormat == "BINARY") && (header.xml) )
    {
        float* floatData = (float*) ( header.dataBuffer + header.dataOffset );
        floatData[header.dataCount++] = data;
    }
    else
    {
		vtkFloatSwap( &data );
		vtkFile.write( (char*) (&data), sizeof(float) );
	}
}

void writeVtkShort(ostream &vtkFile,stVtkHeader& header, short data)
{
	// Verify if the format of the file is ASCII or BINARY and read the data accordingly
	if (header.fileFormat == "ASCII")
		vtkFile << data << " ";
    else if ( (header.fileFormat == "BINARY") && (header.xml) )
    {
        short* shortData = (short*) ( header.dataBuffer + header.dataOffset );
        shortData[header.dataCount++] = data;
    }
    else
	{
		vtkShortSwap( &data );
		vtkFile.write( (char*) (&data), sizeof(short) );
	}
}

bool isLittleEndian()
{
    short int number = 0x1;
    char *numPtr = (char*)&number;
    return (numPtr[0] == 1);
}

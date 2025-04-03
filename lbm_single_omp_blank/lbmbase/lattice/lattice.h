/**
 *  @file   lattice.h
 *  @author Diogo Nardelli Siebert
 *  @brief  This is dummy header, used only to set some macros and choose the lattice
 */

#ifndef __LATTICE_H_INCLUDED__   // if datain.h hasn't been included yet...
#define __LATTICE_H_INCLUDED__

#include "d3q15.h"
#include "d3q19.h"

/* The lattice used is choosen by the namespace set below */
using namespace d3q19;

#define mapMemory(n) (n)*NUM_OF_VEL
#define mapMemoryDir(n,i) (n)*NUM_OF_VEL + (i)

#endif

/**
 *  @file   streaming.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Header of functions necesssary for the streaming step.
 */

#ifndef __STREAMING_H_INCLUDED__   // if header hasn't been included yet...
#define __STREAMING_H_INCLUDED__

#include "../geometry.h"
#include "stdint.h"

struct stStreaming
{
    int64_t cnt;     /*!< Total number of memory positions that will be swap */
    unsigned int *buffer;     /*!< Total number of memory positions that will be swap */
};

/**
 *  @brief  Records the positions that must be swap in the distribution function array during the streaming (for the swap streaming method)
 *  @param geo  A reference to the stGeometry struct containing the information about the medium/domain.
 *  @return A stStreaming struct with the memory positions pairs that must be swap during the second step of the swap streaming method.
 *
 *  It performs a pre-processing recording the memory position that will be exchange during the streaming process using the swap method [1].
 *
 *  [1] Mattila, Keijo, et al. "An efficient swap algorithm for the lattice Boltzmann method." Computer Physics Communications 176.3 (2007): 200-210.
 *
 *  @image html swap.png
 **/

stStreaming defineStreaming(stGeometry& geo);

/**
 *  @brief  Performs the second step of the streaming using the swap method.
 *  @param  streamInfo      A stStreaming struct with the memory positions pairs that must be swapped during the second step of the swap streaming method (generated by defineStreaming function).
 *  @param  iniN            Pointer to the place where the distributions functions are stored in memory.
 *
 *  Performs the second step of the streaming using the swap method. This functions that  the first step was either performed together with collision (faster) or with the preStreaming function.
 *
 *  [1] Mattila, Keijo, et al. "An efficient swap algorithm for the lattice Boltzmann method." Computer Physics Communications 176.3 (2007): 200-210.
 *
 *  @image html swap.png
 */

void streaming(stStreaming &streamInfo, double* iniN);

/**
 *  @brief  Performs the first step of the streaming for the swap method.
 *  @param  iniN                Pointer to the place where the distributions functions are stored in memory.
 *  @param  numberOfPoints      Number of fluid sites in the domain. Note that the size of the iniN array will be numberOfPoints*NUM_OF_VEL
 *
 *  It performs the swap of pairs of distribuntion functions in the same site. It is the first step of swap streaming method. More eficient alternative to this function is to perform this inside
 *  the collision function (functions that have the sufix Pre, ex: collisionBGKPre)
 *
 *  [1] Mattila, Keijo, et al. "An efficient swap algorithm for the lattice Boltzmann method." Computer Physics Communications 176.3 (2007): 200-210.
 *
 *  @image html swap.png
 */

void preStreaming(double* iniN, int numberOfPoints);

#endif

/**
 *  @file   guo.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Declaration of constant and functions for the lattice D3Q19
 *
 *  This file contains explicite formulas for the force term proposed by Guo, Zheng and Shi.
 *
 */

#ifndef __GUO_H_INCLUDED__
#define __GUO_H_INCLUDED__

namespace d3q19
{

namespace guoForce
{
	
/**
 * @brief  Computes the force ther proposed by G
 *
 * @param  F             Point to the array which will stored the force term for each direction in the lattice.
 * @param  Fx            Referece for the value of the x componente of the force
 * @param  Fy            Referece for the value of the y componente of the force
 * @param  Fz            Referece for the value of the z componente of the force
 * @param  vx            macroscopic velocity in the x direction
 * @param  vy            macroscopic velocity in the y direction
 * @param  vz            macroscopic velocity in the z direction
 * @param  tau           The standart relaxation time which is related with the viscosity
 *
 * Computes the inverse of the relaxation time for each of the 19 moments used in the MRT for D3Q19 lattice. For details see [1].
 *
 * [1]  Guo, Zhaoli, Chuguang Zheng, and Baochang Shi. 2002. “Discrete Lattice Effects on the Forcing Term in the Lattice Boltzmann Method.” Physical Review E - Statistical Physics, Plasmas, Fluids, and Related Interdisciplinary Topics 65 (4): 6.
 * https://doi.org/10.1103/PhysRevE.65.046308.
 *
 **/

void forceTermSym(double *F, double& Fx, double&Fy, double&Fz, double& vx, double& vy, double& vz, double& tau);

void forceTermAnt(double *F, double& Fx, double&Fy, double&Fz, double& vx, double& vy, double& vz, double& tau);

void forceTerm(double *F, double& Fx, double&Fy, double&Fz, double& vx, double& vy, double& vz, double& tau);

void calculateMacroEquilibrium( double *f, double& vx, double& vy, double& vz, double& rho , double& Fx, double &Fy, double& Fz);

void calculateMacroOutput( double *f, double& vx, double& vy, double& vz, double& rho, double& Fx, double &Fy, double& Fz);

}

}
#endif

/**
 *  @file   d3q15.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Functions implementations for the D3Q15 lattice
 *
 *  Implementation of functions to compute the discrite equilibrium  distribuition and the macroscopic variables for the D3Q15 lattice
 */

#include "d3q15.h"

namespace d3q15
{
void equilibriumDistribution ( double vx, double vy, double vz, double rho, double* feq)
{

	double cu; // Stores the value of dot product between the fluid velocity (v) and the lattice velocity vector (c_i)
	double rho1 = w1 * rho;  // Stores the weigth times the mass density to avoid recalculation for all directions that shares the same weight
	double rho2 = w2 * rho ;
	// Explicite computes the symmetric  equilibrium  distribution without multiplying by the weight times the mass density
	// since this partial values will be reused for optimization. This multiplication is performed in the final part of the function.

	double a2vx = 3.0 * vx;
	double a2vy = 3.0 * vy;
	double a2vz = 3.0 * vz;

	double vx2 = vx*vx;
	double vy2 = vy*vy;
	double vz2 = vz*vz;

	feq[0] =  ( 1.0 - 1.5 * ( vx2 + vy2 + vz2 ) );

	feq[1] = feq[0] + 4.5 * vx2;
	feq[2] = rho1 * (feq[1] - a2vx);
	feq[1] = rho1 * (feq[1] + a2vx);

	feq[3] = feq[0] + 4.5 * vy2;
	feq[4] = rho1 * (feq[3] - a2vy);
	feq[3] = rho1 * (feq[3] + a2vy);

	feq[5] = feq[0] + 4.5 * vz2;
	feq[6] = rho1 * (feq[5] - a2vz);
	feq[5] = rho1 * (feq[5] + a2vz);

// Stores the weigth times the mass density to avoid recalculation for all directions that shares the same weight

	cu = + a2vx + a2vy + a2vz;
	feq[7] = feq[0] + 0.5 * cu * cu;
	feq[8] = rho2* (feq[7] -  cu) ;
	feq[7] = rho2* (feq[7] +  cu) ;

	cu = - a2vx + a2vy + a2vz;
	feq[9] = feq[0] + 0.5 * cu * cu;
	feq[10] = rho2* (feq[9] -  cu) ;
	feq[9]  = rho2* (feq[9] +  cu) ;

	cu = + a2vx - a2vy + a2vz;
	feq[11] = feq[0] + 0.5 * cu * cu;
	feq[12] = rho2* (feq[11] - cu) ;
	feq[11] = rho2* (feq[11] + cu) ;

	cu = + a2vx + a2vy - a2vz;
	feq[13] = feq[0] + 0.5 * cu * cu;
	feq[14] = rho2* (feq[13] - cu) ;
	feq[13] = rho2* (feq[13] + cu) ;

	// Multiply the the distribution for the 0 direction by the density times the weight.

	feq[0] = w0 * rho * feq[0];


}

void equilibriumDistributionSym( double vx, double vy, double vz, double rho, double* feq)
{
	// Explicite computes the symmetric equilibrium  distribution for the 0 direction of the lattice
	// whithout multiplying by the weight and density.

	feq[0] =  ( 2.0 - 3.0 * ( vx * vx + vy * vy + vz * vz ) );

	double cu; // Stores the value of dot product between the fluid velocity (v) and the lattice velocity vector (c_i)
	double rho1 = w1 * rho;

	// Explicite computes the symmetric  equilibrium  distribution for the 1,3,5,7,9,11,13 direction of the lattice
	// Even directions do not need to the explicit computed since they can be determined from the odd directions
	// by simmetry (Ex. f_sim[2] = f_sim[1]

	feq[1] = rho1*(feq[0] + 9.0 * vx * vx);
	feq[3] = rho1*(feq[0] + 9.0 * vy * vy);
	feq[5] = rho1*(feq[0] + 9.0 * vz * vz);

	double rho2 = w2 * rho;

	cu = vx + vy + vz;
	feq[7] = rho2*(feq[0] + 9.0 * cu * cu);

	cu = -vx + vy + vz;
	feq[9] = rho2*(feq[0] + 9.0 * cu * cu);

	cu = vx - vy + vz;
	feq[11] = rho2*(feq[0] + 9.0 * cu * cu);

	cu = vx + vy - vz;
	feq[13] = rho2*(feq[0] + 9.0 * cu * cu);

	// Multiply the the distribution for the 0 direction by the density times the weight.

	feq[0] = w0 * rho * feq[0];
}

void equilibriumDistributionAnt( double vx, double vy, double vz, double rho, double* feq)
{
	// Explicite computes the antisymmetric equilibrium  distribution for the 0 direction of the lattice
	// whithout multiplying by the weight and density.

	feq[0] =  0;

	double cu; // Stores the value of dot product between the fluid velocity (v) and the lattice velocity vector (c_i)
	double rho1 = w1 * rho ; // Stores the weigth times the mass density to avoid recalculation for all directions that shares the same weight

	// Explicite computes the antisymmetric  equilibrium  distribution for the 1,3,5,7,9,11,13 direction of the lattice
	// Even directions do not need to the explicit computed since they can be determined from the odd directions
	// by simmetry (Ex. f_sim[2] = - f_sim[1]

	feq[1] = rho1 * ( 6.0 * vx );
	feq[3] = rho1 * ( 6.0 * vy );
	feq[5] = rho1 * ( 6.0 * vz );

	double rho2 = w2 * rho;  // Stores the weigth times the mass density to avoid recalculation for all directions that shares the same weight

	cu = + vx + vy + vz;
	feq[7] = rho2 * 6.0 * cu ;

	cu = - vx + vy + vz;
	feq[9]  =  rho2 * 6.0 * cu ;

	cu = + vx - vy + vz;
	feq[11]  = rho2 * 6.0 * cu ;

	cu = + vx + vy - vz;
	feq[13] =  rho2 * 6.0 * cu ;
}

double calculateMomentumX ( double* f)
{
	return f[1] - f[2] + f[7] - f[8] - f[9] + f[10] + f[11] - f[12] + f[13] - f[14];
}

double calculateMass ( double* f )
{
	double rho = 0.0;

	for ( int i = 0 ; i < NUM_OF_VEL; i++ )
	{
		rho = rho + f[i];
	}

	return rho;
}

void calculateMacro( double *f, double& vx, double& vy, double& vz, double& rho )
{
	// Computes the macroscopic variables using explicit formulas for optimization.

	rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14];
	double oneOverRho = 1.0 / rho;

	double mx =  f[1] - f[2] + f[7] - f[8] - f[9] + f[10] + f[11] - f[12] + f[13] - f[14];
	double my =  f[3] - f[4] + f[7] - f[8] + f[9] - f[10] - f[11] + f[12] + f[13] - f[14];
	double mz =  f[5] - f[6] + f[7] - f[8] + f[9] - f[10] + f[11] - f[12] - f[13] + f[14];

	vx = mx * oneOverRho;
	vy = my * oneOverRho;
	vz = mz * oneOverRho;
}

void calculateMoment2nd( double *f, double* m)
{
	m[0] = f[0] + f[1] + f[10] + f[11] + f[12] + f[13] + f[14] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9];
	m[1] = f[1] + f[10] + f[11] - f[12] + f[13] - f[14] - f[2] + f[7] - f[8] - f[9];
	m[2] = -f[10] - f[11] + f[12] + f[13] - f[14] + f[3] - f[4] + f[7] - f[8] + f[9];
	m[3] = -f[10] + f[11] - f[12] - f[13] + f[14] + f[5] - f[6] + f[7] - f[8] + f[9];
	m[4] = f[1] + f[10] + f[11] + f[12] + f[13] + f[14] + f[2] + f[7] + f[8] + f[9];
	m[5] = f[10] + f[11] + f[12] + f[13] + f[14] + f[3] + f[4] + f[7] + f[8] + f[9];
	m[6] = f[10] + f[11] + f[12] + f[13] + f[14] + f[5] + f[6] + f[7] + f[8] + f[9];
	m[7] = -f[10] - f[11] - f[12] + f[13] + f[14] + f[7] + f[8] - f[9];
	m[8] = -f[10] + f[11] + f[12] - f[13] - f[14] + f[7] + f[8] - f[9];
	m[9] = f[10] - f[11] - f[12] - f[13] - f[14] + f[7] + f[8] + f[9];
}

void reconstructedDistribution2nd( double *f, double* m)
{
	f[0] = ((5*m[0] - 3*(m[4] + m[5] + m[6]))*w0)/2.;
	f[1] = (m[0] + 3*m[1] + (6*m[4] - 3*(m[5] + m[6]))/2.)*w1;
	f[2] = ((2*m[0] - 3*(2*m[1] - 2*m[4] + m[5] + m[6]))*w1)/2.;
	f[3] = (m[0] + 3*m[2] - (3*(m[4] - 2*m[5] + m[6]))/2.)*w1;
	f[4] = ((2*m[0] - 3*(2*m[2] + m[4] - 2*m[5] + m[6]))*w1)/2.;
	f[5] = (m[0] + 3*m[3] - (3*(m[4] + m[5] - 2*m[6]))/2.)*w1;
	f[6] = ((2*m[0] - 3*(2*m[3] + m[4] + m[5] - 2*m[6]))*w1)/2.;
	f[7] = (-2*m[0] + 3*(m[1] + m[2] + m[3] + m[4] + m[5] + m[6] + 3*m[7] + 3*m[8] + 3*m[9]))*w2;
	f[8] = -((2*m[0] + 3*(m[1] + m[2] + m[3] - m[4] - m[5] - m[6] - 3*m[7] - 3*m[8] - 3*m[9]))*w2);
	f[9] = -((2*m[0] + 3*m[1] - 3*m[2] - 3*m[3] - 3*m[4] - 3*m[5] - 3*m[6] + 9*m[7] + 9*m[8] - 9*m[9])*w2);
	f[10] = -((2*m[0] - 3*(m[1] - m[2] - m[3] + m[4] + m[5] + m[6] - 3*m[7] - 3*m[8] + 3*m[9]))*w2);
	f[11] = -((2*m[0] - 3*(m[1] - m[2] + m[3] + m[4] + m[5] + m[6] - 3*m[7] + 3*m[8] - 3*m[9]))*w2);
	f[12] = -((2*m[0] + 3*(m[1] - m[2] + m[3] - m[4] - m[5] - m[6] + 3*m[7] - 3*m[8] + 3*m[9]))*w2);
	f[13] = -((2*m[0] - 3*(m[1] + m[2] - m[3] + m[4] + m[5] + m[6] + 3*m[7] - 3*m[8] - 3*m[9]))*w2);
	f[14] = -((2*m[0] + 3*(m[1] + m[2] - m[3] - m[4] - m[5] - m[6] - 3*m[7] + 3*m[8] + 3*m[9]))*w2);
}

void reconstructedDistributionMRT(double* f,double* m)
{
	f[0] = (3*m[0] - 5*m[4] + 2*m[5])/45.;
	f[1] = (6*m[0] + 9*m[1] - 5*m[4] - m[5] - 9*m[6] + 15*m[9])/90.;
	f[2] = (6*m[0] - 9*m[1] - 5*m[4] - m[5] + 9*m[6] + 15*m[9])/90.;
	f[3] = (6*m[0] + 45*m[10] + 9*m[2] - 5*m[4] - m[5] - 9*m[7] + 15*m[9])/90.;
	f[4] = (6*m[0] + 45*m[10] - 9*m[2] - 5*m[4] - m[5] + 9*m[7] + 15*m[9])/90.;
	f[5] = (6*m[0] - 45*m[10] + 9*m[3] - 5*m[4] - m[5] - 9*m[8] - 30*m[9])/90.;
	f[6] = (6*m[0] - 45*m[10] - 9*m[3] - 5*m[4] - m[5] + 9*m[8] - 30*m[9])/90.;
	f[7] = (24*m[0] + 36*m[1] + 45*m[11] + 45*m[12] + 45*m[13] + 45*m[14] + 36*m[2] + 36*m[3] + 20*m[4] + m[5] + 9*m[6] + 9*m[7] + 9*m[8])/360.;
	f[8] = (24*m[0] - 36*m[1] + 45*m[11] + 45*m[12] + 45*m[13] - 45*m[14] - 36*m[2] - 36*m[3] + 20*m[4] + m[5] - 9*m[6] - 9*m[7] - 9*m[8])/360.;
	f[9] = (24*m[0] - 36*m[1] - 45*m[11] + 45*m[12] - 45*m[13] - 45*m[14] + 36*m[2] + 36*m[3] + 20*m[4] + m[5] - 9*m[6] + 9*m[7] + 9*m[8])/360.;
	f[10] = (24*m[0] + 36*m[1] - 45*m[11] + 45*m[12] - 45*m[13] + 45*m[14] - 36*m[2] - 36*m[3] + 20*m[4] + m[5] + 9*m[6] - 9*m[7] - 9*m[8])/360.;
	f[11] = (24*m[0] + 36*m[1] - 45*m[11] - 45*m[12] + 45*m[13] - 45*m[14] - 36*m[2] + 36*m[3] + 20*m[4] + m[5] + 9*m[6] - 9*m[7] + 9*m[8])/360.;
	f[12] = (24*m[0] - 36*m[1] - 45*m[11] - 45*m[12] + 45*m[13] + 45*m[14] + 36*m[2] - 36*m[3] + 20*m[4] + m[5] - 9*m[6] + 9*m[7] - 9*m[8])/360.;
	f[13] = (24*m[0] + 36*m[1] + 45*m[11] - 45*m[12] - 45*m[13] - 45*m[14] + 36*m[2] - 36*m[3] + 20*m[4] + m[5] + 9*m[6] + 9*m[7] - 9*m[8])/360.;
	f[14] = (24*m[0] - 36*m[1] + 45*m[11] - 45*m[12] - 45*m[13] + 45*m[14] - 36*m[2] + 36*m[3] + 20*m[4] + m[5] - 9*m[6] - 9*m[7] + 9*m[8])/360.;
}

void calculateMomentsMRT(double* f,double* m)
{
	m[0] = f[0] + f[1] + f[10] + f[11] + f[12] + f[13] + f[14] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9];
	m[1] = f[1] + f[10] + f[11] - f[12] + f[13] - f[14] - f[2] + f[7] - f[8] - f[9];
	m[2] = -f[10] - f[11] + f[12] + f[13] - f[14] + f[3] - f[4] + f[7] - f[8] + f[9];
	m[3] = -f[10] + f[11] - f[12] - f[13] + f[14] + f[5] - f[6] + f[7] - f[8] + f[9];
	m[4] = -2*f[0] - f[1] + f[10] + f[11] + f[12] + f[13] + f[14] - f[2] - f[3] - f[4] - f[5] - f[6] + f[7] + f[8] + f[9];
	m[5] = 16*f[0] - 4*f[1] + f[10] + f[11] + f[12] + f[13] + f[14] - 4*f[2] - 4*f[3] - 4*f[4] - 4*f[5] - 4*f[6] + f[7] + f[8] + f[9];
	m[6] = -4*f[1] + f[10] + f[11] - f[12] + f[13] - f[14] + 4*f[2] + f[7] - f[8] - f[9];
	m[7] = -f[10] - f[11] + f[12] + f[13] - f[14] - 4*f[3] + 4*f[4] + f[7] - f[8] + f[9];
	m[8] = -f[10] + f[11] - f[12] - f[13] + f[14] - 4*f[5] + 4*f[6] + f[7] - f[8] + f[9];
	m[9] = 2*f[1] + 2*f[2] - f[3] - f[4] - f[5] - f[6];
	m[10] = -f[1] - f[2] + f[3] + f[4];
	m[11] = -f[10] - f[11] - f[12] + f[13] + f[14] + f[7] + f[8] - f[9];
	m[12] = f[10] - f[11] - f[12] - f[13] - f[14] + f[7] + f[8] + f[9];
	m[13] = -f[10] + f[11] + f[12] - f[13] - f[14] + f[7] + f[8] - f[9];
	m[14] = f[10] - f[11] + f[12] - f[13] + f[14] + f[7] - f[8] - f[9];
}

void equilibriumMomentsMRT(double* meq,double* m)
{
	double rho = m[0];
	double vx = m[1]/m[0];
	double vy = m[2]/m[0];
	double vz = m[3]/m[0];
	double vy2 = vy*vy;
	double vz2 = vz*vz;
	double vv = vx*vx + vy2 + vz2;

	meq[0] = m[0];
	meq[1] = m[1];
	meq[2] = m[2];
	meq[3] = m[3];
	meq[4] = rho*(-1 + vv);
	meq[5] = -rho; // meq[5] = rho*(1 - 5*vv);
	meq[6] = - 2.33333333333333333333333333333333333 * m[1];
	meq[7] = - 2.33333333333333333333333333333333333 * m[2];
	meq[8] = - 2.33333333333333333333333333333333333 * m[3];
	meq[9] = rho*(2*vv - 3*(vy2 + vz2) );
	meq[10] = rho*(-vv + 2*vy2 + vz2);
	meq[11] = rho*vx*vy;
	meq[12] = rho*vy*vz;
	meq[13] = rho*vx*vz;
	meq[14] = 0;
}

double* relaxationMRT(double tau)
{
	double *s = new double[NUM_OF_VEL] ();
	s[4] = 1.6;
	s[5] = 1.2;
	s[6] = 1.6;
	s[7] = 1.6;
	s[8] = 1.6;
	s[9] = 1./tau;
	s[10] = 1./tau;
	s[11] = 1./tau;
	s[12] = 1./tau;
	s[13] = 1./tau;
	s[14] = 1.2;
	return s;
}

void calculateMomentum ( double* f , double* m)
{
    m[0] =  f[1] - f[2] + f[7] - f[8] - f[9] + f[10] + f[11] - f[12] + f[13] - f[14];
    m[1] =  f[3] - f[4] + f[7] - f[8] + f[9] - f[10] - f[11] + f[12] + f[13] - f[14];
    m[2] =  f[5] - f[6] + f[7] - f[8] + f[9] - f[10] + f[11] - f[12] - f[13] + f[14];
}

}

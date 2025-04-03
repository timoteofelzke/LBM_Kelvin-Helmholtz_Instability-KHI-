/**
 *  @file   d3q19.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Implementation of Functions for the D3Q19 lattice
 *
 *  Implementation of functions to compute the discrite equilibrium  distribuition and the macroscopic variables for the D3Q19 lattice
 */

#include "d3q19.h"

namespace d3q19
{
/*	
void equilibriumDistribution (double vx,double vy,double vz,double rho,double* feq)
{
	double a2vx = 3.0 * vx;
	double a2vy = 3.0 * vy;
	double a2vz = 3.0 * vz;

	double vx2 = vx*vx;
	double vy2 = vy*vy;
	double vz2 = vz*vz;

	double rho1 = rho*w1;
	double rho2 = rho*w2;

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

	double cv = (a2vx + a2vy);
	feq[ 7] = feq[ 0] + 0.5*cv*cv;
	feq[ 8] = rho2 * (feq[ 7] - cv);
	feq[ 7] = rho2 * (feq[ 7] + cv);

	cv = (a2vx - a2vy);
	feq[ 9]  = feq[ 0] + 0.5*cv*cv;
	feq[10] = rho2 * (feq[ 9] - cv);
	feq[ 9] = rho2 * (feq[ 9] + cv);

	cv = (a2vx + a2vz);
	feq[11] = feq[ 0] + 0.5*cv*cv;
	feq[12] = rho2 * (feq[11] - cv);
	feq[11] = rho2 * (feq[11] + cv);

	cv = (a2vx - a2vz);
	feq[13] = (feq[0] + 0.5*cv*cv);
	feq[14] = rho2* (feq[13] - cv);
	feq[13] = rho2* (feq[13] + cv);

	cv = (-a2vy - a2vz);
	feq[15] = feq[ 0] + 0.5*cv*cv;
	feq[16] = rho2* (feq[15] - cv);
	feq[15] = rho2* (feq[15] + cv);

	cv = (-a2vy + a2vz);
	feq[17] = feq[ 0] + 0.5*cv*cv;
	feq[18] = rho2 * (feq[17] - cv);
	feq[17] = rho2 * (feq[17] + cv);

	feq[0] = w0 *rho * feq[0];

}
*/
void equilibriumDistribution (double vx, double vy, double vz, double rho, double* feq)
{
    double vx2 = vx * vx;
    double vy2 = vy * vy;
    double vz2 = vz * vz;
    double usq = 1.5 * (vx2 + vy2 + vz2);  // |u|^2 * 3/2

    feq[0] = w0 * rho * (1.0 - usq);

    feq[1] = w1 * rho * (1.0 + 3.0 * vx + 4.5 * vx2 - usq);
    feq[2] = w1 * rho * (1.0 - 3.0 * vx + 4.5 * vx2 - usq);

    feq[3] = w1 * rho * (1.0 + 3.0 * vy + 4.5 * vy2 - usq);
    feq[4] = w1 * rho * (1.0 - 3.0 * vy + 4.5 * vy2 - usq);

    feq[5] = w1 * rho * (1.0 + 3.0 * vz + 4.5 * vz2 - usq);
    feq[6] = w1 * rho * (1.0 - 3.0 * vz + 4.5 * vz2 - usq);

    double vxvy = 3.0 * (vx + vy);
    double vxvy2 = 4.5 * (vx + vy) * (vx + vy);
    feq[7]  = w2 * rho * (1.0 + vxvy + vxvy2 - usq);
    feq[8]  = w2 * rho * (1.0 - vxvy + vxvy2 - usq);

    double vxmvy = 3.0 * (vx - vy);
    double vxmvy2 = 4.5 * (vx - vy) * (vx - vy);
    feq[9]  = w2 * rho * (1.0 + vxmvy + vxmvy2 - usq);
    feq[10] = w2 * rho * (1.0 - vxmvy + vxmvy2 - usq);

    double vxvz = 3.0 * (vx + vz);
    double vxvz2 = 4.5 * (vx + vz) * (vx + vz);
    feq[11] = w2 * rho * (1.0 + vxvz + vxvz2 - usq);
    feq[12] = w2 * rho * (1.0 - vxvz + vxvz2 - usq);

    double vxmvz = 3.0 * (vx - vz);
    double vxmvz2 = 4.5 * (vx - vz) * (vx - vz);
    feq[13] = w2 * rho * (1.0 + vxmvz + vxmvz2 - usq);
    feq[14] = w2 * rho * (1.0 - vxmvz + vxmvz2 - usq);

    double vymvz = 3.0 * (vy - vz);
    double vymvz2 = 4.5 * (vy - vz) * (vy - vz);
    feq[15] = w2 * rho * (1.0 + vymvz + vymvz2 - usq);
    feq[16] = w2 * rho * (1.0 - vymvz + vymvz2 - usq);

    double vyvz = 3.0 * (vy + vz);
    double vyvz2 = 4.5 * (vy + vz) * (vy + vz);
    feq[17] = w2 * rho * (1.0 + vyvz + vyvz2 - usq);
    feq[18] = w2 * rho * (1.0 - vyvz + vyvz2 - usq);
}

void equilibriumDistributionSym( double vx, double vy, double vz, double rho, double* feq)
{
	// Explicite computes the symmetric equilibrium  distribution for the 0 direction of the lattice
	// whithout multiplying by the weight and density.

	feq[0] =  ( 2.0 - 3.0 * ( vx * vx + vy * vy + vz * vz ) );

	double cu; // Stores the value of dot product between the fluid velocity (v) and the lattice velocity vector (c_i)
	double rho1 = w1 * rho;

	// Explicitly computes the symmetric  equilibrium  distribution for the 1,3,5,7,9,11,13,15,17 direction of the lattice
	// Even directions do not need to the explicit computed since they can be determined from the odd directions
	// by simmetry (Ex. f_sim[2] = f_sim[1]

	feq[1] = rho1*(feq[0] + 9.0 * vx * vx);
	feq[3] = rho1*(feq[0] + 9.0 * vy * vy);
	feq[5] = rho1*(feq[0] + 9.0 * vz * vz);

	double rho2 = w2 * rho;

	cu = + vx + vy;
	feq[7] = rho2*(feq[0] + 9.0 * cu * cu);

	cu = + vx - vy;
	feq[9] = rho2*(feq[0] + 9.0 * cu * cu);

	cu = + vx + vz;
	feq[11] = rho2*(feq[0] + 9.0 * cu * cu);

	cu = + vx - vz;
	feq[13] = rho2*(feq[0] + 9.0 * cu * cu);

	cu = - vy - vz;
	feq[15] = rho2*(feq[0] + 9.0 * cu * cu);

	cu =  -vy + vz;
	feq[17] = rho2*(feq[0] + 9.0 * cu * cu);


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

	cu = + vx + vy;
	feq[7] = rho2 * 6.0 * cu ;

	cu = + vx - vy;
	feq[9] = rho2 * 6.0 * cu ;

	cu = + vx + vz;
	feq[11] = rho2 * 6.0 * cu ;

	cu = + vx - vz;
	feq[13] = rho2 * 6.0 * cu ;

	cu = - vy - vz;
	feq[15] = rho2 * 6.0 * cu ;

	cu =  -vy + vz;
	feq[17] = rho2 * 6.0 * cu ;
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

	rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18];

	double oneOverRho = 1.0 / rho;


	double mx =  f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + f[11] - f[12] + f[13] - f[14];
	double my =  f[3] - f[4] + f[7] - f[8] - f[9] + f[10] - f[15] + f[16] - f[17] + f[18];
	double mz =  f[5] - f[6] +f[11] -f[12]- f[13] + f[14] - f[15] + f[16] + f[17] - f[18];

	vx = mx * oneOverRho;
	vy = my * oneOverRho;
	vz = mz * oneOverRho;
}

void calculateMoment2nd( double *f, double* m)
{
	m[0] = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18] ;
	m[1] = f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + f[11] - f[12] + f[13] - f[14];
	m[2] = f[3] - f[4] + f[7] - f[8] - f[9] + f[10] - f[15] + f[16] - f[17] + f[18];
	m[3] = f[5] - f[6] + f[11] -f[12]- f[13]+ f[14] - f[15] + f[16] + f[17] - f[18];
	m[4] = f[1] + f[2] + f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14];
	m[5] = f[3] + f[4] + f[7] + f[8] + f[9] + f[10] + f[15] + f[16] + f[17] + f[18];
	m[6] = f[5] + f[6] + f[11]+ f[12] +f[13]+ f[14] + f[15] + f[16] + f[17] + f[18];
	m[7] = f[7] + f[8] - f[9] - f[10];
	m[8] = f[11] + f[12] - f[13] - f[14];
	m[9] = f[15] + f[16] - f[17] - f[18];
}

void reconstructedDistribution2nd( double *f, double* m)
{
	f[0] = (2.5*m[0] - 1.5*m[4] - 1.5*m[5] - 1.5*m[6])*w0;
	f[1] = (m[0] + 3.*m[1]  + 3.*m[4] - 1.5*m[5] - 1.5*m[6])*w1;
	f[2] = (m[0] - 3.*m[1]  + 3.*m[4] - 1.5*m[5] - 1.5*m[6])*w1;
	f[3] = (m[0] + 3.*m[2] - 1.5*m[4] + 3.0*m[5] - 1.5*m[6])*w1;
	f[4] = (m[0] - 3.*m[2] - 1.5*m[4] + 3.0*m[5] - 1.5*m[6])*w1;
	f[5] = (m[0] + 3.*m[3] - 1.5*m[4] - 1.5*m[5] + 3.0*m[6])*w1;
	f[6] = (m[0] - 3.*m[3] - 1.5*m[4] - 1.5*m[5] + 3.0*m[6])*w1;
	f[7] =  (-0.5*m[0] + 3.*m[1] + 3.*m[2] + 3.*m[4] + 3.*m[5] - 1.5*m[6] + 9.*m[7])*w2;
	f[8] =  (-0.5*m[0] - 3.*m[1] - 3.*m[2] + 3.*m[4] + 3.*m[5] - 1.5*m[6] + 9.*m[7])*w2;
	f[9] =  (-0.5*m[0] + 3.*m[1] - 3.*m[2] + 3.*m[4] + 3.*m[5] - 1.5*m[6] - 9.*m[7])*w2;
	f[10] = (-0.5*m[0] - 3.*m[1] + 3.*m[2] + 3.*m[4] + 3.*m[5] - 1.5*m[6] - 9.*m[7])*w2;
	f[11] = (-0.5*m[0] + 3.*m[1] + 3.*m[3] + 3.*m[4] - 1.5*m[5] + 3.*m[6] + 9.*m[8])*w2;
	f[12] = (-0.5*m[0] - 3.*m[1] - 3.*m[3] + 3.*m[4] - 1.5*m[5] + 3.*m[6] + 9.*m[8])*w2;
	f[13] = (-0.5*m[0] + 3.*m[1] - 3.*m[3] + 3.*m[4] - 1.5*m[5] + 3.*m[6] - 9.*m[8])*w2;
	f[14] = (-0.5*m[0] - 3.*m[1] + 3.*m[3] + 3.*m[4] - 1.5*m[5] + 3.*m[6] - 9.*m[8])*w2;
	f[15] = (-0.5*m[0] - 3.*m[2] - 3.*m[3] - 1.5*m[4] + 3.*m[5] + 3.*m[6] + 9.*m[9])*w2;
	f[16] = (-0.5*m[0] + 3.*m[2] + 3.*m[3] - 1.5*m[4] + 3.*m[5] + 3.*m[6] + 9.*m[9])*w2;
	f[17] = (-0.5*m[0] - 3.*m[2] + 3.*m[3] - 1.5*m[4] + 3.*m[5] + 3.*m[6] - 9.*m[9])*w2;
	f[18] = (-0.5*m[0] + 3.*m[2] - 3.*m[3] - 1.5*m[4] + 3.*m[5] + 3.*m[6] - 9.*m[9])*w2;
}

void calculateMomentsMRT(double* f,double* m)
{
	m[0] = f[0] + f[1] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9];
	m[1] = f[1] - f[10] + f[11] - f[12] + f[13] - f[14] - f[2] + f[7] - f[8] + f[9];
	m[2] = f[10] - f[15] + f[16] - f[17] + f[18] + f[3] - f[4] + f[7] - f[8] - f[9];
	m[3] = f[11] - f[12] - f[13] + f[14] - f[15] + f[16] + f[17] - f[18] + f[5] - f[6];
	m[4] = -30*f[0] - 11*f[1] + 8*f[10] + 8*f[11] + 8*f[12] + 8*f[13] + 8*f[14] + 8*f[15] + 8*f[16] + 8*f[17] + 8*f[18] - 11*f[2] - 11*f[3] - 11*f[4] - 11*f[5] - 11*f[6] + 8*f[7] + 8*f[8] + 8*f[9];
	m[5] = 12*f[0] - 4*f[1] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18] - 4*f[2] - 4*f[3] - 4*f[4] - 4*f[5] - 4*f[6] + f[7] + f[8] + f[9];
	m[6] = -4*f[1] - f[10] + f[11] - f[12] + f[13] - f[14] + 4*f[2] + f[7] - f[8] + f[9];
	m[7] = f[10] - f[15] + f[16] - f[17] + f[18] - 4*f[3] + 4*f[4] + f[7] - f[8] - f[9];
	m[8] = f[11] - f[12] - f[13] + f[14] - f[15] + f[16] + f[17] - f[18] - 4*f[5] + 4*f[6];
	m[9] = 2*f[1] + f[10] + f[11] + f[12] + f[13] + f[14] - 2*f[15] - 2*f[16] - 2*f[17] - 2*f[18] + 2*f[2] - f[3] - f[4] - f[5] - f[6] + f[7] + f[8] + f[9];
	m[10] = -4*f[1] + f[10] + f[11] + f[12] + f[13] + f[14] - 2*f[15] - 2*f[16] - 2*f[17] - 2*f[18] - 4*f[2] + 2*f[3] + 2*f[4] + 2*f[5] + 2*f[6] + f[7] + f[8] + f[9];
	m[11] = f[10] - f[11] - f[12] - f[13] - f[14] + f[3] + f[4] - f[5] - f[6] + f[7] + f[8] + f[9];
	m[12] = f[10] - f[11] - f[12] - f[13] - f[14] - 2*f[3] - 2*f[4] + 2*f[5] + 2*f[6] + f[7] + f[8] + f[9];
	m[13] = -f[10] + f[7] + f[8] - f[9];
	m[14] = f[15] + f[16] - f[17] - f[18];
	m[15] = f[11] + f[12] - f[13] - f[14];
	m[16] = -f[10] - f[11] + f[12] - f[13] + f[14] + f[7] - f[8] + f[9];
	m[17] = -f[10] - f[15] + f[16] - f[17] + f[18] - f[7] + f[8] + f[9];
	m[18] = f[11] - f[12] - f[13] + f[14] + f[15] - f[16] - f[17] + f[18];
}

void reconstructedDistributionMRT(double* f,double* m)
{
	f[0] = (21*m[0] - 5*m[4] + 19*m[5])/399.;
	f[1] = (630*m[0] + 1197*m[1] - 665*m[10] - 55*m[4] - 190*m[5] - 1197*m[6] + 665*m[9])/11970.;
	f[2] = (630*m[0] - 1197*m[1] - 665*m[10] - 55*m[4] - 190*m[5] + 1197*m[6] + 665*m[9])/11970.;
	f[3] = m[0]/19. + m[10]/36. + m[11]/12. - m[12]/12. + m[2]/10. - (11*m[4])/2394. - m[5]/63. - m[7]/10. - m[9]/36.;
	f[4] = m[0]/19. + m[10]/36. + m[11]/12. - m[12]/12. - m[2]/10. - (11*m[4])/2394. - m[5]/63. + m[7]/10. - m[9]/36.;
	f[5] = m[0]/19. + m[10]/36. - m[11]/12. + m[12]/12. + m[3]/10. - (11*m[4])/2394. - m[5]/63. - m[8]/10. - m[9]/36.;
	f[6] = m[0]/19. + m[10]/36. - m[11]/12. + m[12]/12. - m[3]/10. - (11*m[4])/2394. - m[5]/63. + m[8]/10. - m[9]/36.;
	f[7] = m[0]/19. + m[1]/10. + m[10]/72. + m[11]/12. + m[12]/24. + m[13]/4. + m[16]/8. - m[17]/8. + m[2]/10. + (4*m[4])/1197. + m[5]/252. + m[6]/40. + m[7]/40. + m[9]/36.;
	f[8] = m[0]/19. - m[1]/10. + m[10]/72. + m[11]/12. + m[12]/24. + m[13]/4. - m[16]/8. + m[17]/8. - m[2]/10. + (4*m[4])/1197. + m[5]/252. - m[6]/40. - m[7]/40. + m[9]/36.;
	f[9] = m[0]/19. + m[1]/10. + m[10]/72. + m[11]/12. + m[12]/24. - m[13]/4. + m[16]/8. + m[17]/8. - m[2]/10. + (4*m[4])/1197. + m[5]/252. + m[6]/40. - m[7]/40. + m[9]/36.;
	f[10] = m[0]/19. - m[1]/10. + m[10]/72. + m[11]/12. + m[12]/24. - m[13]/4. - m[16]/8. - m[17]/8. + m[2]/10. + (4*m[4])/1197. + m[5]/252. - m[6]/40. + m[7]/40. + m[9]/36.;
	f[11] = m[0]/19. + m[1]/10. + m[10]/72. - m[11]/12. - m[12]/24. + m[15]/4. - m[16]/8. + m[18]/8. + m[3]/10. + (4*m[4])/1197. + m[5]/252. + m[6]/40. + m[8]/40. + m[9]/36.;
	f[12] = m[0]/19. - m[1]/10. + m[10]/72. - m[11]/12. - m[12]/24. + m[15]/4. + m[16]/8. - m[18]/8. - m[3]/10. + (4*m[4])/1197. + m[5]/252. - m[6]/40. - m[8]/40. + m[9]/36.;
	f[13] = m[0]/19. + m[1]/10. + m[10]/72. - m[11]/12. - m[12]/24. - m[15]/4. - m[16]/8. - m[18]/8. - m[3]/10. + (4*m[4])/1197. + m[5]/252. + m[6]/40. - m[8]/40. + m[9]/36.;
	f[14] = m[0]/19. - m[1]/10. + m[10]/72. - m[11]/12. - m[12]/24. - m[15]/4. + m[16]/8. + m[18]/8. + m[3]/10. + (4*m[4])/1197. + m[5]/252. - m[6]/40. + m[8]/40. + m[9]/36.;
	f[15] = m[0]/19. - m[10]/36. + m[14]/4. - m[17]/8. + m[18]/8. - m[2]/10. - m[3]/10. + (4*m[4])/1197. + m[5]/252. - m[7]/40. - m[8]/40. - m[9]/18.;
	f[16] = m[0]/19. - m[10]/36. + m[14]/4. + m[17]/8. - m[18]/8. + m[2]/10. + m[3]/10. + (4*m[4])/1197. + m[5]/252. + m[7]/40. + m[8]/40. - m[9]/18.;
	f[17] = m[0]/19. - m[10]/36. - m[14]/4. - m[17]/8. - m[18]/8. - m[2]/10. + m[3]/10. + (4*m[4])/1197. + m[5]/252. - m[7]/40. + m[8]/40. - m[9]/18.;
	f[18] = m[0]/19. - m[10]/36. - m[14]/4. + m[17]/8. + m[18]/8. + m[2]/10. - m[3]/10. + (4*m[4])/1197. + m[5]/252. + m[7]/40. - m[8]/40. - m[9]/18.;
}

void equilibriumMomentsMRT(double* meq,double* m)
{
	double rho = m[0];
	double vx  = m[1]/m[0];
	double vy  = m[2]/m[0];
	double vz  = m[3]/m[0];
//	double vx2 = vx*vx;
	double vy2 = vy*vy;
	double vz2 = vz*vz;
	double vv = vx*vx + vy2 + vz2;
	meq[0] = rho;
	meq[1] = rho*vx;
	meq[2] = rho*vy;
	meq[3] = rho*vz;
	meq[4] = rho*(-11 + 19*vv);
	meq[5] = -475./63. * rho*vv ; // rho*(3 - (11*vv)/2.);
	meq[6] = (-2*rho*vx)/3.;
	meq[7] = (-2*rho*vy)/3.;
	meq[8] = (-2*rho*vz)/3.;
	meq[9] = rho*(2*vv - 3*(vy2 + vz2));
	meq[10] = 0.;// (rho*(-2*vv + 3*(vy2 + vz2)))/2.;
	meq[11] = rho*(vy2 - vz2);
	meq[12] = 0.; // (rho*(-vy2 + vz2))/2.;
	meq[13] = rho*vx*vy;
	meq[14] = rho*vy*vz;
	meq[15] = rho*vx*vz;
	meq[16] = 0;
	meq[17] = 0;
	meq[18] = 0;
}

double* relaxationMRT(double tau)
{
	double *s = new double[NUM_OF_VEL] ();
	s[4] =  1.19;
	s[5] =  1.4;
	s[6] =  1.2;
	s[7] =  1.2;
	s[8] =  1.2;
	s[9] =  1./tau;
	s[10] = 1.4;
	s[11] = 1./tau;
	s[12] = 1.4;
	s[13] = 1./tau;
	s[14] = 1./tau;
	s[15] = 1./tau;
	s[16] = 1.98;
	s[17] = 1.98;
	s[18] = 1.98;
	return s;
}

void calculateMomentum ( double* f , double* m)
{
    m[0] = f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + f[11] - f[12] + f[13] - f[14];
    m[1] = f[3] - f[4] + f[7] - f[8] - f[9] + f[10] - f[15] + f[16] - f[17] + f[18];
    m[2] = f[5] - f[6] +f[11] -f[12]- f[13] + f[14] - f[15] + f[16] + f[17] - f[18];
}

}

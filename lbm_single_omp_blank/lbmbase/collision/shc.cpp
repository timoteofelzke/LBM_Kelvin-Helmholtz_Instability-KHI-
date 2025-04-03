/**
 *  @file   shc.cpp
 *  @author Diogo Nardelli Siebert
 *
 *  Declaration of  functions needed in the collision step.
 */

#include <iostream>
#include <string>
#include <iomanip>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <cmath>

#include "../lattice/lattice.h"
#include "../force/force.h"

#include "shc.h"

using namespace  std;

namespace shcTrt
{

void generateWr(double *ini_R, double* ini_B, double* ini_W,int numberOfPoints)
{
    #pragma omp parallel for
    for (int n = 0; n < numberOfPoints ; n++)
    {
        double* R = ini_R + n * NUM_OF_VEL;
        double* B = ini_B + n * NUM_OF_VEL;
        double* W = ini_W + n * NUM_OF_VEL;

        double rhoR = calculateMass(R);
        double rhoB = calculateMass(B);

        double wr = rhoR / (rhoR + rhoB);

        for (int i = 0; i < NUM_OF_VEL ; i++ )
        {
            W[i] = wr;
        }
    }
}

stCollision collisionParameters(std::string filename)
{
    using boost::property_tree::ptree;
    stCollision prm;

    ptree pt;
    read_ini(filename.c_str(),pt);

    prm.tauSymRed      =  pt.get<double>("SHCTRT.tau_symmetric_red");
    prm.tauSymBlue     =  pt.get<double>("SHCTRT.tau_symmetric_blue");

    prm.tauAntRed  = pt.get<double>("SHCTRT.tau_antisymmetric_red",  ( 8.0 -1.0 / prm.tauSymRed ) / ( 8.0 * (  2.0 - 1.0 / prm.tauSymRed ) ) );
    prm.tauAntBlue = pt.get<double>("SHCTRT.tau_antisymmetric_blue", ( 8.0 -1.0 / prm.tauSymBlue )/ ( 8.0 * (  2.0 - 1.0 / prm.tauSymBlue) ) );

    prm.alpha  = pt.get<double>("SHCTRT.alpha" );
    prm.beta   = pt.get<double>("SHCTRT.beta" );

    prm.wrWall   = pt.get<double>("SHCTRT.wall_concentration" );

    prm.coeffRed  = prm.alpha * one_over_c_s2 / prm.tauSymRed;
    prm.coeffBlue = prm.alpha * one_over_c_s2 / prm.tauSymBlue;

    prm.alphaSymRed = 0.5 / prm.tauSymRed;
    prm.alphaAntRed = 0.5 / prm.tauAntRed;

    prm.alphaSymBlue = 0.5 / prm.tauSymBlue;
    prm.alphaAntBlue = 0.5 / prm.tauAntBlue;

    prm.kinematicViscosityRed  = c_s2 * (prm.tauSymRed - 0.5);
    prm.kinematicViscosityBlue = c_s2 * (prm.tauSymBlue - 0.5);

    prm.surfaceTension = 2 * c_s2 * prm.alpha;

    return prm;
}

void reportCollision(stCollision& prm)
{
    cout << endl;
    cout << "Collision Model: SHC/TRT";
    cout << endl;    
    cout << left << setw(40) << "Relaxation Time Red (Sym):";
    cout << right << setw(30) << prm.tauSymRed << endl;
    cout << left << setw(40) << "Relaxation Time Red (Ant):";
    cout << right << setw(30) << prm.tauAntRed << endl;
    cout << left << setw(40) << "Relaxation Time Blue (Sym):";
    cout << right << setw(30) << prm.tauSymBlue << endl;
    cout << left << setw(40) << "Relaxation Time Blue (Ant):";
    cout << right << setw(30) << prm.tauAntBlue << endl;
    cout << left << setw(40) << "Alpha:";
    cout << right << setw(30) << prm.alpha << endl;
    cout << left << setw(40) << "Beta:";
    cout << right << setw(30) << prm.beta << endl;

    cout << left << setw(40) << "Red concentration at the wall";
    cout << right << setw(30) << prm.wrWall << endl;

    cout << left << setw(40) << "Viscosity Red:";
    cout << right << setw(30) << prm.kinematicViscosityRed << endl;
    cout << left << setw(40) << "Viscosity Blue:";
    cout << right << setw(30) << prm.kinematicViscosityBlue << endl;
    cout << left << setw(40) << "Surface Tension:";
    cout << right << setw(30) <<"Not Implemented" << endl;
    cout << left << setw(40) << "Interface length";
    cout << right << setw(30) <<"Not Implemented" << endl;

    cout << endl;
}

void collisionSwap(double* ini_R, double* ini_B, double* ini_W,  stForce& forceRed, stForce& forceBlue,  stCollision& prm, int numberOfPoints)
{
	#pragma omp parallel for
	for (int n = 0; n < numberOfPoints; n++)
	{	
		double* R = ini_R + n * NUM_OF_VEL;
		double* B = ini_B + n * NUM_OF_VEL;
		double* W = ini_W + n * NUM_OF_VEL;

        double symTau , antTau, coeff, alpha_sim, alpha_ant;

		if (W[0] >= 0.5)
		{
            alpha_sim = prm.alphaSymRed;
            alpha_ant = prm.alphaAntRed;
            coeff  = prm.coeffRed;
            symTau = prm.tauSymRed;
            antTau = prm.tauAntRed;
		}
		else
		{
            alpha_sim = prm.alphaSymBlue;
            alpha_ant = prm.alphaAntBlue;
            coeff  = prm.coeffBlue;
            symTau = prm.tauSymBlue;
            antTau = prm.tauAntBlue;
		}
		
		double grad[3];		
		
		// Compute the Red Concentration Gradient;
		getGradient(grad,W);	
		
		// Computes the norm of the gradient and the normal vector	
		double gradNorm = sqrt( grad[0]*grad[0] + grad[1]*grad[1] + grad[2] * grad[2] );

		double normal[3] = { 0.0, 0.0 , 0.0} ;
		
		if (gradNorm > 1E-10)
		{
			normal[0] = grad[0] / gradNorm;
			normal[1] = grad[1] / gradNorm;
			normal[2] = grad[2] / gradNorm;
		}
		
        double Fshc[ NUM_OF_VEL];
		double f[ NUM_OF_VEL];
		
		for (int i = 0; i < NUM_OF_VEL ; i++ ) f[i] = R[i] + B[i];

		// Computes the symmetric and antisymmetric distribution functions.
		double vx, vy, vz, rho;

        double wB = 1.0 - W[0];

        double* Fred  =   forceRed.constant ?  forceRed.constantForce :  forceRed.forceField + n * NUM_OF_DIM;
        double* Fblue =  forceBlue.constant ? forceBlue.constantForce : forceBlue.forceField + n * NUM_OF_DIM;

        double Fblind[3] = {  W[0]*Fred[0] + wB*Fblue[0] ,W[0]*Fred[1] + wB*Fblue[1]  , W[0]*Fred[2] + wB*Fblue[2]  };

        calculateMacroEquilibrium( f, vx, vy, vz, rho, Fblind[0] ,  Fblind[1]   , Fblind[2] );

        double f_eq_sim[NUM_OF_VEL];
		double f_eq_ant[NUM_OF_VEL];
 
        equilibriumDistributionSym( vx , vy  , vz , rho, f_eq_sim);
        equilibriumDistributionAnt( vx , vy  , vz , rho, f_eq_ant);

		double trt_sim = alpha_sim *  ( f_eq_sim[0]- 2.0 * f[0]);
		double trt_ant = 0;

        shcForce(Fshc, normal ,rho,  (coeff*gradNorm)  );

        double F_sim[NUM_OF_VEL];
        double F_ant[NUM_OF_VEL];

        forceTermSym(F_sim, Fblind[0] ,  Fblind[1] ,  Fblind[2]  ,vx,vy,vz,symTau);
        forceTermAnt(F_ant, Fblind[0] ,  Fblind[1] ,  Fblind[2]  ,vx,vy,vz,antTau);

        f[0] = f[0] + trt_sim + Fshc[0] + F_sim[0] ;
;
		for (int i = 1; i< NUM_OF_VEL; i = i + 2)
		{
			trt_sim = alpha_sim *  ( f_eq_sim[i]   - f[i] - f[i+1] );
			trt_ant = alpha_ant *  ( f_eq_ant[i]   - f[i] + f[i+1] );
            f[i] = f[i] + trt_sim + trt_ant + Fshc[i] + F_sim[i] + F_ant[i];
            f[i+1] = f[i+1] + trt_sim - trt_ant + Fshc[i] + F_sim[i] - F_ant[i];
		}
		
        shcRecoloring(f, R, B, normal, rho, W[0], prm.beta);
	}	
}

}

namespace d3q19
{
	
void getGradient(double* grad, double* f)
{
	grad[0] = one_over_c_s2 * ( w1*(-f[1]+f[2]) + w2*(-f[7]+f[8]-f[9]+f[10]-f[11]+f[12]-f[13]+f[14]) );
	grad[1] = one_over_c_s2 * ( w1*(-f[3]+f[4]) + w2*(-f[7]+f[8]+f[9]-f[10]+f[15]-f[16]+f[17]-f[18]) );
	grad[2] = one_over_c_s2 * ( w1*(-f[5]+f[6]) + w2*(-f[11]+f[12]+f[13]-f[14]+f[15]-f[16]-f[17]+f[18]) );
}

void shcRecoloring(double *f, double *R, double *B, double* normal, double rho, double w_r, double beta)
{
	double w_b = 1.0 - w_r;
	double coeff = beta * rho * w_r * w_b;
	
	R[0] = w_r * f[0];
	B[0] = w_b * f[0];
	
	double df = w1 * coeff * normal[0];
	R[2] = w_r * f[1] + df;
	B[2] = w_b * f[1] - df;
	R[1] = w_r * f[2] - df;
	B[1] = w_b * f[2] + df;
	
	df = w1 * coeff * normal[1];
	R[4] = w_r * f[3] + df;
	B[4] = w_b * f[3] - df;
	R[3] = w_r * f[4] - df;
	B[3] = w_b * f[4] + df;

	df = w1 * coeff * normal[2];
	R[6] = w_r * f[5] + df;
	B[6] = w_b * f[5] - df;
	R[5] = w_r * f[6] - df;
	B[5] = w_b * f[6] + df;
	
	df = w2 * coeff * (normal[0] + normal[1]);
	R[8] = w_r * f[7] + df;
	B[8] = w_b * f[7] - df;
	R[7] = w_r * f[8] - df;
	B[7] = w_b * f[8] + df;

	df = w2 * coeff * (normal[0] - normal[1]);
	R[10] = w_r * f[9] + df;
	B[10] = w_b * f[9] - df;
	R[9] = w_r * f[10] - df;
	B[9] = w_b * f[10] + df;

	df = w2 * coeff * (normal[0] + normal[2]);
	R[12] = w_r * f[11] + df;
	B[12] = w_b * f[11] - df;
	R[11] = w_r * f[12] - df;
	B[11] = w_b * f[12] + df;
	
	df = w2 * coeff * (normal[0] - normal[2]);
	R[14] = w_r * f[13] + df;
	B[14] = w_b * f[13] - df;
	R[13] = w_r * f[14] - df;
	B[13] = w_b * f[14] + df;
	
	df = w2 * coeff * (-normal[1] - normal[2]);
	R[16] = w_r * f[15] + df;
	B[16] = w_b * f[15] - df;
	R[15] = w_r * f[16] - df;
	B[15] = w_b * f[16] + df;
	
	df = w2 * coeff * (-normal[1] + normal[2]);
	R[18] = w_r * f[17] + df;
	B[18] = w_b * f[17] - df;
	R[17] = w_r * f[18] - df;
	B[17] = w_b * f[18] + df;
}	
	
void shcForce(double* F, double* normal,double& rho,  double coeff )
{
	double nc = 0;
	F[0] = rho * w0 * coeff  * (2./3.) ;

	double rho1 = rho * w1 * coeff ;
	nc = normal[0];
	F[1] = rho1 * (nc*nc - 1./3.);

	nc = normal[1];
	F[3] = rho1 *  (nc*nc - 1./3.);

	nc = normal[2];
	F[5] = rho1 *  (nc*nc - 1./3.);
	
	double rho2 = rho * w2 * coeff ;

	nc = normal[0] + normal[1];
	F[7] = rho2 * (nc*nc - 4./3.);

	nc = normal[0] - normal[1];
	F[9] = rho2 * (nc*nc - 4./3.);

	nc = normal[0] + normal[2];
	F[11] = rho2 * (nc*nc - 4./3.);

	nc = normal[0] - normal[2];
	F[13] = rho2 * (nc*nc - 4./3.);

	nc = -normal[1] - normal[2];
	F[15] = rho2 * (nc*nc - 4./3.);

	nc = -normal[1] + normal[2];
	F[17] = rho2 * (nc*nc - 4./3.);	
}

}


/**
 *  @file   zouhe.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Header of functions to impose boundary conditions (pressure and velocity) using the Zou & He method.
 *
 *  @image html boundary.png
 */

#ifndef __ZOUHE_H_INCLUDED__
#define __ZOUHE_H_INCLUDED__

struct stGeometry;

namespace d3q15
{

/**
 *  Sets a defined density (or Pressure since is an ideal gas) in a given site using the
 *  Zou & He boundary condition for the d3q15 lattice. This function is used when the boundary
 *  is a YZ plane (constant X coordinate) located at the left of the fluid domain, i. e., the first YZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/ovt6ote3h0wtlp5/ZouHe_D3Q15_Pressure.nb?dl=0
 *
 *  [1] ZOU, Q.; HE, X. On pressure and velocity flow boundary conditions and bounceback for the lattice Boltzmann BGK model. v. 9, n. 6, p. 18, 1996.
 *
 *  @brief  Set Fixed Pressure Value on the Left Boundary using Zou & He method for the d3q15 lattice.
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  rho      Value of the density
 *
 */

void zouHePressureLeft(double *f, double rho);

/**
 *  Sets a defined density (or Pressure since is an ideal gas) in a given site using the
 *  Zou & He boundary condition for the d3q15 lattice. This function is used when the boundary
 *  is a YZ plane (constant X coordinate) located at the right of the fluid domain, i. e., the last YZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/ovt6ote3h0wtlp5/ZouHe_D3Q15_Pressure.nb?dl=0
 *
 *  [1] ZOU, Q.; HE, X. On pressure and velocity flow boundary conditions and bounceback for the lattice Boltzmann BGK model. v. 9, n. 6, p. 18, 1996.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  rho      Value of the density
 */

void zouHePressureRight(double *f, double rho);

/**
 *  Sets a defined density (or Pressure since is an ideal gas) in a given site using the
 *  Zou & He boundary condition for the d3q15 lattice. This function is used when the boundary
 *  is a XZ plane (constant Y coordinate) located at the Front of the fluid domain, i. e., the first XZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/ovt6ote3h0wtlp5/ZouHe_D3Q15_Pressure.nb?dl=0
 *
 *  [1] ZOU, Q.; HE, X. On pressure and velocity flow boundary conditions and bounceback for the lattice Boltzmann BGK model. v. 9, n. 6, p. 18, 1996.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  rho      Value of the density
 */

void zouHePressureFront(double *f, double rho);

/**
 *  Sets a defined density (or Pressure since is an ideal gas) in a given site using the
 *  Zou & He boundary condition for the d3q15 lattice. This function is used when the boundary
 *  is a XZ plane (constant Y coordinate) located at the back of the fluid domain, i. e., the last XZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/ovt6ote3h0wtlp5/ZouHe_D3Q15_Pressure.nb?dl=0
 *
 *  [1] ZOU, Q.; HE, X. On pressure and velocity flow boundary conditions and bounceback for the lattice Boltzmann BGK model. v. 9, n. 6, p. 18, 1996.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  rho      Value of the density
 */

void zouHePressureBack(double *f, double rho);

/**
 *  Sets a defined density (or Pressure since is an ideal gas) in a given site using the
 *  Zou & He boundary condition for the d3q15 lattice. This function is used when the boundary
 *  is a XY plane (constant Z coordinate) located at the bottom of the fluid domain, i. e., the first XY plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/ovt6ote3h0wtlp5/ZouHe_D3Q15_Pressure.nb?dl=0
 *
 *  [1] ZOU, Q.; HE, X. On pressure and velocity flow boundary conditions and bounceback for the lattice Boltzmann BGK model. v. 9, n. 6, p. 18, 1996.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  rho      Value of the density
 */

void zouHePressureBottom(double *f, double rho);

/**
 *  Sets a defined density (or Pressure since is an ideal gas) in a given site using the
 *  Zou & He boundary condition for the d3q15 lattice. This function is used when the boundary
 *  is a XY plane (constant Z coordinate) located at the top of the fluid domain, i. e., the last XY plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/ovt6ote3h0wtlp5/ZouHe_D3Q15_Pressure.nb?dl=0
 *
 *  [1] ZOU, Q.; HE, X. On pressure and velocity flow boundary conditions and bounceback for the lattice Boltzmann BGK model. v. 9, n. 6, p. 18, 1996.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  rho      Value of the density
 */

void zouHePressureTop(double *f, double rho);

/**
 *  Sets a defined velocity (vx,vy,vz) in a given site using the Zou & He boundary
 *  condition for the d3q15 lattice. This function is used when the boundary is a YZ
 *  plane (constant X coordinate) located at the left of the fluid domain, i. e., the first YZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/567e8xhzzq6qegz/ZouHe_D3Q15.nb?dl=0
 *
 *  [1] ZOU, Q.; HE, X. On pressure and velocity flow boundary conditions and bounceback for the lattice Boltzmann BGK model. v. 9, n. 6, p. 18, 1996.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  vx       X component of the velocity that will be set at site.
 *  @param  vy       Y component of the velocity that will be set at site.
 *  @param  vz       Z component of the velocity that will be set at site.
 */

void zouHeVelocityLeft(double *f, double vx, double vy, double vz);

/**
 *  Sets a defined velocity (vx,vy,vz) in a given site using the Zou & He boundary
 *  condition for the d3q15 lattice. This function is used when the boundary is a YZ
 *  plane (constant X coordinate) located at the right of the fluid domain, i. e., the last YZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/567e8xhzzq6qegz/ZouHe_D3Q15.nb?dl=0
 *
 *  [1] ZOU, Q.; HE, X. On pressure and velocity flow boundary conditions and bounceback for the lattice Boltzmann BGK model. v. 9, n. 6, p. 18, 1996.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  vx       X component of the velocity that will be set at site.
 *  @param  vy       Y component of the velocity that will be set at site.
 *  @param  vz       Z component of the velocity that will be set at site.
 */

void zouHeVelocityRight(double *f, double vx, double vy, double vz);

/**
 *  Sets a defined velocity (vx,vy,vz) in a given site using the Zou & He boundary
 *  condition for the d3q15 lattice. This function is used when the boundary is a XZ
 *  plane (constant Y coordinate) located at the front of the fluid domain, i. e., the first XZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/567e8xhzzq6qegz/ZouHe_D3Q15.nb?dl=0
 *
 *  [1] ZOU, Q.; HE, X. On pressure and velocity flow boundary conditions and bounceback for the lattice Boltzmann BGK model. v. 9, n. 6, p. 18, 1996.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  vx       X component of the velocity that will be set at site.
 *  @param  vy       Y component of the velocity that will be set at site.
 *  @param  vz       Z component of the velocity that will be set at site.
 */

void zouHeVelocityFront(double *f, double vx, double vy, double vz);

/**
 *  Sets a defined velocity (vx,vy,vz) in a given site using the Zou & He boundary
 *  condition for the d3q15 lattice. This function is used when the boundary is a XZ
 *  plane (constant Y coordinate) located at the back of the fluid domain, i. e., the last XZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/567e8xhzzq6qegz/ZouHe_D3Q15.nb?dl=0
 *
 *  [1] ZOU, Q.; HE, X. On pressure and velocity flow boundary conditions and bounceback for the lattice Boltzmann BGK model. v. 9, n. 6, p. 18, 1996.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  vx       X component of the velocity that will be set at site.
 *  @param  vy       Y component of the velocity that will be set at site.
 *  @param  vz       Z component of the velocity that will be set at site.
 */

void zouHeVelocityBack(double *f, double vx, double vy, double vz);

/**
 *  Sets a defined velocity (vx,vy,vz) in a given site using the Zou & He boundary
 *  condition for the d3q15 lattice. This function is used when the boundary is a XY
 *  plane (constant Z coordinate) located at the bottom of the fluid domain, i. e., the first XY plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/567e8xhzzq6qegz/ZouHe_D3Q15.nb?dl=0
 *
 *  [1] ZOU, Q.; HE, X. On pressure and velocity flow boundary conditions and bounceback for the lattice Boltzmann BGK model. v. 9, n. 6, p. 18, 1996.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  vx       X component of the velocity that will be set at site.
 *  @param  vy       Y component of the velocity that will be set at site.
 *  @param  vz       Z component of the velocity that will be set at site.
 */

void zouHeVelocityBottom(double *f, double vx, double vy, double vz);

/**
 *  Sets a defined velocity (vx,vy,vz) in a given site using the Zou & He boundary
 *  condition for the d3q15 lattice. This function is used when the boundary is a XY
 *  plane (constant Z coordinate) located at the top of the fluid domain, i. e., the last XY plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/567e8xhzzq6qegz/ZouHe_D3Q15.nb?dl=0
 *
 *  [1] ZOU, Q.; HE, X. On pressure and velocity flow boundary conditions and bounceback for the lattice Boltzmann BGK model. v. 9, n. 6, p. 18, 1996.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  vx       X component of the velocity that will be set at site.
 *  @param  vy       Y component of the velocity that will be set at site.
 *  @param  vz       Z component of the velocity that will be set at site.
 */

void zouHeVelocityTop(double *f, double vx, double vy, double vz);
}

namespace d3q19
{

/**
 *  Sets a defined velocity (vx,vy,vz) in a given site using the Zou & He boundary
 *  condition for the d3q19 lattice. This function is used when the boundary is a YZ
 *  plane (constant X coordinate) located at the left of the fluid domain, i. e., the first YZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/nvvcokmdh5mi29n/ZouHe_D3Q19.nb?dl=0
 *
 *  [1] Han, Y. (2015). An Analytical Boundary Condition for D3Q19 Lattice Boltzmann Model.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  vx       X component of the velocity that will be set at site.
 *  @param  vy       Y component of the velocity that will be set at site.
 *  @param  vz       Z component of the velocity that will be set at site.
 */

void zouHeVelocityLeft(double *f, double vx, double vy, double vz);

/**
 *  Sets a defined velocity (vx,vy,vz) in a given site using the Zou & He boundary
 *  condition for the d3q19 lattice. This function is used when the boundary is a YZ
 *  plane (constant X coordinate) located at the right of the fluid domain, i. e., the last YZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/nvvcokmdh5mi29n/ZouHe_D3Q19.nb?dl=0
 *
 *  [1] Han, Y. (2015). An Analytical Boundary Condition for D3Q19 Lattice Boltzmann Model.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  vx       X component of the velocity that will be set at site.
 *  @param  vy       Y component of the velocity that will be set at site.
 *  @param  vz       Z component of the velocity that will be set at site.
 */

void zouHeVelocityRight(double* f, double vx, double vy, double vz);

/**
 *  Sets a defined velocity (vx,vy,vz) in a given site using the Zou & He boundary
 *  condition for the d3q19 lattice. This function is used when the boundary is a XZ
 *  plane (constant Y coordinate) located at the front of the fluid domain, i. e., the first XZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/nvvcokmdh5mi29n/ZouHe_D3Q19.nb?dl=0
 *
 *  [1] Han, Y. (2015). An Analytical Boundary Condition for D3Q19 Lattice Boltzmann Model.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  vx       X component of the velocity that will be set at site.
 *  @param  vy       Y component of the velocity that will be set at site.
 *  @param  vz       Z component of the velocity that will be set at site.
 */

void zouHeVelocityFront(double* f, double vx, double vy, double vz);

/**
 *  Sets a defined velocity (vx,vy,vz) in a given site using the Zou & He boundary
 *  condition for the d3q19 lattice. This function is used when the boundary is a XZ
 *  plane (constant Y coordinate) located at the back of the fluid domain, i. e., the last XZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/nvvcokmdh5mi29n/ZouHe_D3Q19.nb?dl=0
 *
 *  [1] Han, Y. (2015). An Analytical Boundary Condition for D3Q19 Lattice Boltzmann Model.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  vx       X component of the velocity that will be set at site.
 *  @param  vy       Y component of the velocity that will be set at site.
 *  @param  vz       Z component of the velocity that will be set at site.
 */

void zouHeVelocityBack(double* f, double vx, double vy, double vz);

/**
 *  Sets a defined velocity (vx,vy,vz) in a given site using the Zou & He boundary
 *  condition for the d3q19 lattice. This function is used when the boundary is a XY
 *  plane (constant Z coordinate) located at the bottom of the fluid domain, i. e., the first XY plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/nvvcokmdh5mi29n/ZouHe_D3Q19.nb?dl=0
 *
 *  [1] Han, Y. (2015). An Analytical Boundary Condition for D3Q19 Lattice Boltzmann Model.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  vx       X component of the velocity that will be set at site.
 *  @param  vy       Y component of the velocity that will be set at site.
 *  @param  vz       Z component of the velocity that will be set at site.
 */

void zouHeVelocityBottom(double* f, double vx, double vy, double vz);

/**
 *  Sets a defined velocity (vx,vy,vz) in a given site using the Zou & He boundary
 *  condition for the d3q19 lattice. This function is used when the boundary is a XY
 *  plane (constant Z coordinate) located at the top of the fluid domain, i. e., the last XY plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/nvvcokmdh5mi29n/ZouHe_D3Q19.nb?dl=0
 *
 *  [1] Han, Y. (2015). An Analytical Boundary Condition for D3Q19 Lattice Boltzmann Model.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  vx       X component of the velocity that will be set at site.
 *  @param  vy       Y component of the velocity that will be set at site.
 *  @param  vz       Z component of the velocity that will be set at site.
 */

void zouHeVelocityTop(double* f, double vx, double vy, double vz);

/**
 *  Sets a defined density (or Pressure since is an ideal gas) in a given site using the
 *  Zou & He boundary condition for the d3q19 lattice. This function is used when the boundary
 *  is a YZ plane (constant X coordinate) located at the left of the fluid domain, i. e., the first YZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/mtjdl3qzkvgkpgu/ZouHe_D3Q19_Pressure.nb?dl=0
 *
 *  [1] Han, Y. (2015). An Analytical Boundary Condition for D3Q19 Lattice Boltzmann Model.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  rho      Value of the density
 */

void zouHePressureLeft(double *f, double rho);

/**
 *  Sets a defined density (or Pressure since is an ideal gas) in a given site using the
 *  Zou & He boundary condition for the d3q19 lattice. This function is used when the boundary
 *  is a YZ plane (constant X coordinate) located at the right of the fluid domain, i. e., the last YZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/mtjdl3qzkvgkpgu/ZouHe_D3Q19_Pressure.nb?dl=0
 *
 *  [1] Han, Y. (2015). An Analytical Boundary Condition for D3Q19 Lattice Boltzmann Model.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  rho      Value of the density
 */

void zouHePressureRight(double *f, double rho);

/**
 *  Sets a defined density (or Pressure since is an ideal gas) in a given site using the
 *  Zou & He boundary condition for the d3q19 lattice. This function is used when the boundary
 *  is a XZ plane (constant Y coordinate) located at the Front of the fluid domain, i. e., the first XZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/mtjdl3qzkvgkpgu/ZouHe_D3Q19_Pressure.nb?dl=0
 *
 *  [1] Han, Y. (2015). An Analytical Boundary Condition for D3Q19 Lattice Boltzmann Model.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  rho      Value of the density
 */

void zouHePressureFront(double *f, double rho);

/**
 *  Sets a defined density (or Pressure since is an ideal gas) in a given site using the
 *  Zou & He boundary condition for the d3q19 lattice. This function is used when the boundary
 *  is a XZ plane (constant Y coordinate) located at the back of the fluid domain, i. e., the last XZ plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/mtjdl3qzkvgkpgu/ZouHe_D3Q19_Pressure.nb?dl=0
 *
 *  [1] Han, Y. (2015). An Analytical Boundary Condition for D3Q19 Lattice Boltzmann Model.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  rho      Value of the density
 */

void zouHePressureBack(double *f, double rho);

/**
 *  Sets a defined density (or Pressure since is an ideal gas) in a given site using the
 *  Zou & He boundary condition for the d3q19 lattice. This function is used when the boundary
 *  is a XY plane (constant Z coordinate) located at the bottom of the fluid domain, i. e., the first XY plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/mtjdl3qzkvgkpgu/ZouHe_D3Q19_Pressure.nb?dl=0
 *
 *  [1] Han, Y. (2015). An Analytical Boundary Condition for D3Q19 Lattice Boltzmann Model.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  rho      Value of the density
 */

void zouHePressureBottom(double *f, double rho);

/**
 *  Sets a defined density (or Pressure since is an ideal gas) in a given site using the
 *  Zou & He boundary condition for the d3q19 lattice. This function is used when the boundary
 *  is a XY plane (constant Z coordinate) located at the top of the fluid domain, i. e., the last XY plane.
 *
 *  The formulas were calculated by Diogo Siebert using the Mathematica Notebook from the link
 *
 *  https://www.dropbox.com/s/mtjdl3qzkvgkpgu/ZouHe_D3Q19_Pressure.nb?dl=0
 *
 *  [1] Han, Y. (2015). An Analytical Boundary Condition for D3Q19 Lattice Boltzmann Model.
 *
 *  @param  f        Pointer where the distribuition function is stored . The values of these distribution will change to match the boundary condition.
 *  @param  rho      Value of the density
 */

void zouHePressureTop(double *f, double rho);

}

#endif

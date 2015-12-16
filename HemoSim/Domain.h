#ifndef DOMAIN_H
#define DOMAIN_H

#include "Element.h"

typedef struct domain {
  // 1D Model vars
  int     nel;          /**< Number of elements in the domain. */
  int     hispts;       /**< Number of history points in the domain. */
  int     **bifur;      /**< Connectivity with the other 2 domains joining in the bifurcation. */

  char    *bctype;      /**< Boundary condition type ('a', 'A', 'd', 'D', 'T', 'W'...) for A & u at both sides of the domain. */
  double  *bcval;       /**< Numerical value of the boundary condition for A & u at both sides of the domain. */
  char    **bcstring;   /**< Character string with the function introduced as boundary condition for A & u at both sides of the domain. */
  char    *solution;
  Element *U,*A;        /**< Flow and area variables. */
  Element **Uf,**Af;    /**< Momentum and mass fluxes. */
  Element *Ht, *Hd;		  /**< Tube & Discharge Hematocrit of the domain. Added by panqing. 2012-1-31 */

  double *hisptsx; /**< Axial coordinate of each history point "hispts". */
  double *oldAhis, *oldUhis, *oldHdhis;
  double Ahis[3], Uhis[3], Hdhis[3];
  double *visc;  /**< Viscosity of the domain. Added by panqing. 2011-07-29 */
  double *gamma; /**< visco part of the viscoelastic property of the vessel wall */
  
  // Metabolic boundaries
  // added by pq. 2015-10-08
  double BHd;   // Boundary Hd
  double BSO2;  // Boundary SO2 input
  double BJm;   // Boundary Jm
  double BJc;   // Boundary Jc

  // Wall Thickness (for structural adaptation)
  double WallTh;

  // Measured Velocity (for calculating Ev)
  double MesVel;

  // Metabolic parameters
  // For PO2
  double JconM[3];  // 0:JconMin, 1:JconMmid, 2:JconMout
  double SO2[3];    // 0:SO2in, 1:SO2mid, 2:SO2out
  double PO2;
  // For Sm
  double Jmmid;
  double Sm;
  // For Sc
  double Jcmid;
  double Sc;
  // For Sp and Stau
  double Sp;
  double Stau;

  // 0D Model vars
  double ZeroD_R, ZeroD_C, ZeroD_I; // Resistance, Capacitor and Inertia of the domain(segment)
  int nodes[2][4];  /**< start and end nodes of the domain(segment):nodes[0/1][0], and their respective connected segments:nodes[0/1][1,2,3] */

  // Steady state Model vars
  double J; // Conductance
} Domain;

#endif

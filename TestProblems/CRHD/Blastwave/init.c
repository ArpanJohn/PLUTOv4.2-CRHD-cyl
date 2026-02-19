/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

int flag = 0;
/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 *
 *********************************************************************** */
{

  if (flag == 0){
    printf("--- Variable Check (NVAR = %d) ---\n", NVAR);
    
    // Print standard primitive/conservative variables
    printf("RHO: %d\n", RHO);
    printf("VX1: %d, VX2: %d, VX3: %d\n", VX1, VX2, VX3);
    #if PHYSICS == MHD
      printf("BX1: %d, BX2: %d, BX3: %d\n", BX1, BX2, BX3);
    #endif
    printf("PRS: %d\n", PRS);
    printf("PCR: %d\n", PCR);

    if (ECR >= NVAR){
      printf("ERROR: ECR index %d >= NVAR %d\n", ECR, NVAR);
      QUIT_PLUTO(1);
    }
    else {
      printf("ECR index %d safe,  NVAR %d\n", ECR, NVAR);
    }


    // // Loop to print ALL variables by index to ensure none are missed
    // for (int i = 0; i < NVAR; i++) {
    //     printf("Variable Index [%d] Value: %e\n", i, i);
    // }
    printf("----------------------------------\n");
    flag = 1;
  }
  double dr, vol, r;
  double rinj = g_inputParam[Rinj]; // im cgs units;
  vol = 4.0/3.0*CONST_PI*rinj*rinj*rinj;

#if GEOMETRY == CARTESIAN
  r = sqrt(x1*x1 + x2*x2 + x3*x3);
#elif GEOMETRY == CYLINDRICAL
  r = sqrt(x1*x1 + x2*x2);
#elif GEOMETRY  == SPHERICAL
  r = x1; // Spherical geometry
#endif

  v[RHO] = g_inputParam[RHO_AMB]/UNIT_DENSITY;
  EXPAND(v[VX1] = 0.0;, v[VX2] = 0.0;, v[VX3] = 0.0;)
  v[PRS] = (g_inputParam[RHO_AMB]/CONST_mH/0.6)*CONST_kB*g_inputParam[TEMP_AMB]/unitPRS;
  
  #if CR_FLUID != NO
  double fcr = 0.1;
  v[PCR] = fcr*v[PRS]; // or any small value (e.g, 1.e-30);
  #endif


  #if PHYSICS == MHD 
    double beta_bck = 2.22;
    double UNIT_MAG = UNIT_VELOCITY * sqrt(4.0 * CONST_PI * UNIT_DENSITY);
    /* setting up a purely toroidal field with a constant plasma beta equal to plasma beta of injected field */
    double  B_phi = 0;
    double B_mag = sqrt( 8 / beta_bck * CONST_PI * (v[PRS]*unitPRS)) ; /* calculing the amount  of injected toroidal field based on beta */
    B_phi = B_mag / UNIT_MAG;

    // if (flag < 10) {
    //   printf("\n\n UNIT_MAG = %.4e \n B_phi
    
#if GEOMETRY == CYLINDRICAL
    v[iBR] = 0;
    v[iBZ] = 0;
    v[iBPHI] = 0;

    if (!isnan(B_phi)) v[iBPHI]= B_phi;
# else
    v[BX1] = 0;
    v[BX2] = 0;
    v[BX3] = 0;

    if (!isnan(B_phi)) v[BX3]= B_phi;
#endif
  #endif


  if (r <= rinj/UNIT_LENGTH ) {
  //  v[PRS] += (g_gamma - 1.0)*(g_inputParam[E_SN]/unitENERGY)/vol;
   v[PRS] += (g_gamma - 1.0)*(g_inputParam[E_SN]/vol/unitPRS);
  #if CR_FLUID != NO
   v[PCR] += fcr * (g_gamma - 1.0)*(g_inputParam[E_SN]/vol/unitPRS); 
  #endif
  } 
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{

}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;


  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){
      if (d->Vc[RHO][k][j][i] < 1e-10){d->Vc[RHO][k][j][i] = 1e-10;}
      if (d->Vc[PRS][k][j][i] < 1e-10){d->Vc[PRS][k][j][i] = 1e-10;}
  #if CR_FLUID != NO
      if(d->Vc[PCR][k][j][i] < 1e-10){d->Vc[PCR][k][j][i] = 1e-10;}
  #endif
    }
  }
}



/*@@
  @file      Boundary2.h
  @date      Tue Sep 26 11:50:46 2000
  @author    Gerd Lanfermann
  @desc
             Prototypes for boundary routines
  @enddesc
  @version   $Header$
@@*/

#ifndef _BOUNDARY2_H_
#define _BOUNDARY2_H_

#include "PreSync.h"

#ifdef __cplusplus
extern "C" {
#endif

/* check boundary width and abort if unlikely large (>100 points) */
void BndSanityCheckWidths2(const cGH *GH, CCTK_INT varindex, CCTK_INT dim,
                          const CCTK_INT *boundary_widths, const char *bcname);

/* prototype for routine registered as providing 'None' boundary condition */
void Bndry_None(const cGH *cctkGH, CCTK_INT num_vars, CCTK_INT *var_indices,
                 CCTK_INT *faces, CCTK_INT *widths,
                 CCTK_INT *table_handles);

/* prototype for routine registered as providing 'Scalar' boundary condition */
void Bndry_Scalar(const cGH *cctkGH, CCTK_INT num_vars, CCTK_INT *var_indices,
                   CCTK_INT *faces, CCTK_INT *widths,
                   CCTK_INT *table_handles);

/* prototype for routine registered as providing 'Copy' boundary condition */
void Bndry_Copy(const cGH *cctkGH, CCTK_INT num_vars, CCTK_INT *var_indices,
                 CCTK_INT *faces, CCTK_INT *widths,
                 CCTK_INT *table_handles);

/* prototype for routine registered as providing 'Static' boundary condition */
void Bndry_Static(const cGH *cctkGH, CCTK_INT num_vars, CCTK_INT *var_indices,
                   CCTK_INT *faces, CCTK_INT *widths,
                   CCTK_INT *table_handles);

/* prototype for routine registered as providing 'Radiative' boundary conditions */
void Bndry_Radiative(const cGH *cctkGH, CCTK_INT num_vars, CCTK_INT *var_indices,
                      CCTK_INT *faces, CCTK_INT *widths,
                      CCTK_INT *table_handles);

/* prototype for routine registered as providing 'Robin' boundary condition */
void Bndry_Robin(const cGH *cctkGH, CCTK_INT num_vars, CCTK_INT *var_indices,
                  CCTK_INT *faces, CCTK_INT *widths,
                  CCTK_INT *table_handles);

/* prototype for routine registered as providing 'Flat' boundary condition */
void Bndry_Flat(const cGH *cctkGH, CCTK_INT num_vars, CCTK_INT *var_indices,
                 CCTK_INT *faces, CCTK_INT *widths,
                 CCTK_INT *table_handles);

#ifdef __cplusplus
}
#endif

#endif /* _BOUNDARY2_H_ */

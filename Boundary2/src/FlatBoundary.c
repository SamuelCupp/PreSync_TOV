/*@@
  @file      FlatBoundary.c
  @date      Mon Mar 15 15:09:00 1999
  @author    Gerd Lanfermann
  @desc
             Routines for applying flat boundary conditions
  @enddesc
  @history
  @hdate     Tue 10 Apr 2001
  @hauthor   Thomas Radke
  @hdesc     BC routines generalized for applying to arbitrary CCTK data types
  @endhistory
  @version   $Id$
@@*/

/*#define DEBUG_BOUNDARY*/

#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "util_Table.h"
#include "util_ErrorCodes.h"
#include "cctk_Parameters.h"
#include "cctk_FortranString.h"
#include "Boundary2.h"

static int ApplyBndFlat(const cGH *GH, CCTK_INT stencil_dir,
                        const CCTK_INT *stencil_alldirs,
                        int dir, CCTK_INT faces,
                        int first_var, int num_vars);

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
/*@@
   @routine    BndFlat
   @date       13 Feb 2003
   @author     David Rideout
   @desc
               Top level function which is registered as handling
               the Flat boundary condition
   @enddesc
   @calls      ApplyBndFlat

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        num_vars
   @vdesc      number of variables passed in through var_indices[]
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @var        var_indices
   @vdesc      array of variable indicies to which to apply this boundary
               condition
   @vtype      CCTK_INT *
   @vio        in
   @endvar
   @var        faces
   @vdesc      array of set of faces to which to apply the bc
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @var        widths
   @vdesc      array of boundary widths for each variable
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @var        table_handles
   @vdesc      array of table handles which hold extra arguments
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @returntype CCTK_INT
   @returndesc
               return code of @seeroutine ApplyBndFlat
               -21 error reading boundary width array from table
               -22 wrong size boundary width array in table
   @endreturndesc
@@*/
CCTK_INT Bndry_Flat(const cGH *GH, const CCTK_INT num_vars, const CCTK_INT *vars,
                 const CCTK_INT *faces, const CCTK_INT *widths, const CCTK_INT *tables) {
  int i, j, k, gi, gdim, max_gdim, err, retval;

  /* variables to pass to ApplyBndFlat */
  CCTK_INT *width_alldirs; /* width of boundary in all directions */
  int dir;                 /* direction in which to apply bc */

  retval = 0;
  width_alldirs = NULL;
  max_gdim = 0;

  /* loop through variables, j at a time */
  for (i = 0; i < num_vars; i += j) {
    /* find other adjacent vars which are selected for identical bcs */
    j = 1;
    /* Since GFs are allowed to have different staggering, the best we
       can do is find variables of the same group which are selected
       for identical bcs.  If all GFs had the same staggering then we
       could group many GFs together. */
    gi = CCTK_GroupIndexFromVarI(vars[i]);
    while (i + j < num_vars && vars[i + j] == vars[i] + j &&
           CCTK_GroupIndexFromVarI(vars[i + j]) == gi &&
           tables[i + j] == tables[i] && faces[i + j] == faces[i] &&
           widths[i + j] == widths[i]) {
      ++j;
    }

    dir = 0; /* apply bc to all faces */

    /* Determine boundary width on all faces */
    /* allocate memory for buffer */
    gdim = CCTK_GroupDimI(gi);
    if (gdim > max_gdim) {
      width_alldirs =
          (CCTK_INT *)realloc(width_alldirs, 2 * gdim * sizeof(CCTK_INT));
      max_gdim = gdim;
    }

    /* fill it with values, either from table or the boundary_width
       parameter */
    if (widths[i] < 0) {
      err = Util_TableGetIntArray(tables[i], 2 * gdim, width_alldirs,
                                  "BOUNDARY_WIDTH");
      if (err < 0) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Error %d when reading boundary width array from table "
                   "for %s",
                   err, CCTK_VarName(vars[i]));
        return -21;
      } else if (err != 2 * gdim) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Boundary width array for %s has %d elements, but %d "
                   "expected",
                   CCTK_VarName(vars[i]), err, 2 * gdim);
        return -22;
      }
    } else {
      for (k = 0; k < 2 * gdim; ++k) {
        width_alldirs[k] = widths[i];
      }
    }

    /* Apply the boundary condition */
    if ((retval = ApplyBndFlat(GH, 0, width_alldirs, dir, faces[i], vars[i],
                               j)) < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "ApplyBndFlat() returned %d", retval);
    }
  }
#ifdef DEBUG
  printf("BndFlat(): returning %d\n", retval);
#endif

  free(width_alldirs);

  return retval;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

/* maximum dimension we can deal with */
#define MAXDIM 3

/* macro to compute the linear index of a 3D point */
#define INDEX_3D(ash, i, j, k) ((i) + (ash)[0] * ((j) + (ash)[1] * (k)))

/*@@
   @routine    FLAT_BOUNDARY
   @date       Tue 10 Apr 2001
   @author     Thomas Radke
   @desc
               Macro to apply flat boundary conditions to a variable
               Currently it is limited up to 3D variables only.
   @enddesc

   @var        doBC
   @vdesc      flag telling whether to apply boundary conditions or not
   @vtype      int
   @vio        in
   @endvar
   @var        iend, jend, kend
   @vdesc      upper ranges for the loopers
   @vtype      int
   @vio        in
   @endvar
   @var        ii_to, jj_to, kk_to
   @vdesc      indices of the current grid point to copy to
   @vtype      int
   @vio        in
   @endvar
   @var        ii_from, jj_from, kk_from
   @vdesc      indices of the current grid point to copy from
   @vtype      int
   @vio        in
   @endvar
@@*/
#define FLAT_BOUNDARY(doBC, iend, jend, kend, ii_to, jj_to, kk_to, ii_from,    \
                      jj_from, kk_from)                                        \
  {                                                                            \
    if (doBC) {                                                                \
      for (k = 0; k < kend; k++) {                                             \
        for (j = 0; j < jend; j++) {                                           \
          for (i = 0; i < iend; i++) {                                         \
            int _index_to, _index_from;                                        \
                                                                               \
            _index_to = INDEX_3D(ash, ii_to, jj_to, kk_to) * vtypesize;        \
            _index_from =                                                      \
                INDEX_3D(ash, ii_from, jj_from, kk_from) * vtypesize;          \
            memcpy((char *)GH->data[var][timelvl] + _index_to,                 \
                   (char *)GH->data[var][timelvl] + _index_from, vtypesize);   \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }

/*@@
   @routine    ApplyBndFlat
   @date       Jul 5 2000
   @author     Gabrielle Allen, Gerd Lanfermann
   @desc
               Apply flat boundary conditions to a group of grid functions
               given by their indices
               This routine is called by the various BndFlatXXX wrappers.

               Although it is currently limited to handle 1D, 2D, or 3D
               variables only it can easily be extended for higher dimensions
               by adapting the appropriate macros.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        width_dir
   @vdesc      boundary width in direction dir
   @vtype      CCTK_INT
   @vio        in
   @endvar
   @var        in_widths
   @vdesc      boundary widths for all directions
   @vtype      CCTK_INT [ dimension of variable(s) ]
   @vio        in
   @endvar
   @var        dir
   @vdesc      direction to set boundaries (0 for setting all directions)
   @vtype      int
   @vio        in
   @endvar
   @var        first_var
   @vdesc      index of first variable to apply boundaries to
   @vtype      int
   @vio        in
   @endvar
   @var        num_vars
   @vdesc      number of variables
   @vtype      int
   @vio        in
   @endvar

   @calls      CCTK_GroupIndexFromVarI
               CCTK_GroupDimI
               CCTK_VarTypeI
               CCTK_GroupStaggerDirArrayGI
               FLAT_BOUNDARY
   @history
   @hdate      Tue 10 Apr 2001
   @hauthor    Thomas Radke
   @hdesc      Merged separate routines for 1D, 2D, and 3D
               into a single generic routine
   @endhistory

   @returntype int
   @returndesc
                0 for success
               -1 if dimension is not supported
               -2 if direction parameter is invalid
               -3 if boundary width array parameter is NULL
   @endreturndesc
@@*/
static int ApplyBndFlat(const cGH *GH, CCTK_INT width_dir,
                        const CCTK_INT *in_widths,
                        int dir, CCTK_INT faces,
                        int first_var, int num_vars) {
  int i, j, k;
  int var, vtypesize, gindex, gdim, timelvl;
  int doBC[2 * MAXDIM], ash[MAXDIM], lsh[MAXDIM];
  CCTK_INT widths[2 * MAXDIM];
  CCTK_INT symtable;
  CCTK_INT symbnd[2 * MAXDIM];
  CCTK_INT is_physical[2 * MAXDIM];
  CCTK_INT ierr;

  /* get the group index of the variables */
  gindex = CCTK_GroupIndexFromVarI(first_var);

  /* get the number of dimensions and the size of the variables' type */
  gdim = CCTK_GroupDimI(gindex);
  vtypesize = CCTK_VarTypeSize(CCTK_VarTypeI(first_var));

  /* make sure we can deal with this number of dimensions */
  if (gdim > MAXDIM) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "ApplyBndFlat: Variable dimension of %d not supported", gdim);
    return (-1);
  }

  /* check the direction parameter */
  if (abs(dir) > gdim) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "ApplyBndFlat: direction %d greater than dimension %d", dir,
               gdim);
    return (-2);
  }

  /* set up boundary width array */
  if (dir) {
    widths[2 * (abs(dir) - 1)] = width_dir;
    widths[2 * (abs(dir) - 1) + 1] = width_dir;
  } else if (in_widths) {
    memcpy(widths, in_widths, 2 * gdim * sizeof *widths);
  } else {
    CCTK_WARN(1, "ApplyBndFlat: NULL pointer passed for boundary width "
                 "array");
    return (-3);
  }

  /* initialize arrays for variables with less dimensions than MAXDIM
     so that we can use the INDEX_3D macro later on */
  for (i = gdim; i < MAXDIM; i++) {
    ash[i] = 1;
    lsh[i] = 1;
  }

  /* get the current timelevel */
  timelvl = 0;

  /* see if we have a physical boundary */
  symtable = SymmetryTableHandleForGrid(GH);
  if (symtable < 0)
    CCTK_WARN(0, "internal error");
  ierr = Util_TableGetIntArray(symtable, 2 * gdim, symbnd, "symmetry_handle");
  if (ierr != 2 * gdim)
    CCTK_WARN(0, "internal error");
  for (i = 0; i < 2 * gdim; i++) {
    is_physical[i] = symbnd[i] < 0;
  }

  /* sanity check on width of boundary,  */
  BndSanityCheckWidths2(GH, first_var, gdim, widths, "Flat");

  /* now loop over all variables */
  for (var = first_var; var < first_var + num_vars; var++) {
    /* Apply condition if:
       + boundary is a physical boundary
       + boundary is an outer boundary
       + have enough grid points
    */
    for (i = 0; i < 2 * gdim; i++) {
      doBC[i] = is_physical[i] && (faces == CCTK_ALL_FACES || (faces & (1<<i)));
    }
    for (i = 0; i < gdim; i++) {
      ash[i] = GH->cctk_ash[i];
      lsh[i] = GH->cctk_lsh[i];
      doBC[i * 2] &= GH->cctk_lsh[i] > widths[i * 2] && GH->cctk_bbox[i * 2];
      doBC[i * 2 + 1] &=
          GH->cctk_lsh[i] > widths[i * 2 + 1] && GH->cctk_bbox[i * 2 + 1];
      if (dir != 0) {
        doBC[i * 2] &= (dir < 0 && (i + 1 == abs(dir)));
        doBC[i * 2 + 1] &= (dir > 0 && (i + 1 == abs(dir)));
      }
    }

    /* now apply the boundaries face by face */
    if (gdim > 0) {
#ifdef DEBUG_BOUNDARY
      if (doBC[0]) {
        printf("Boundary: Applying lower x flat boundary condition\n");
      }
      if (doBC[1]) {
        printf("Boundary: Applying upper x flat boundary condition\n");
      }
#endif /* DEBUG_BOUNDARY */
      /* lower x */
      FLAT_BOUNDARY(doBC[0], widths[0], lsh[1], lsh[2], i, j, k, widths[0], j,
                    k);
      /* upper x */
      FLAT_BOUNDARY(doBC[1], widths[1], lsh[1], lsh[2], lsh[0] - i - 1, j, k,
                    lsh[0] - widths[1] - 1, j, k);
    }
    if (gdim > 1) {
#ifdef DEBUG_BOUNDARY
      if (doBC[2]) {
        printf("Boundary: Applying lower y flat boundary condition\n");
      }
      if (doBC[3]) {
        printf("Boundary: Applying upper y flat boundary condition\n");
      }
#endif /* DEBUG_BOUNDARY */
      /* lower y */
      FLAT_BOUNDARY(doBC[2], lsh[0], widths[2], lsh[2], i, j, k, i, widths[2],
                    k);
      /* upper y */
      FLAT_BOUNDARY(doBC[3], lsh[0], widths[3], lsh[2], i, lsh[1] - j - 1, k, i,
                    lsh[1] - widths[3] - 1, k);
    }
    if (gdim > 2) {
#ifdef DEBUG_BOUNDARY
      if (doBC[4]) {
        printf("Boundary: Applying lower z flat boundary condition\n");
      }
      if (doBC[5]) {
        printf("Boundary: Applying upper z flat boundary condition\n");
      }
#endif /* DEBUG_BOUNDARY */
      /* lower z */
      FLAT_BOUNDARY(doBC[4], lsh[0], lsh[1], widths[4], i, j, k, i, j,
                    widths[4]);
      /* upper z */
      FLAT_BOUNDARY(doBC[5], lsh[0], lsh[1], widths[5], i, j, lsh[2] - k - 1, i,
                    j, lsh[2] - widths[5] - 1);
    }
  }

  return (0);
}

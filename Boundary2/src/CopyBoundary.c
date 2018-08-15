/*@@
  @file      CopyBoundary.c
  @date      Mon Mar 15 15:09:00 1999
  @author    Gerd Lanfermann, Gabrielle Allen
  @desc
             Routines for applying copying-boundary conditions
  @enddesc
  @history
  @hdate     Sun 25 Feb 2001
  @hauthor   Thomas Radke
  @hdesc     BC routines generalized for applying to arbitrary CCTK data types
  @endhistory
  @version   $Id$
@@*/

#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "util_Table.h"
#include "util_ErrorCodes.h"
#include "cctk_FortranString.h"

static int ApplyBndCopy(const cGH *GH, CCTK_INT stencil_dir,
                        const CCTK_INT *stencil_alldirs,
                        int dir, CCTK_INT faces,
                        int first_var_to, int first_var_from, int num_vars);

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
/*@@
   @routine    BndCopy
   @date       13 Feb 2003
   @author     David Rideout
   @desc
               Top level function which is registered as handling
               the Copy boundary condition
   @enddesc
   @calls      ApplyBndCopy
               CCTK_GroupDimFromVarI
               Util_TableGetIntArray
               Util_TableQueryValueInfo
               CCTK_VWarn
               Util_TableGetString
               CCTK_VarIndex
               Util_TableGetInt

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
               return code of @seeroutine ApplyBndCopy
               -11 invalid table handle
               -12 no "COPY_FROM" key in table
               -21 error reading boundary width array from table
               -22 wrong size boundary width array in table
   @endreturndesc
@@*/
CCTK_INT Bndry_Copy(const cGH *GH, CCTK_INT num_vars, CCTK_INT *vars,
                 CCTK_INT *faces, CCTK_INT *widths, CCTK_INT *tables) {
  int i, j, k, gi, gdim, max_gdim, err, retval;
  CCTK_INT value_type, value_size;
  char *copy_from_name;

  /* variables to pass to ApplyBndCopy */
  CCTK_INT *width_alldirs; /* width of boundary on each face */
  int dir;                 /* direction in which to apply bc */
  CCTK_INT
      copy_from; /* variable (index) from which to copy the boundary data */

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
       could groups many GFs together. */
    gi = CCTK_GroupIndexFromVarI(vars[i]);
    while (i + j < num_vars && vars[i + j] == vars[i] + j &&
           CCTK_GroupIndexFromVarI(vars[i + j]) == gi &&
           tables[i + j] == tables[i] && faces[i + j] == faces[i] &&
           widths[i + j] == widths[i]) {
      ++j;
    }

    /* Check to see if faces specification is valid */
    if (faces[i] != CCTK_ALL_FACES) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Faces specification %d for Copy boundary conditions on "
                 "%s is not implemented yet.  "
                 "Applying Copy bcs to all (external) faces.",
                 (int)faces[i], CCTK_VarName(vars[i]));
    }
    dir = 0; /* apply bc to all faces */

    /* Look on table for copy-from variable */
    err = Util_TableQueryValueInfo(tables[i], &value_type, &value_size,
                                   "COPY_FROM");
    if (err == UTIL_ERROR_BAD_HANDLE) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Invalid table handle passed for Copy boundary "
                 "conditions for %s.  Name or index of variable to copy from "
                 "must be provided via key \"COPY_FROM\".  Aborting.",
                 CCTK_VarName(vars[i]));
      return -11;
    } else if (err == 1) {
      if (value_type == CCTK_VARIABLE_STRING) {
        copy_from_name = malloc(value_size * sizeof(char));
        Util_TableGetString(tables[i], value_size, copy_from_name, "COPY_FROM");
        copy_from = CCTK_VarIndex(copy_from_name);
        free(copy_from_name);
      } else if (value_type == CCTK_VARIABLE_INT) {
        Util_TableGetInt(tables[i], &copy_from, "COPY_FROM");
      } else {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Invalid data type for key \"COPY_FROM\" "
                   "Please use CCTK_STRING for the variable name, "
                   "or CCTK_INT for the variable index.");
      }
    } else {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "No key \"COPY_FROM\" provided in table.  Please enter the "
                 "name or index of variable to copy from into the table "
                 "under this key.  Aborting.");
      return -12;
    }

    /* Determine boundary width on all faces */
    /* (re-)allocate memory for buffer */
    gdim = CCTK_GroupDimI(gi);
    if (gdim > max_gdim) {
      width_alldirs =
          (CCTK_INT *)realloc(width_alldirs, 2 * gdim * sizeof(CCTK_INT));
      max_gdim = gdim;
    }

    /* fill it with values, either from table or the boundary_width
       parameter */
    if (widths[i] < 0) {
      err = Util_TableGetIntArray(tables[i], gdim, width_alldirs,
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
    if (!retval &&
        (retval = ApplyBndCopy(GH, 0, width_alldirs, dir, vars[i], faces[i],
                               copy_from, j)) < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "ApplyBndCopy() returned %d", retval);
    }
  }
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
   @routine    COPY_BOUNDARY
   @date       Sat 20 Jan 2001
   @author     Thomas Radke
   @desc
               Macro to apply copy boundary conditions to a variable
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
   @var        ii, jj, kk
   @vdesc      indices of the current grid point
   @vtype      int
   @vio        in
   @endvar
@@*/
#define COPY_BOUNDARY(doBC, iend, jend, kend, ii, jj, kk)                      \
  {                                                                            \
    if (doBC) {                                                                \
      for (k = 0; k < kend; k++) {                                             \
        for (j = 0; j < jend; j++) {                                           \
          for (i = 0; i < iend; i++) {                                         \
            int _index;                                                        \
                                                                               \
            _index = INDEX_3D(ash, ii, jj, kk) * vtypesize;                    \
            memcpy((char *)GH->data[var_to][timelvl_to] + _index,              \
                   (char *)GH->data[var_from][timelvl_from] + _index,          \
                   vtypesize);                                                 \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }

/*@@
   @routine    ApplyBndCopy
   @date       Thu Mar  2 11:02:10 2000
   @author     Gerd Lanfermann
   @desc
               Apply copy boundary conditions to a group of grid functions
               given by their indices
               This routine is called by the various BndCopyXXX wrappers.

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
   @vdesc      direction to copy boundaries (0 for copying all directions)
   @vtype      int
   @vio        in
   @endvar
   @var        first_var_to
   @vdesc      index of first variable to copy boundaries to
   @vtype      int
   @vio        in
   @endvar
   @var        first_var_from
   @vdesc      index of first variable to copy boundaries from
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
               COPY_BOUNDARY
   @history
   @hdate      Sat 20 Jan 2001
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
static int ApplyBndCopy(const cGH *GH, CCTK_INT width_dir,
                        const CCTK_INT *in_widths, int dir, CCTK_INT faces,
                        int first_var_to, int first_var_from, int num_vars) {
  int i, j, k;
  int timelvl_to, timelvl_from;
  int gindex, gdim;
  int var_to, var_from, vtypesize;
  int doBC[2 * MAXDIM], ash[MAXDIM], lsh[MAXDIM];
  CCTK_INT widths[2 * MAXDIM];
  CCTK_INT symtable;
  CCTK_INT symbnd[2 * MAXDIM];
  CCTK_INT is_physical[2 * MAXDIM];
  CCTK_INT ierr;

  /* get the group index of the target variable */
  gindex = CCTK_GroupIndexFromVarI(first_var_to);

  /* get the number of dimensions and the size of the variable's type */
  gdim = CCTK_GroupDimI(gindex);
  vtypesize = CCTK_VarTypeSize(CCTK_VarTypeI(first_var_to));

  /* make sure we can deal with this number of dimensions */
  if (gdim > MAXDIM) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Variable dimension of %d not supported", gdim);
    return (-1);
  }

  /* check the direction parameter */
  if (abs(dir) > gdim) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "ApplyBndCopy: direction %d greater than dimension %d", dir,
               gdim);
    return (-2);
  }

  /* set up stencil width array */
  if (dir) {
    widths[2 * (abs(dir) - 1)] = width_dir;
    widths[2 * (abs(dir) - 1) + 1] = width_dir;
  } else if (in_widths) {
    memcpy(widths, in_widths, 2 * gdim * sizeof *widths);
  } else {
    CCTK_WARN(1, "ApplyBndCopy: NULL pointer passed for boundary width "
                 "array");
    return (-3);
  }

  /* sanity check on width of boundary,  */
  BndSanityCheckWidths2(GH, first_var_to, gdim, widths, "Copy");

  /* initialize arrays for variables with less dimensions than MAXDIM
     so that we can use the INDEX_3D macro later on */
  for (i = gdim; i < MAXDIM; i++) {
    ash[i] = 1;
    lsh[i] = 1;
  }

  /* get the current timelevel */
  timelvl_to = 0;
  timelvl_from = 0;

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

  /* now loop over all variables */
  for (var_to = first_var_to, var_from = first_var_from;
       var_to < first_var_to + num_vars; var_to++, var_from++) {
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

    /* now copy the boundaries face by face */
    if (gdim > 0) {
      /* lower x */
      COPY_BOUNDARY(doBC[0], widths[0], lsh[1], lsh[2], i, j, k);
      /* upper x */
      COPY_BOUNDARY(doBC[1], widths[1], lsh[1], lsh[2], lsh[0] - i - 1, j, k);
    }
    if (gdim > 1) {
      /* lower y */
      COPY_BOUNDARY(doBC[2], lsh[0], widths[2], lsh[2], i, j, k);
      /* upper y */
      COPY_BOUNDARY(doBC[3], lsh[0], widths[3], lsh[2], i, lsh[1] - j - 1, k);
    }
    if (gdim > 2) {
      /* lower z */
      COPY_BOUNDARY(doBC[4], lsh[0], lsh[1], widths[4], i, j, k);
      /* upper z */
      COPY_BOUNDARY(doBC[5], lsh[0], lsh[1], widths[5], i, j, lsh[2] - k - 1);
    }
  }

  return (0);
}

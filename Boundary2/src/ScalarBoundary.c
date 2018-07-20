/*@@
  @file      ScalarBoundary.c
  @date      Mon Mar 15 15:09:00 1999
  @author    Gabrielle Allen, Gerd Lanfermann
  @desc
             Routines for applying scalar boundary conditions
  @enddesc
  @history
  @hdate     Tue 10 Apr 2001
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
#include "cctk_Parameters.h"
#include "cctk_FortranString.h"

#include "Boundary2.h"

static int ApplyBndScalar(const cGH *GH,
                          CCTK_INT width_dir, const CCTK_INT *in_widths,
                          int dir, CCTK_INT faces,
                          CCTK_REAL scalar, int first_var, int num_vars);

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/

/*@@
   @routine    BndScalar
   @date       13 Feb 2003
   @author     David Rideout
   @desc
               Top level function which is registered as handling
               the Scalar boundary condition
   @enddesc
   @calls      ApplyBndScalar

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        num_vars
   @vdesc      number of variables passed in through var_indices[]
   @vtype      int
   @vio        in
   @endvar
   @var        var_indices
   @vdesc      array of variable indicies to which to apply this boundary
               condition
   @vtype      int *
   @vio        in
   @endvar
   @var        faces
   @vdesc      array of set of faces to which to apply the bc
   @vtype      int
   @vio        in
   @endvar
   @var        widths
   @vdesc      array of boundary widths for each variable
   @vtype      int
   @vio        in
   @endvar
   @var        table_handles
   @vdesc      array of table handles which hold extra arguments
   @vtype      int
   @vio        in
   @endvar
   @returntype int
   @returndesc
               return code of @seeroutine ApplyBndScalar
               -21 error reading boundary width array from table
               -22 wrong size boundary width array in table
   @endreturndesc
@@*/
CCTK_INT Bndry_Scalar(const cGH *GH, CCTK_INT num_vars, CCTK_INT *vars,
                   CCTK_INT *faces, CCTK_INT *widths, CCTK_INT *tables) {
  int i, j, k, gi, gdim, max_gdim, err, retval;

  /* variables to pass to ApplyBndScalar */
  CCTK_INT *width_alldirs; /* width of stencil in all directions */
  int dir;                 /* direction in which to apply bc */
  CCTK_REAL scalar;

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
                 "Faces specification %d for Scalar boundary conditions on "
                 "%s is not implemented yet.  "
                 "Applying Scalar bcs to all (external) faces.",
                 (int)faces[i], CCTK_VarName(vars[i]));
    }
    dir = 0; /* apply bc to all faces */

    /* Set up default arguments for ApplyBndScalar */
    scalar = 0.;

    /* Look on table for possible non-default arguments
     * (If any of these table look-ups fail, the value will be unchanged
     * from its default value)
     */
    /* Scalar value */
    err = Util_TableGetReal(tables[i], &scalar, "SCALAR");
    if (err == UTIL_ERROR_BAD_HANDLE) {
      CCTK_VWarn(5, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Invalid table handle passed for Scalar boundary "
                 "conditions for %s.  Using all default values.",
                 CCTK_VarName(vars[i]));
    }

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
    if ((retval = ApplyBndScalar(GH, 0, width_alldirs, dir, faces[i], scalar,
                                 vars[i], j)) < 0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "ApplyBndScalar() returned %d", retval);
    }
  }
#ifdef DEBUG
  printf("BndScalar(): returning %d\n", err);
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

/* an empty macro */
#define NOTHING

/*@@
   @routine    SCALAR_BOUNDARY_TYPED
   @date       Sat 20 Jan 2001
   @author     Thomas Radke
   @desc
               Macro to apply scalar boundary conditions to a variable
               of given datatype
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
   @var        cctk_type
   @vdesc      CCTK datatype of the variable
   @vtype      <cctk_type>
   @vio        in
   @endvar
@@*/
#define SCALAR_BOUNDARY_TYPED(doBC, iend, jend, kend, ii, jj, kk,              \
                              left_cctk_type, right_cctk_type)                 \
  {                                                                            \
    if (doBC) {                                                                \
      for (k = 0; k < kend; k++) {                                             \
        for (j = 0; j < jend; j++) {                                           \
          for (i = 0; i < iend; i++) {                                         \
            int _index;                                                        \
                                                                               \
            _index = INDEX_3D(ash, ii, jj, kk);                                \
            ((left_cctk_type *)GH->data[var][timelvl])[_index] =               \
                (right_cctk_type)scalar;                                       \
          }                                                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }

/*@@
   @routine    SCALAR_BOUNDARY
   @date       Sat 20 Jan 2001
   @author     Thomas Radke
   @desc
               Macro to apply scalar boundary conditions to a variable
               of a given datatype in all directions
               Currently it is limited up to 3D variables only.
   @enddesc

   @var        cctk_type
   @vdesc      CCTK datatype of the variable
   @vtype      <cctk_type>
   @vio        in
   @endvar
@@*/
#define SCALAR_BOUNDARY(left_cctk_type, right_cctk_type)                       \
  {                                                                            \
    /* now set the boundaries face by face */                                  \
    if (gdim > 0) {                                                            \
      /* lower x */                                                            \
      SCALAR_BOUNDARY_TYPED(doBC[0], widths[0], lsh[1], lsh[2], i, j, k,       \
                            left_cctk_type, right_cctk_type);                  \
      /* upper x */                                                            \
      SCALAR_BOUNDARY_TYPED(doBC[1], widths[1], lsh[1], lsh[2],                \
                            lsh[0] - i - 1, j, k, left_cctk_type,              \
                            right_cctk_type);                                  \
    }                                                                          \
    if (gdim > 1) {                                                            \
      /* lower y */                                                            \
      SCALAR_BOUNDARY_TYPED(doBC[2], lsh[0], widths[2], lsh[2], i, j, k,       \
                            left_cctk_type, right_cctk_type);                  \
      /* upper y */                                                            \
      SCALAR_BOUNDARY_TYPED(doBC[3], lsh[0], widths[3], lsh[2], i,             \
                            lsh[1] - j - 1, k, left_cctk_type,                 \
                            right_cctk_type);                                  \
    }                                                                          \
    if (gdim > 2) {                                                            \
      /* lower z */                                                            \
      SCALAR_BOUNDARY_TYPED(doBC[4], lsh[0], lsh[1], widths[4], i, j, k,       \
                            left_cctk_type, right_cctk_type);                  \
      /* upper z */                                                            \
      SCALAR_BOUNDARY_TYPED(doBC[5], lsh[0], lsh[1], widths[5], i, j,          \
                            lsh[2] - k - 1, left_cctk_type, right_cctk_type);  \
    }                                                                          \
  }

/*@@
  @routine    ApplyBndScalar
  @date       Tue Jul 18 18:10:33 2000
  @author     Gerd Lanfermann
  @desc
              Apply scalar boundary conditions to a group of grid functions
              given by their indices
              This routine is called by the various BndScalarXXX wrappers.

              Although it is currently limited to handle 3D variables only
              it can easily be extended for other dimensions
              by adapting the appropriate macros.
  @enddesc
  @calls      CCTK_VarTypeI
              CCTK_GroupDimFromVarI
              SCALAR_BOUNDARY

  @var        GH
  @vdesc      Pointer to CCTK grid hierarchy
  @vtype      const cGH *
  @vio        in
  @endvar
  @var        width_dir
  @vdesc      boundary width in direction dir
  @vtype      int
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
  @var        scalar
  @vdesc      scalar value to set the boundaries to
  @vtype      CCTK_REAL
  @vio        in
  @endvar
  @var        first_var
  @vdesc      index of first variable to apply boundary conditions to
  @vtype      int
  @vio        in
  @endvar
  @var        num_vars
  @vdesc      number of variables
  @vtype      int
  @vio        in
  @endvar

  @history
  @hdate      Tue 10 Apr 2001
  @hauthor    Thomas Radke
  @hdesc      Merged separate routines for 1D, 2D, and 3D
              into a single generic routine
  @endhistory

  @returntype int
  @returndesc
               0 for success
              -1 if abs(direction) is greater than variables' dimension
              -2 if variable dimension is not supported
              -3 if NULL pointer passed as boundary width array
              -4 if variable type is not supported
  @endreturndesc
@@*/
static int ApplyBndScalar(const cGH *GH,
                          CCTK_INT width_dir, const CCTK_INT *in_widths,
                          int dir, CCTK_INT faces,
                          CCTK_REAL scalar, int first_var, int num_vars) {
  int ierr;
  int i, j, k;
  int gindex, gdim;
  int var, timelvl;
  int doBC[2 * MAXDIM], ash[MAXDIM], lsh[MAXDIM];
  CCTK_INT widths[2 * MAXDIM];
  CCTK_INT symtable;
  CCTK_INT symbnd[2 * MAXDIM];
  CCTK_INT is_physical[2 * MAXDIM];

  /* check the direction parameter */
  if (abs(dir) > MAXDIM) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "ApplyBndScalar: direction %d is greater than maximum "
               "dimension %d",
               dir, MAXDIM);
    return (-1);
  }

  /* get the group index and dimension */
  gindex = CCTK_GroupIndexFromVarI(first_var);
  gdim = CCTK_GroupDimI(gindex);

  /* check the dimension */
  if (gdim > MAXDIM) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "ApplyBndScalar: variable dimension of %d not supported", gdim);
    return (-2);
  }

  /* set up boundary width array */
  if (dir) {
    widths[2 * (abs(dir) - 1)] = width_dir;
    widths[2 * (abs(dir) - 1) + 1] = width_dir;
  } else if (in_widths) {
    memcpy(widths, in_widths, 2 * gdim * sizeof *widths);
  } else {
    CCTK_WARN(1, "ApplyBndScalar: NULL pointer passed "
                 "for boundary width array");
    return (-3);
  }

  /* sanity check on width of boundary,  */
  BndSanityCheckWidths2(GH, first_var, gdim, widths, "Scalar");

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

    switch (CCTK_VarTypeI(var)) {
    /* FIXME: can't pass an empty preprocessor constant as a macro argument
              on some systems (e.g. MacOS X), so we have to define it outside */
    case CCTK_VARIABLE_BYTE:
      SCALAR_BOUNDARY(CCTK_BYTE, CCTK_BYTE);
      break;

    case CCTK_VARIABLE_INT:
      SCALAR_BOUNDARY(CCTK_INT, CCTK_INT);
      break;

    case CCTK_VARIABLE_REAL:
      SCALAR_BOUNDARY(CCTK_REAL, CCTK_REAL);
      break;

#ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      SCALAR_BOUNDARY(CCTK_INT1, CCTK_INT1);
      break;
#endif

#ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      SCALAR_BOUNDARY(CCTK_INT2, CCTK_INT2);
      break;
#endif

#ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      SCALAR_BOUNDARY(CCTK_INT4, CCTK_INT4);
      break;
#endif

#ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      SCALAR_BOUNDARY(CCTK_INT8, CCTK_INT8);
      break;
#endif

#ifdef HAVE_CCTK_INT16
    case CCTK_VARIABLE_INT16:
      SCALAR_BOUNDARY(CCTK_INT16, CCTK_INT16);
      break;
#endif

#ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      SCALAR_BOUNDARY(CCTK_REAL4, CCTK_REAL4);
      break;
#endif

#ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      SCALAR_BOUNDARY(CCTK_REAL8, CCTK_REAL8);
      break;
#endif

#ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      SCALAR_BOUNDARY(CCTK_REAL16, CCTK_REAL16);
      break;
#endif

    case CCTK_VARIABLE_COMPLEX:
      SCALAR_BOUNDARY(CCTK_COMPLEX, CCTK_REAL);
      break;

    default:
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Unsupported variable type %d for variable '%s'",
                 CCTK_VarTypeI(var), CCTK_VarName(var));
      return (-4);
    }
  }

  return (0);
}

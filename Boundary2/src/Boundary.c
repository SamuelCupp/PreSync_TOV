/*@@
  @file      Boundary.c
  @date      Sat Oct 26 22:39:40 CEST 2002
  @author    David Rideout
  @desc
             Implements the new boundary specification.
  @enddesc
  @version   $Header$
@@*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include "util_String.h"
#include "Boundary.h"

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/
/* #define DEBUG 1 */

/* Linked list, called a var list, for holding variables selected for a bc:
 * Entries are sorted in the order they appear in struct BCVAR,
 * i.e. the var index varies more rapidly that table handle.
 * (Currently no sorting is done on faces specification.)
 *
 * There will be one such linked list for each type of boundary
 * condition selected (i.e. one for each bc_name).
 */
struct BCVAR {
  struct BCVAR *next; /* pointer to next entry in list */
  int faces;          /* set of faces for this application of bc */
  int width;          /* width of the boundary, if it is equal for all faces */
  int table;          /* table handle holding extra arguments */
  int var;            /* index of grid variable to which to apply the bc */
};

/*
 * Linked list of var lists, one for each type of requested bc
 * (i.e. one for each bc_name).
 *
 * Here is also recorded how many of each bc type have
 * been selected so far, so that the GetSelectedBCs doesn't have to be
 * run through twice; once simply to get the number of selected vars.
 *
 * This list is sorted by bc_name.  Alternatively one could sort it by
 * associated function pointer, but this seems a bit obtuse.
 */
struct BCDATA {
  struct BCDATA *next;    /* pointer to next element of this list */
  struct BCVAR *var_list; /* pointer to first element of a var list */
  const char *bc_name;    /* name of bc */
  int num;                /* number of entries for this bc in var list */
};

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/* Table for holding function pointers associated with each boundary condition:
 * This table has
 *        key   = boundary condition name (eg "Radiation")
 *        value = a CCTK_FPOINTER pointing to a function to implement that BC
 */
static int physbc_table_handle = -1;

/* Linked list for storing data associated with selections list itself */
static struct BCDATA *bcdata_list = NULL;

/* 'The' GH, i.e. to check that there is not more than one... */
static CCTK_POINTER_TO_CONST theGH = NULL;

/*@@
  @routine    Bdry_Boundary_SelectedGVs
  @date       Sun Nov  3 19:51:37 CET 2002
  @author     David Rideout
  @desc
	      Returns list of variable indices and table handles of
	      variables selected for boundary conditions.
  @enddesc
  @calls
  @history
  @endhistory
  @var        _GH
  @vdesc      cctkGH *

  @vtype      CCTK_POINTER_TO_CONST
  @vio        in
  @endvar
  @var        array_size
  @vdesc      size of arrays to which var_indices and table_handles point
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        var_indices
  @vdesc      array into which selected variable indices will be placed
  @vtype      CCTK_INT
  @vio        out
  @endvar
  @var        faces
  @vdesc      array into which a set of selected faces for variables selected
              for bc will be placed
  @vtype      CCTK_INT
  @vio        out
  @endvar
  @var        widths
  @vdesc      array into which boundary widths of selected variables will be
              placed
  @vtype      CCTK_INT
  @vio        in
  @endvar
  @var        table_handles
  @vdesc      array into which table_handles for variables selected for bc
              will be placed
  @vtype      CCTK_INT
  @vio        out
  @endvar
  @var        bc_name
  @vdesc      name of bc for which to get the selected vars,
              NULL returns all selected vars for all bcs
  @vtype      CCTK_STRING
  @vio        in
  @endvar
  @returntype CCTK_INT
  @returndesc
  -1 no boundary condition registered under bc_name
  -5 new value passed for GH
   number of variables selected for bc_name
  @endreturndesc
@@*/
CCTK_INT SelectedGVs(CCTK_POINTER_TO_CONST _GH, CCTK_INT array_size, 
                                   CCTK_INT *var_indices, CCTK_INT *faces,
                                   CCTK_INT *widths, CCTK_INT *table_handles,
                                   CCTK_STRING bc_name) {
  int retval, i, j;
  struct BCVAR *current;
  struct BCDATA *current_bcdata;
  const cGH *GH = _GH;

  current = NULL;
  retval = 0;

  /* Check to see if this is a new GH */
  if (!theGH) /* This is the first valid GH passed to a Boundary routine */
  {
    theGH = GH;
  } else if (GH != theGH) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "New GH passed to Boundary2_SelectedGVs.  "
               "Thorn CactusBase/Boundary does not yet handle multiple GHs "
               "properly.");
    retval = -5;
  }

#ifdef DEBUG
  printf("Boundary2_SelectedGVs: called with bc_name=\"%s\" array_size=%d\n",
         bc_name, array_size);
  fflush(stdout);
#endif

  i = 0; /* i indexes the location in the returned arrays */

  /* Step through bcdata list */
  for (current_bcdata = bcdata_list; current_bcdata;
       current_bcdata = current_bcdata->next) {
#ifdef DEBUG
    printf("  looping through bcdata list, at bcdata entry for %s bc\n",
           current_bcdata->bc_name);
#endif

    if (!bc_name || CCTK_Equals(current_bcdata->bc_name, bc_name)) {

      /* Add these selected vars to return value */
      retval += current_bcdata->num;
      current = current_bcdata->var_list;

      /* Loop through var list */
      for (j = 0; /* j counts the bcs to check internal consistency */
           i < array_size && current; current = current->next, ++i, ++j) {

#ifdef DEBUG
        printf("    looping through selected vars, at current->var_index = "
               "%d\n",
               current->var);
        printf("      current->next is %p\n", current->next);
#endif

        if (faces) {
          faces[i] = current->faces;
        }
        if (widths) {
          widths[i] = current->width;
        }
        if (table_handles) {
          table_handles[i] = current->table;
        }
        if (var_indices) {
          var_indices[i] = current->var;
        }
      }
      if (j > current_bcdata->num)
        CCTK_WARN(0, "internal error");
      if (i != array_size && j != current_bcdata->num)
        CCTK_WARN(0, "internal error");
    }
  }

  /* Warn if there is no bc registered under this name */
  if (bc_name &&
      !Util_TableQueryValueInfo(physbc_table_handle, NULL, NULL, bc_name)) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "There is no boundary condition registered under the name %s",
               bc_name);
    retval = -1;
  }

  return retval;
}

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <iostream>
#include <algorithm>
#include "Boundary.h"

typedef void (*boundary_function)(
  const cGH *cctkGH,
  int num_vars,
  int *var_indices,
  int *faces,
  int *widths,
  int *table_handles);

// From PreSync
extern "C" void Carpet_RegisterPhysicalBC(
    const cGH *cctkGH,
    boundary_function func,
    const char *bc_name,
    int before);

extern "C"
void Carpet_SelectGroupForBC(
    const cGH *cctkGH,
    int faces,
    int width,
    int table_handle,
    const char *group_name,
    const char *bc_name);

extern "C"
void presync_registerboundary(CCTK_ARGUMENTS)
{
  _DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  int err;
//  CCTK_INT group, rhs;

  std::cout << "Register Boundary Conditions" << std::endl;

  if (register_scalar) {
    Carpet_RegisterPhysicalBC(cctkGH, BndScalar, "Scalar",1);
  }

  if (register_flat) {
    Carpet_RegisterPhysicalBC(cctkGH, BndFlat, "Flat",1);
  }

  if (register_radiation) {
    Carpet_RegisterPhysicalBC(cctkGH, BndRadiative, "Radiation",1);
  }

  if (register_copy) {
    Carpet_RegisterPhysicalBC(cctkGH, BndCopy, "Copy",1);
  }

  if (register_robin) {
    Carpet_RegisterPhysicalBC(cctkGH, BndRobin, "Robin",1);
  }

  if (register_static) {
    Carpet_RegisterPhysicalBC(cctkGH, BndStatic, "Static",1);
  }

  if (register_none) {
    Carpet_RegisterPhysicalBC(cctkGH, BndNone, "None",1);
  }

  int w = 1;

  Carpet_SelectGroupForBC(cctkGH,
    CCTK_ALL_FACES, w,
   -1 /* no table */, "PresyncWave::evo_vars",
   "zero");

  Carpet_SelectGroupForBC(cctkGH,
    CCTK_ALL_FACES, w,
   -1 /* no table */, "PresyncWave::rhs_vars",
   "zero");
}

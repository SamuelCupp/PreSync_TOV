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
int Carpet_SelectGroupForBC(
    const cGH *cctkGH,
    int faces,
    int width,
    int table_handle,
    const char *group_name,
    const char *bc_name);

extern "C"
void Boundary2_RegisterBCs(CCTK_ARGUMENTS)
{
  _DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  std::cout << "Register Boundary Conditions" << std::endl;

  if (register_scalar) {
    Carpet_RegisterPhysicalBC(cctkGH, BndScalar, "scalar",1);
  }

  if (register_flat) {
    Carpet_RegisterPhysicalBC(cctkGH, BndFlat, "flat",1);
  }

  if (register_radiation) {
    Carpet_RegisterPhysicalBC(cctkGH, BndRadiative, "radiation",1);
  }

  if (register_copy) {
    Carpet_RegisterPhysicalBC(cctkGH, BndCopy, "copy",1);
  }

  if (register_robin) {
    Carpet_RegisterPhysicalBC(cctkGH, BndRobin, "robin",1);
  }

  if (register_static) {
    Carpet_RegisterPhysicalBC(cctkGH, BndStatic, "static",1);
  }

  if (register_none) {
    Carpet_RegisterPhysicalBC(cctkGH, BndNone, "none",1);
  }
}

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <algorithm>
#include "Boundary2.h"

void Boundary2_RegisterBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS

  CCTK_INFO("Registering Boundary Conditions");

  if (register_scalar) {
    int err = 0;
    err = Boundary_RegisterPhysicalBC(cctkGH, (boundary_function)Bndry_Scalar, "scalar");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Scalar\" "
                 "boundary condition",
                 err);
    }
  }

  if (register_flat) {
    int err = 0;
    err = Boundary_RegisterPhysicalBC(cctkGH, (boundary_function)Bndry_Flat, "flat");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Flat\" "
                 "boundary condition",
                 err);
    }
  }

  if (register_radiation) {
    int err = 0;
    err = Boundary_RegisterPhysicalBC(cctkGH, (boundary_function)Bndry_Radiative, "radiation");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Radiation\" "
                 "boundary condition",
                 err);
    }
  }

  if (register_copy) {
    int err = 0;
    err = Boundary_RegisterPhysicalBC(cctkGH, (boundary_function)Bndry_Copy, "copy");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Copy\" "
                 "boundary condition",
                 err);
    }
  }

  if (register_robin) {
    int err = 0;
    err = Boundary_RegisterPhysicalBC(cctkGH, (boundary_function)Bndry_Robin, "robin");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Robin\" "
                 "boundary condition",
                 err);
    }
  }

  if (register_static) {
    int err = 0;
    err = Boundary_RegisterPhysicalBC(cctkGH, (boundary_function)Bndry_Static, "static");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Static\" "
                 "boundary condition",
                 err);
    }
  }

  if (register_none) {
    int err = 0;
    err = Boundary_RegisterPhysicalBC(cctkGH, (boundary_function)Bndry_None, "none");
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"None\" "
                 "boundary condition",
                 err);
    }
  }
}

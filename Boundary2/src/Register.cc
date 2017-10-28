#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Arguments_Boundary2.h"
#include "cctk_Parameters.h"
#include <algorithm>
#include "Boundary2.h"

void Boundary2_RegisterBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_Boundary2_RegisterBCs
  DECLARE_CCTK_PARAMETERS

  int err = 0;
  CCTK_INFO("Registering Boundary Conditions");

  if (register_scalar) {
    RegisterPhysicalBC(cctkGH, Bndry_Scalar, "scalar", 1);
/*    err = Carpet_RegisterPhysicalBC(cctkGH, BndScalar, "scalar", 1);
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Scalar\" "
                 "boundary condition",
                 err);
    }
*/  }

  if (register_flat) {
    RegisterPhysicalBC(cctkGH, Bndry_Flat, "flat", 1);
/*    err = Carpet_RegisterPhysicalBC(cctkGH, BndFlat, "flat", 1);
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Flat\" "
                 "boundary condition",
                 err);
    }
*/  }

  if (register_radiation) {
    RegisterPhysicalBC(cctkGH, Bndry_Radiative, "radiation", 1);
/*    err = Carpet_RegisterPhysicalBC(cctkGH, BndRadiative, "radiation", 1);
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Radiation\" "
                 "boundary condition",
                 err);
    }
*/  }

  if (register_copy) {
    RegisterPhysicalBC(cctkGH, Bndry_Copy, "copy", 1);
/*    err = Carpet_RegisterPhysicalBC(cctkGH, BndCopy, "copy", 1);
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Copy\" "
                 "boundary condition",
                 err);
    }
*/  }

  if (register_robin) {
    RegisterPhysicalBC(cctkGH, Bndry_Robin, "robin", 1);
/*    err = Carpet_RegisterPhysicalBC(cctkGH, BndRobin, "robin", 1);
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Robin\" "
                 "boundary condition",
                 err);
    }
*/  }

  if (register_static) {
    RegisterPhysicalBC(cctkGH, Bndry_Static, "static", 1);
/*    err = Carpet_RegisterPhysicalBC(cctkGH, BndStatic, "static", 1);
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"Static\" "
                 "boundary condition",
                 err);
    }
*/  }

  if (register_none) {
    RegisterPhysicalBC(cctkGH, Bndry_None, "none", 1);
/*    err = Carpet_RegisterPhysicalBC(cctkGH, BndNone, "none", 1);
    if (err) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error %d when registering routine to handle \"None\" "
                 "boundary condition",
                 err);
    }
*/  }
}

#include <set>
#include <cstring>
#include <map>
#include <vector>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <math.h>
#include <array>
#include <iostream>
#include "PreSync.h"

namespace Carpet {

inline void cctk_assert_(int line,const char *file,const char *thorn,const char *str) {
  std::ostringstream msg; 
  msg << "Assertion Failed: " << str;
  CCTK_Error(line,file,thorn,msg.str().c_str());
}

#define CCTK_ASSERT(X) if(!(X)) cctk_assert_(__LINE__,__FILE__,CCTK_THORNSTRING,#X);

struct Bound {
  std::string bc_name;
  int faces;
  int width;
  int table_handle;
};

struct Func {
  boundary_function func;
  int before;
};

struct SymFunc {
  boundary_function func;
  int handle;
  int faces;
  int width[6];
};

std::map<std::string,Func> boundary_functions;
std::map<std::string,SymFunc> symmetry_functions;
/**
 * The index into the array is the same as the "before"
 * argument defined when registering a BC.
 */
std::array<std::map<int,std::vector<Bound>>,2> boundary_conditions;

extern "C"
CCTK_INT Bdry2_Boundary_RegisterPhysicalBC(
    const cGH *cctkGH,
    boundary_function func,
    const char *bc_name) {
    int before = 1;
  if(before != 0) before = 1;
  if(NULL==func) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,  
               "Physical Boundary condition '%s' points to NULL.", bc_name);
  }
  Func& f = boundary_functions[bc_name];
  f.func = func;
  f.before = before;
  return 0;
}

extern "C"
void Bdry2_Boundary_RegisterSymmetryBC(
    const cGH *cctkGH,
    boundary_function func,
    int handle,
    int faces,
    int width,
    const char *bc_name) {
  if(NULL==func) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,  
               "Symmetry Boundary condition '%s' points to NULL.", bc_name);
  }
  SymFunc& f = symmetry_functions[bc_name];
  f.func = func;
  f.handle = handle;
  f.faces = faces;
//  &f.width = width;
//  std::cout << "Register Width of " << f.width[0] << " for " << bc_name << std::endl;
}

void Boundary_SelectVarForBCI(
    const cGH *cctkGH,
    int faces,
    int width,
    int table_handle,
    int var_index,
    const char *bc_name) {
  if(!boundary_functions.count(bc_name)) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,  
               "Requested BC '%s' not found.", bc_name);
  }
  Func& f = boundary_functions.at(bc_name);
  CCTK_ASSERT(var_index != 0);
  std::vector<Bound>& bv = boundary_conditions[f.before][var_index];
  Bound b;
  b.faces = faces;
  b.width = width;
  b.table_handle = table_handle;
  b.bc_name = bc_name;
  bv.push_back(b);
}

extern "C"
void Bdry2_Boundary_SelectVarForBC(
    const cGH *cctkGH,
    int faces,
    int width,
    int table_handle,
    const char *var_name,
    const char *bc_name) {
  int i = CCTK_VarIndex(var_name);
  Boundary_SelectVarForBCI(cctkGH,faces,width,table_handle,i,bc_name);
}

extern "C"
CCTK_INT Bdry2_Boundary_SelectGroupForBC(
    const cGH *cctkGH,
    int faces,
    int width,
    int table_handle,
    const char *group_name,
    const char *bc_name) {
  int group = CCTK_GroupIndex(group_name);
  int vstart = CCTK_FirstVarIndexI(group);
  int vnum   = CCTK_NumVarsInGroupI(group);
  for(int i=vstart;i<vstart+vnum;i++) {
    Boundary_SelectVarForBCI(cctkGH,faces,width,table_handle,i,bc_name);
  }
  return 0;
}

extern "C"
void Boundary_ClearBCForVarI(
    const cGH *cctkGH,
    int var_index) {
  CCTK_ASSERT(var_index != 0);
  boundary_conditions[0][var_index].resize(0);
  boundary_conditions[1][var_index].resize(0);
}

}

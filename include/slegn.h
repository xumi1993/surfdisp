#pragma once
#include "surfdisp.h"   // for surfdisp::Layer, surfdisp::EarthModel
#include <vector>

namespace slegn {

struct EigenResult {
    double cp;                   // phase velocity (km/s)
    double cg;                   // group velocity (km/s)
    std::vector<double> disp;   // displacement eigenfunction [nlayer]
    std::vector<double> stress; // stress eigenfunction [nlayer]
    std::vector<double> dc2db;  // dc/dVs [nlayer]
    std::vector<double> dc2dh;  // dc/dh  [nlayer]
    std::vector<double> dc2dr;  // dc/drho [nlayer]
};

struct GroupKernelResult : EigenResult {
    std::vector<double> du2db;  // dU/dVs [nlayer]
    std::vector<double> du2dh;  // dU/dh  [nlayer]
    std::vector<double> du2dr;  // dU/drho [nlayer]
};

// Love-wave eigenfunctions + phase-velocity kernels at one period.
// phase_vel_km_s: phase velocity at period_s (from surfdisp::dispersion)
EigenResult slegn96(
    const std::vector<surfdisp::Layer>& model,
    double period_s,
    double phase_vel_km_s,
    surfdisp::EarthModel earth = surfdisp::EarthModel::Flat
);

// Array-based overload: thk/vs/rho are flat float arrays of length nlayer.
EigenResult slegn96(
    int nlayer,
    const float* thk, const float* vs, const float* rho,
    double period_s,
    double phase_vel_km_s,
    surfdisp::EarthModel earth = surfdisp::EarthModel::Flat
);

// Love-wave phase + group velocity kernels.
// t/cp: the target period; t1/cp1, t2/cp2: bracket for finite-difference group kernel
GroupKernelResult slegnpu(
    const std::vector<surfdisp::Layer>& model,
    double t,  double cp,
    double t1, double cp1,
    double t2, double cp2,
    surfdisp::EarthModel earth = surfdisp::EarthModel::Flat
);

// Array-based overload
GroupKernelResult slegnpu(
    int nlayer,
    const float* thk, const float* vs, const float* rho,
    double t,  double cp,
    double t1, double cp1,
    double t2, double cp2,
    surfdisp::EarthModel earth = surfdisp::EarthModel::Flat
);

} // namespace slegn

// Fortran bind(C) symbols for comparison tests.
// Argument order matches the Fortran subroutine definitions exactly;
// iflsph comes LAST (not after nlayer) in both routines.
extern "C" {

/* slegn96: thk,vs,rhom,nlayer, t,cp,cg, disp,stress, dc2db,dc2dh,dc2dr, iflsph */
void slegn96_(const float* thk, const float* vs, const float* rhom,
              int nlayer,
              const double* t, double* cp, double* cg,
              double* disp, double* stress,
              double* dc2db, double* dc2dh, double* dc2dr,
              int iflsph);

/* slegnpu: thk,vs,rhom,nlayer, t,cp,cg, disp,stress,
            t1,cp1, t2,cp2, dc2db,dc2dh,dc2dr, du2db,du2dh,du2dr, iflsph */
void slegnpu_(const float* thk, const float* vs, const float* rhom,
              int nlayer,
              const double* t,  double* cp,  double* cg,
              double* disp, double* stress,
              double* t1, double* cp1,
              double* t2, double* cp2,
              double* dc2db, double* dc2dh, double* dc2dr,
              double* du2db, double* du2dh, double* du2dr,
              int iflsph);

}

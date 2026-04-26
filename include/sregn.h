#pragma once
#include "surfdisp.h"
#include <vector>

namespace sregn {

struct EigenResult {
    double cp, cg;
    std::vector<double> dispu;   // ur: horizontal displacement [nlayer]
    std::vector<double> dispw;   // uz: vertical displacement [nlayer]
    std::vector<double> stressu; // tr: horizontal stress [nlayer]
    std::vector<double> stressw; // tz: vertical stress [nlayer]
    std::vector<double> dc2da;   // dc/dVp [nlayer]
    std::vector<double> dc2db;   // dc/dVs [nlayer]
    std::vector<double> dc2dh;   // dc/dh  [nlayer]
    std::vector<double> dc2dr;   // dc/drho [nlayer]
};

struct HTIResult : EigenResult {
    std::vector<double> dc2dgc;  // dc/dgc [nlayer]
    std::vector<double> dc2dgs;  // dc/dgs [nlayer]
};

struct GroupKernelResult : EigenResult {
    std::vector<double> du2da, du2db, du2dh, du2dr;
};

// vector<Layer> overloads
EigenResult sregn96(const std::vector<surfdisp::Layer>& model,
    double period_s, double phase_vel_km_s,
    surfdisp::EarthModel earth = surfdisp::EarthModel::Flat);

HTIResult sregn96_hti(const std::vector<surfdisp::Layer>& model,
    double period_s, double phase_vel_km_s,
    surfdisp::EarthModel earth = surfdisp::EarthModel::Flat);

GroupKernelResult sregnpu(const std::vector<surfdisp::Layer>& model,
    double t, double cp, double t1, double cp1, double t2, double cp2,
    surfdisp::EarthModel earth = surfdisp::EarthModel::Flat);

// Array overloads
EigenResult sregn96(int nlayer,
    const float* thk, const float* vp, const float* vs, const float* rho,
    double period_s, double phase_vel_km_s,
    surfdisp::EarthModel earth = surfdisp::EarthModel::Flat);

HTIResult sregn96_hti(int nlayer,
    const float* thk, const float* vp, const float* vs, const float* rho,
    double period_s, double phase_vel_km_s,
    surfdisp::EarthModel earth = surfdisp::EarthModel::Flat);

GroupKernelResult sregnpu(int nlayer,
    const float* thk, const float* vp, const float* vs, const float* rho,
    double t, double cp, double t1, double cp1, double t2, double cp2,
    surfdisp::EarthModel earth = surfdisp::EarthModel::Flat);

} // namespace sregn

// Fortran bind(C) symbols for test comparison
extern "C" {
void sregn96_(const float* thk, const float* vp, const float* vs, const float* rhom,
              int nlayer, double* t, double* cp, double* cg,
              double* dispu, double* dispw, double* stressu, double* stressw,
              double* dc2da, double* dc2db, double* dc2dh, double* dc2dr,
              int iflsph);

void sregn96_hti_(const float* thk, const float* vp, const float* vs, const float* rhom,
                  int nlayer, double* t, double* cp, double* cg,
                  double* dispu, double* dispw, double* stressu, double* stressw,
                  double* dc2da, double* dc2db, double* dc2dh, double* dc2dr,
                  double* dc2dgc, double* dc2dgs, int iflsph);

void sregnpu_(const float* thk, const float* vp, const float* vs, const float* rhom,
              int nlayer, double* t, double* cp, double* cg,
              double* dispu, double* dispw, double* stressu, double* stressw,
              double* t1, double* cp1, double* t2, double* cp2,
              double* dc2da, double* dc2db, double* dc2dh, double* dc2dr,
              double* du2da, double* du2db, double* du2dh, double* du2dr,
              int iflsph);
}

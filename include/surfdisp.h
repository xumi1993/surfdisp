#pragma once
#include <vector>
#include <stdexcept>

namespace surfdisp {

enum class EarthModel  { Flat, Spherical };
enum class WaveType    { Love, Rayleigh };
enum class VelocityType{ Phase, Group };

struct Layer {
    float thickness;  /* km  -- set to 0 for the halfspace (last layer) */
    float vp;         /* km/s */
    float vs;         /* km/s */
    float density;    /* g/cm^3 */
};

/*
 * dispersion()
 *
 * Compute surface-wave phase or group velocities for a 1-D layered model.
 *
 *   model      : layer stack, last entry is the halfspace (thickness ignored)
 *   periods_s  : periods in seconds (any order, any count)
 *   wave_type  : Love or Rayleigh
 *   mode       : 0 = fundamental mode, 1 = first overtone, …
 *   vel_type   : Phase or Group velocity
 *   earth      : Flat (default) or Spherical (applies earth-flattening)
 *
 * Returns a vector of velocities (km/s) matching periods_s in order.
 * Entries are 0 if a mode could not be found at that period.
 *
 * Throws std::invalid_argument for bad input.
 * Throws std::runtime_error if the fundamental mode search fails entirely.
 */
std::vector<double> dispersion(
    const std::vector<Layer>& model,
    const std::vector<double>& periods_s,
    WaveType    wave_type,
    int         mode     = 0,
    VelocityType vel_type = VelocityType::Phase,
    EarthModel  earth    = EarthModel::Flat
);

} // namespace surfdisp

/* ------------------------------------------------------------------ */
/* Legacy C-compatible wrappers (kept for the Fortran comparison test) */
/* ------------------------------------------------------------------ */
#ifdef __cplusplus
extern "C" {
#endif

void surfdisp96_cpp(const float *thkm, const float *vpm, const float *vsm,
                    const float *rhom, int nlayer, int iflsph, int iwave,
                    int mode, int igr, int kmax,
                    const double *t, double *cg);

void surfdisp96_f(const float *thkm, const float *vpm, const float *vsm,
                  const float *rhom, int nlayer, int iflsph, int iwave,
                  int mode, int igr, int kmax,
                  const double *t, double *cg);

#ifdef __cplusplus
}
#endif

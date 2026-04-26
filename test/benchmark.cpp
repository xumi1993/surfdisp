/*
 * benchmark.cpp
 * Wall-clock timing: surfdisp::dispersion (C++) vs Fortran surfdisp96
 *
 * Each implementation is called NREP times for every test case.
 * The first call is a warm-up and excluded from timing.
 */

#include "surfdisp.h"
#include <chrono>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <numeric>

using Clock = std::chrono::steady_clock;
using Ms    = std::chrono::duration<double, std::milli>;

static const std::vector<surfdisp::Layer> MODEL = {
    {  3.0f, 4.0f, 2.3f, 2.3f },
    {  7.0f, 5.8f, 3.3f, 2.7f },
    { 12.0f, 6.3f, 3.6f, 2.8f },
    { 15.0f, 6.8f, 3.9f, 2.9f },
    { 10.0f, 7.2f, 4.1f, 3.1f },
    { 10.0f, 7.5f, 4.3f, 3.2f },
    { 20.0f, 7.8f, 4.4f, 3.3f },
    { 30.0f, 8.0f, 4.5f, 3.4f },
    { 60.0f, 8.1f, 4.6f, 3.4f },
    {  0.0f, 8.2f, 4.7f, 3.4f },
};

static std::vector<double> log_periods(int n, double tmin, double tmax)
{
    std::vector<double> t(n);
    for (int i = 0; i < n; i++)
        t[i] = tmin * std::pow(tmax / tmin, (double)i / (n - 1));
    return t;
}

/* Build flat arrays once so Fortran wrapper overhead is not counted */
struct FlatModel {
    std::vector<float> thk, vp, vs, rho;
    int nlayer;
    FlatModel(const std::vector<surfdisp::Layer>& m) : nlayer((int)m.size()) {
        for (auto& l : m) {
            thk.push_back(l.thickness); vp.push_back(l.vp);
            vs.push_back(l.vs);         rho.push_back(l.density);
        }
    }
};

struct Case {
    const char*          label;
    surfdisp::WaveType   wave;
    int                  mode;   /* 0-based */
    surfdisp::VelocityType vtype;
};

static const Case CASES[] = {
    { "Love  mode=0 phase",     surfdisp::WaveType::Love,     0, surfdisp::VelocityType::Phase },
    { "Love  mode=0 group",     surfdisp::WaveType::Love,     0, surfdisp::VelocityType::Group },
    { "Love  mode=1 phase",     surfdisp::WaveType::Love,     1, surfdisp::VelocityType::Phase },
    { "Rayleigh mode=0 phase",  surfdisp::WaveType::Rayleigh, 0, surfdisp::VelocityType::Phase },
    { "Rayleigh mode=0 group",  surfdisp::WaveType::Rayleigh, 0, surfdisp::VelocityType::Group },
    { "Rayleigh mode=1 phase",  surfdisp::WaveType::Rayleigh, 1, surfdisp::VelocityType::Phase },
};

int main(int argc, char* argv[])
{
    const int  NPER = 200;          /* periods per call  */
    const int  NREP = 2000;         /* repetitions       */

    auto periods   = log_periods(NPER, 1.0, 100.0);
    FlatModel flat(MODEL);

    const int iflsph = 0;
    const int nlayer = flat.nlayer;
    const int kmax   = NPER;

    printf("Benchmark: %d periods, %d repetitions per case\n\n", NPER, NREP);
    printf("%-28s  %10s  %10s  %8s\n",
           "Case", "C++ (ms)", "F90 (ms)", "ratio F/C");
    printf("%s\n", std::string(64, '-').c_str());

    for (auto& c : CASES) {
        const int iwave  = (c.wave  == surfdisp::WaveType::Rayleigh) ? 2 : 1;
        const int igr    = (c.vtype == surfdisp::VelocityType::Group) ? 1 : 0;
        const int imode  = c.mode + 1;   /* Fortran 1-based */

        std::vector<double> cg_cpp(NPER), cg_f(NPER);

        /* ---- warm-up ---- */
        surfdisp::dispersion(MODEL, periods, c.wave, c.mode, c.vtype);
        surfdisp96_f(flat.thk.data(), flat.vp.data(), flat.vs.data(),
                     flat.rho.data(), nlayer, iflsph, iwave, imode, igr,
                     kmax, periods.data(), cg_f.data());

        /* ---- time C++ ---- */
        auto t0 = Clock::now();
        for (int r = 0; r < NREP; r++)
            surfdisp::dispersion(MODEL, periods, c.wave, c.mode, c.vtype);
        double ms_cpp = Ms(Clock::now() - t0).count();

        /* ---- time Fortran ---- */
        t0 = Clock::now();
        for (int r = 0; r < NREP; r++)
            surfdisp96_f(flat.thk.data(), flat.vp.data(), flat.vs.data(),
                         flat.rho.data(), nlayer, iflsph, iwave, imode, igr,
                         kmax, periods.data(), cg_f.data());
        double ms_f90 = Ms(Clock::now() - t0).count();

        printf("%-28s  %10.2f  %10.2f  %8.3f\n",
               c.label, ms_cpp, ms_f90, ms_f90 / ms_cpp);
    }

    /* Total */
    printf("%s\n", std::string(64, '-').c_str());
    printf("Times are total for %d calls. Lower is better.\n", NREP);
    return 0;
}

/*
 * benchmark_slegn.cpp
 * Wall-clock timing: slegn::slegn96 / slegn::slegnpu (C++)
 *                 vs Fortran slegn96_ / slegnpu_
 *
 * Each implementation is called NREP times at several periods.
 * The first call is a warm-up and excluded from timing.
 */

#include "slegn.h"
#include "surfdisp.h"
#include <chrono>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

using Clock = std::chrono::steady_clock;
using Ms    = std::chrono::duration<double, std::milli>;

static const std::vector<surfdisp::Layer> MODEL = {
    {  3.0f,  4.0f, 2.3f, 2.3f },
    {  7.0f,  5.8f, 3.3f, 2.7f },
    { 12.0f,  6.3f, 3.6f, 2.8f },
    { 15.0f,  6.8f, 3.9f, 2.9f },
    { 10.0f,  7.2f, 4.1f, 3.1f },
    { 10.0f,  7.5f, 4.3f, 3.2f },
    { 20.0f,  7.8f, 4.4f, 3.3f },
    { 30.0f,  8.0f, 4.5f, 3.4f },
    { 60.0f,  8.1f, 4.6f, 3.4f },
    {  0.0f,  8.2f, 4.7f, 3.4f },
};

static void unpack_model(const std::vector<surfdisp::Layer>& m,
                         std::vector<float>& thk,
                         std::vector<float>& vs,
                         std::vector<float>& rho)
{
    int n = static_cast<int>(m.size());
    thk.resize(n); vs.resize(n); rho.resize(n);
    for (int i = 0; i < n; i++) {
        thk[i] = m[i].thickness;
        vs[i]  = m[i].vs;
        rho[i] = m[i].density;
    }
}

/* Log-spaced periods from tmin to tmax */
static std::vector<double> log_periods(int n, double tmin, double tmax)
{
    std::vector<double> t(n);
    for (int i = 0; i < n; i++)
        t[i] = tmin * std::pow(tmax / tmin, (double)i / (n - 1));
    return t;
}

int main()
{
    using namespace surfdisp;

    const int NPERIODS = 50;
    const int NREP     = 5000;

    auto periods = log_periods(NPERIODS, 5.0, 100.0);

    int nlayer = static_cast<int>(MODEL.size());
    int iflsph = 0;

    std::vector<float> thk, vs, rho;
    unpack_model(MODEL, thk, vs, rho);

    /* Pre-compute phase velocities for all periods (not counted in timing) */
    std::vector<double> cp(NPERIODS), cp1(NPERIODS), cp2(NPERIODS);
    std::vector<double> t1(NPERIODS), t2(NPERIODS);
    for (int i = 0; i < NPERIODS; i++) {
        double T = periods[i];
        double dT = T * 0.05;   // 5 % bracket for group kernel
        t1[i] = T - dT;
        t2[i] = T + dT;
    }
    {
        std::vector<double> all_t;
        for (int i = 0; i < NPERIODS; i++) {
            all_t.push_back(periods[i]);
            all_t.push_back(t1[i]);
            all_t.push_back(t2[i]);
        }
        auto all_cp = dispersion(MODEL, all_t, WaveType::Love, 0,
                                 VelocityType::Phase, EarthModel::Flat);
        for (int i = 0; i < NPERIODS; i++) {
            cp[i]  = all_cp[3*i];
            cp1[i] = all_cp[3*i+1];
            cp2[i] = all_cp[3*i+2];
        }
    }

    /* Scratch arrays for Fortran calls */
    std::vector<double> f_disp(nlayer), f_stress(nlayer);
    std::vector<double> f_dc2db(nlayer), f_dc2dh(nlayer), f_dc2dr(nlayer);
    std::vector<double> f_du2db(nlayer), f_du2dh(nlayer), f_du2dr(nlayer);

    printf("Benchmark: %d periods × %d repetitions each\n\n", NPERIODS, NREP);
    printf("%-16s  %10s  %10s  %8s\n", "Case", "C++ (ms)", "F90 (ms)", "ratio F/C");
    printf("%s\n", std::string(52, '-').c_str());

    /* ==============================================================
     * slegn96
     * ============================================================== */
    {
        /* warm-up */
        slegn::slegn96(MODEL, periods[0], cp[0], EarthModel::Flat);
        {
            double f_cp = cp[0], f_cg = 0.0;
            slegn96_(thk.data(), vs.data(), rho.data(), nlayer,
                     &periods[0], &f_cp, &f_cg,
                     f_disp.data(), f_stress.data(),
                     f_dc2db.data(), f_dc2dh.data(), f_dc2dr.data(),
                     iflsph);
        }

        /* C++ */
        auto t0 = Clock::now();
        for (int r = 0; r < NREP; r++)
            for (int i = 0; i < NPERIODS; i++)
                slegn::slegn96(MODEL, periods[i], cp[i], EarthModel::Flat);
        double ms_cpp = Ms(Clock::now() - t0).count();

        /* Fortran */
        t0 = Clock::now();
        for (int r = 0; r < NREP; r++)
            for (int i = 0; i < NPERIODS; i++) {
                double f_cp = cp[i], f_cg = 0.0;
                slegn96_(thk.data(), vs.data(), rho.data(), nlayer,
                         &periods[i], &f_cp, &f_cg,
                         f_disp.data(), f_stress.data(),
                         f_dc2db.data(), f_dc2dh.data(), f_dc2dr.data(),
                         iflsph);
            }
        double ms_f90 = Ms(Clock::now() - t0).count();

        printf("%-16s  %10.2f  %10.2f  %8.3f\n",
               "slegn96", ms_cpp, ms_f90, ms_f90 / ms_cpp);
    }

    /* ==============================================================
     * slegnpu
     * ============================================================== */
    {
        /* warm-up */
        slegn::slegnpu(MODEL, periods[0], cp[0], t1[0], cp1[0], t2[0], cp2[0],
                       EarthModel::Flat);
        {
            double f_cp = cp[0], f_cg = 0.0;
            double f_t1 = t1[0], f_cp1 = cp1[0];
            double f_t2 = t2[0], f_cp2 = cp2[0];
            slegnpu_(thk.data(), vs.data(), rho.data(), nlayer,
                     &periods[0], &f_cp, &f_cg,
                     f_disp.data(), f_stress.data(),
                     &f_t1, &f_cp1, &f_t2, &f_cp2,
                     f_dc2db.data(), f_dc2dh.data(), f_dc2dr.data(),
                     f_du2db.data(), f_du2dh.data(), f_du2dr.data(),
                     iflsph);
        }

        /* C++ */
        auto t0 = Clock::now();
        for (int r = 0; r < NREP; r++)
            for (int i = 0; i < NPERIODS; i++)
                slegn::slegnpu(MODEL, periods[i], cp[i],
                               t1[i], cp1[i], t2[i], cp2[i],
                               EarthModel::Flat);
        double ms_cpp = Ms(Clock::now() - t0).count();

        /* Fortran */
        t0 = Clock::now();
        for (int r = 0; r < NREP; r++)
            for (int i = 0; i < NPERIODS; i++) {
                double f_cp = cp[i], f_cg = 0.0;
                double f_t1 = t1[i], f_cp1 = cp1[i];
                double f_t2 = t2[i], f_cp2 = cp2[i];
                slegnpu_(thk.data(), vs.data(), rho.data(), nlayer,
                         &periods[i], &f_cp, &f_cg,
                         f_disp.data(), f_stress.data(),
                         &f_t1, &f_cp1, &f_t2, &f_cp2,
                         f_dc2db.data(), f_dc2dh.data(), f_dc2dr.data(),
                         f_du2db.data(), f_du2dh.data(), f_du2dr.data(),
                         iflsph);
            }
        double ms_f90 = Ms(Clock::now() - t0).count();

        printf("%-16s  %10.2f  %10.2f  %8.3f\n",
               "slegnpu", ms_cpp, ms_f90, ms_f90 / ms_cpp);
    }

    printf("%s\n", std::string(52, '-').c_str());
    printf("Times are total for %d × %d calls. Lower is better.\n",
           NREP, NPERIODS);
    return 0;
}

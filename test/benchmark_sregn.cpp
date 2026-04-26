/*
 * benchmark_sregn.cpp
 * Benchmark sregn::sregn96 / sregn::sregnpu (C++) vs Fortran at 50 periods.
 */

#include "sregn.h"
#include "surfdisp.h"
#include <cstdio>
#include <ctime>
#include <vector>

static const std::vector<surfdisp::Layer> MODEL = {
    {  3.0f,   4.0f, 2.3f, 2.3f },
    {  7.0f,   5.8f, 3.3f, 2.7f },
    { 12.0f,   6.3f, 3.6f, 2.8f },
    { 15.0f,   6.8f, 3.9f, 2.9f },
    { 10.0f,   7.2f, 4.1f, 3.1f },
    { 10.0f,   7.5f, 4.3f, 3.2f },
    { 20.0f,   7.8f, 4.4f, 3.3f },
    { 30.0f,   8.0f, 4.5f, 3.4f },
    { 60.0f,   8.1f, 4.6f, 3.4f },
    {  0.0f,   8.2f, 4.7f, 3.4f },
};

static const int NPERIODS = 50;
static const int NREP     = 5000;

int main()
{
    using namespace surfdisp;

    int nlayer = static_cast<int>(MODEL.size());

    std::vector<double> periods(NPERIODS);
    for (int i = 0; i < NPERIODS; i++)
        periods[i] = 1.0 + i * 1.0;

    std::vector<double> cph = dispersion(MODEL, periods, WaveType::Rayleigh, 0,
                                         VelocityType::Phase, EarthModel::Flat);

    std::vector<float> thk(nlayer), vp(nlayer), vs(nlayer), rho(nlayer);
    for (int i = 0; i < nlayer; i++) {
        thk[i] = MODEL[i].thickness;
        vp[i]  = MODEL[i].vp;
        vs[i]  = MODEL[i].vs;
        rho[i] = MODEL[i].density;
    }

    int iflsph = 0;
    std::vector<double> f_dispu(nlayer), f_dispw(nlayer);
    std::vector<double> f_stressu(nlayer), f_stressw(nlayer);
    std::vector<double> f_dc2da(nlayer), f_dc2db(nlayer);
    std::vector<double> f_dc2dh(nlayer), f_dc2dr(nlayer);
    double dummy_cg;

    /* ----------------------------------------------------------------
     * sregn96
     * ---------------------------------------------------------------- */
    {
        clock_t t0 = clock();
        for (int r = 0; r < NREP; r++) {
            for (int k = 0; k < NPERIODS; k++) {
                double T = periods[k];
                double cp = cph[k];
                double f_cp = cp;
                sregn96_(thk.data(), vp.data(), vs.data(), rho.data(),
                         nlayer, &T, &f_cp, &dummy_cg,
                         f_dispu.data(), f_dispw.data(),
                         f_stressu.data(), f_stressw.data(),
                         f_dc2da.data(), f_dc2db.data(),
                         f_dc2dh.data(), f_dc2dr.data(),
                         iflsph);
            }
        }
        clock_t t1 = clock();
        double ms_f = 1000.0 * (t1 - t0) / CLOCKS_PER_SEC;

        t0 = clock();
        for (int r = 0; r < NREP; r++) {
            for (int k = 0; k < NPERIODS; k++) {
                volatile auto res = sregn::sregn96(MODEL, periods[k], cph[k],
                                                   EarthModel::Flat);
                (void)res;
            }
        }
        t1 = clock();
        double ms_c = 1000.0 * (t1 - t0) / CLOCKS_PER_SEC;

        printf("sregn96  (%d periods x %d reps):\n", NPERIODS, NREP);
        printf("  Fortran: %.1f ms  C++: %.1f ms  speedup: %.2fx\n",
               ms_f, ms_c, ms_f / ms_c);
    }

    /* ----------------------------------------------------------------
     * sregnpu
     * ---------------------------------------------------------------- */
    {
        std::vector<double> periods2(NPERIODS), periods3(NPERIODS);
        std::vector<double> cph2(NPERIODS), cph3(NPERIODS);
        for (int i = 0; i < NPERIODS; i++) {
            periods2[i] = periods[i] - 0.5;
            periods3[i] = periods[i] + 0.5;
        }
        cph2 = dispersion(MODEL, periods2, WaveType::Rayleigh, 0, VelocityType::Phase);
        cph3 = dispersion(MODEL, periods3, WaveType::Rayleigh, 0, VelocityType::Phase);

        std::vector<double> f_du2da(nlayer), f_du2db(nlayer);
        std::vector<double> f_du2dh(nlayer), f_du2dr(nlayer);

        clock_t t0 = clock();
        for (int r = 0; r < NREP; r++) {
            for (int k = 0; k < NPERIODS; k++) {
                double T = periods[k], t1 = periods2[k], t2 = periods3[k];
                double cp = cph[k], cp1 = cph2[k], cp2 = cph3[k];
                double f_cp = cp, f_cg = 0.0;
                double f_t1 = t1, f_cp1 = cp1, f_t2 = t2, f_cp2 = cp2;
                sregnpu_(thk.data(), vp.data(), vs.data(), rho.data(),
                         nlayer, &T, &f_cp, &f_cg,
                         f_dispu.data(), f_dispw.data(),
                         f_stressu.data(), f_stressw.data(),
                         &f_t1, &f_cp1, &f_t2, &f_cp2,
                         f_dc2da.data(), f_dc2db.data(),
                         f_dc2dh.data(), f_dc2dr.data(),
                         f_du2da.data(), f_du2db.data(),
                         f_du2dh.data(), f_du2dr.data(),
                         iflsph);
            }
        }
        clock_t t1_c = clock();
        double ms_f = 1000.0 * (t1_c - t0) / CLOCKS_PER_SEC;

        t0 = clock();
        for (int r = 0; r < NREP; r++) {
            for (int k = 0; k < NPERIODS; k++) {
                volatile auto res = sregn::sregnpu(MODEL,
                    periods[k], cph[k],
                    periods2[k], cph2[k],
                    periods3[k], cph3[k],
                    EarthModel::Flat);
                (void)res;
            }
        }
        t1_c = clock();
        double ms_c = 1000.0 * (t1_c - t0) / CLOCKS_PER_SEC;

        printf("sregnpu  (%d periods x %d reps):\n", NPERIODS, NREP);
        printf("  Fortran: %.1f ms  C++: %.1f ms  speedup: %.2fx\n",
               ms_f, ms_c, ms_f / ms_c);
    }

    return 0;
}

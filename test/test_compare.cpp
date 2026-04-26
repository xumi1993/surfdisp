/*
 * test_compare.cpp
 * Compare surfdisp::dispersion (C++) against the original Fortran surfdisp96
 * using a 10-layer crustal/mantle model.
 */

#include "surfdisp.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

/* ------------------------------------------------------------------ */
/* Earth model: PREM-like crust + upper mantle                        */
/* ------------------------------------------------------------------ */
static const std::vector<surfdisp::Layer> MODEL = {
    /*  thickness   vp    vs   density  */
    {  3.0f,   4.0f, 2.3f, 2.3f },
    {  7.0f,   5.8f, 3.3f, 2.7f },
    { 12.0f,   6.3f, 3.6f, 2.8f },
    { 15.0f,   6.8f, 3.9f, 2.9f },
    { 10.0f,   7.2f, 4.1f, 3.1f },
    { 10.0f,   7.5f, 4.3f, 3.2f },
    { 20.0f,   7.8f, 4.4f, 3.3f },
    { 30.0f,   8.0f, 4.5f, 3.4f },
    { 60.0f,   8.1f, 4.6f, 3.4f },
    {  0.0f,   8.2f, 4.7f, 3.4f },  /* halfspace: thickness = 0 */
};

/* 50 log-spaced periods 1–100 s */
static std::vector<double> make_periods(int n, double tmin, double tmax)
{
    std::vector<double> t(n);
    for (int i = 0; i < n; i++)
        t[i] = tmin * std::pow(tmax / tmin, (double)i / (n - 1));
    return t;
}

/* ------------------------------------------------------------------ */
/* Helpers                                                             */
/* ------------------------------------------------------------------ */
static void write_csv(const char *fname,
                      const std::vector<double>& t,
                      const std::vector<double>& cpp,
                      const std::vector<double>& f90)
{
    FILE *fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }
    fprintf(fp, "period_s,cpp,fortran\n");
    for (size_t i = 0; i < t.size(); i++)
        fprintf(fp, "%.6f,%.8f,%.8f\n", t[i], cpp[i], f90[i]);
    fclose(fp);
    printf("  wrote %s\n", fname);
}

static void compare(const char *label, const char *csv,
                    surfdisp::WaveType   wave,
                    int                  mode,   /* 0-based */
                    surfdisp::VelocityType vtype,
                    const std::vector<double>& periods)
{
    using namespace surfdisp;

    /* New C++ API */
    auto cg_cpp = dispersion(MODEL, periods, wave, mode, vtype);

    /* Fortran reference via legacy wrapper */
    const int nlayer  = static_cast<int>(MODEL.size());
    const int kmax    = static_cast<int>(periods.size());
    const int iflsph  = 0;
    const int iwave   = (wave  == WaveType::Rayleigh)       ? 2 : 1;
    const int igr     = (vtype == VelocityType::Group)      ? 1 : 0;
    const int imode_f = mode + 1;   /* Fortran 1-based */

    std::vector<float> thk(nlayer), vp(nlayer), vs(nlayer), rho(nlayer);
    for (int i = 0; i < nlayer; i++) {
        thk[i] = MODEL[i].thickness;
        vp[i]  = MODEL[i].vp;
        vs[i]  = MODEL[i].vs;
        rho[i] = MODEL[i].density;
    }

    std::vector<double> cg_f(kmax, 0.0);
    surfdisp96_f(thk.data(), vp.data(), vs.data(), rho.data(),
                 nlayer, iflsph, iwave, imode_f, igr, kmax,
                 periods.data(), cg_f.data());

    /* Statistics */
    double max_rel = 0.0, sum_rel = 0.0;
    int cnt = 0;
    for (int i = 0; i < kmax; i++) {
        if (cg_cpp[i] == 0.0 && cg_f[i] == 0.0) continue;
        double ref = std::fabs(cg_f[i]) > 1e-10 ? cg_f[i] : cg_cpp[i];
        if (std::fabs(ref) < 1e-10) continue;
        double rel = std::fabs(cg_cpp[i] - cg_f[i]) / std::fabs(ref);
        if (rel > max_rel) max_rel = rel;
        sum_rel += rel;
        cnt++;
    }
    double mean_rel = cnt > 0 ? sum_rel / cnt : 0.0;

    printf("%-42s  max_rel=%.2e  mean_rel=%.2e  [%s]\n",
           label, max_rel, mean_rel, max_rel < 1e-4 ? "PASS" : "WARN");

    write_csv(csv, periods, cg_cpp, cg_f);
}

int main()
{
    using namespace surfdisp;

    auto periods = make_periods(50, 1.0, 100.0);

    printf("=== surfdisp: C++ (new API) vs Fortran comparison ===\n\n");

    compare("Love  mode=0 (fundamental)  phase", "love_fund_phase.csv",
            WaveType::Love, 0, VelocityType::Phase, periods);

    compare("Love  mode=0 (fundamental)  group", "love_fund_group.csv",
            WaveType::Love, 0, VelocityType::Group, periods);

    compare("Love  mode=1 (1st overtone) phase", "love_1st_phase.csv",
            WaveType::Love, 1, VelocityType::Phase, periods);

    compare("Rayleigh mode=0 (fundamental) phase", "ray_fund_phase.csv",
            WaveType::Rayleigh, 0, VelocityType::Phase, periods);

    compare("Rayleigh mode=0 (fundamental) group", "ray_fund_group.csv",
            WaveType::Rayleigh, 0, VelocityType::Group, periods);

    compare("Rayleigh mode=1 (1st overtone) phase", "ray_1st_phase.csv",
            WaveType::Rayleigh, 1, VelocityType::Phase, periods);

    printf("\nDone.\n");
    return 0;
}

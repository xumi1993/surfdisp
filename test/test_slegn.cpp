/*
 * test_slegn.cpp
 * Compare slegn::slegn96 / slegn::slegnpu (C++) against the original Fortran
 * slegn96_ / slegnpu_ using a 10-layer crustal/mantle model.
 */

#include "slegn.h"
#include "surfdisp.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

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
    {  0.0f,   8.2f, 4.7f, 3.4f },
};

/* ------------------------------------------------------------------ */
/* Helper: relative error                                              */
/* ------------------------------------------------------------------ */
static double rel_err(double a, double b)
{
    double ref = std::fabs(b) > 1.0e-10 ? std::fabs(b) : std::fabs(a);
    if (ref < 1.0e-30) return 0.0;
    return std::fabs(a - b) / ref;
}

/* ------------------------------------------------------------------ */
/* Write comparison CSV                                                 */
/* ------------------------------------------------------------------ */
static void write_csv(const char *fname,
                      const std::vector<double>& depth_km,
                      const std::vector<double>& cpp_vals,
                      const std::vector<double>& fort_vals)
{
    FILE *fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }
    fprintf(fp, "depth_km,cpp_value,fortran_value\n");
    for (size_t i = 0; i < depth_km.size(); i++)
        fprintf(fp, "%.4f,%.10e,%.10e\n", depth_km[i], cpp_vals[i], fort_vals[i]);
    fclose(fp);
    printf("  wrote %s\n", fname);
}

/* ------------------------------------------------------------------ */
/* Build cumulative depth array (midpoint of each layer)               */
/* ------------------------------------------------------------------ */
static std::vector<double> make_depths()
{
    std::vector<double> d;
    double cum = 0.0;
    for (size_t i = 0; i < MODEL.size(); i++) {
        d.push_back(cum + MODEL[i].thickness * 0.5);
        cum += MODEL[i].thickness;
    }
    return d;
}

/* ------------------------------------------------------------------ */
/* Max relative error over all layers                                   */
/* ------------------------------------------------------------------ */
static double max_rel_err(const std::vector<double>& a,
                          const std::vector<double>& b)
{
    double mx = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        double e = rel_err(a[i], b[i]);
        if (e > mx) mx = e;
    }
    return mx;
}

/* ------------------------------------------------------------------ */
/* Build flat model arrays for Fortran call                             */
/* ------------------------------------------------------------------ */
static void unpack_model(std::vector<float>& thk,
                         std::vector<float>& vs,
                         std::vector<float>& rho)
{
    int n = static_cast<int>(MODEL.size());
    thk.resize(n); vs.resize(n); rho.resize(n);
    for (int i = 0; i < n; i++) {
        thk[i] = MODEL[i].thickness;
        vs[i]  = MODEL[i].vs;
        rho[i] = MODEL[i].density;
    }
}

/* ================================================================== */
int main()
{
    using namespace surfdisp;

    int nlayer = static_cast<int>(MODEL.size());
    int iflsph = 0;   // flat earth

    std::vector<float> thk, vs, rho;
    unpack_model(thk, vs, rho);

    auto depths = make_depths();

    printf("=== slegn C++ vs Fortran comparison ===\n\n");

    /* ----------------------------------------------------------------
     * slegn96 test at T=10 s
     * ---------------------------------------------------------------- */
    {
        double T = 10.0;

        // Get phase velocity from surfdisp::dispersion
        std::vector<double> pv = dispersion(MODEL, {T}, WaveType::Love, 0,
                                            VelocityType::Phase, EarthModel::Flat);
        double cp_init = pv[0];
        printf("slegn96 test: T=%.1f s, cp=%.6f km/s\n", T, cp_init);

        /* --- Fortran call --- */
        std::vector<double> f_disp(nlayer), f_stress(nlayer);
        std::vector<double> f_dc2db(nlayer), f_dc2dh(nlayer), f_dc2dr(nlayer);
        double f_cp = cp_init, f_cg = 0.0;

        slegn96_(thk.data(), vs.data(), rho.data(),
                 nlayer,
                 &T, &f_cp, &f_cg,
                 f_disp.data(), f_stress.data(),
                 f_dc2db.data(), f_dc2dh.data(), f_dc2dr.data(),
                 iflsph);

        printf("  Fortran: cp=%.8f cg=%.8f\n", f_cp, f_cg);

        /* --- C++ call --- */
        auto res = slegn::slegn96(MODEL, T, cp_init, EarthModel::Flat);
        printf("  C++    : cp=%.8f cg=%.8f\n", res.cp, res.cg);

        /* --- Scalar comparison --- */
        double cp_err = rel_err(res.cp, f_cp);
        double cg_err = rel_err(res.cg, f_cg);
        printf("  cp_rel_err=%.2e  cg_rel_err=%.2e\n", cp_err, cg_err);

        /* --- Array comparisons --- */
        double err_disp   = max_rel_err(res.disp,   f_disp);
        double err_stress = max_rel_err(res.stress,  f_stress);
        double err_dc2db  = max_rel_err(res.dc2db,   f_dc2db);
        double err_dc2dh  = max_rel_err(res.dc2dh,   f_dc2dh);
        double err_dc2dr  = max_rel_err(res.dc2dr,   f_dc2dr);

        auto thresh = 1.0e-4;
        printf("  disp   max_rel=%.2e [%s]\n", err_disp,   err_disp   < thresh ? "PASS" : "WARN");
        printf("  stress max_rel=%.2e [%s]\n", err_stress, err_stress < thresh ? "PASS" : "WARN");
        printf("  dc2db  max_rel=%.2e [%s]\n", err_dc2db,  err_dc2db  < thresh ? "PASS" : "WARN");
        printf("  dc2dh  max_rel=%.2e [%s]\n", err_dc2dh,  err_dc2dh  < thresh ? "PASS" : "WARN");
        printf("  dc2dr  max_rel=%.2e [%s]\n", err_dc2dr,  err_dc2dr  < thresh ? "PASS" : "WARN");

        /* --- CSV output --- */
        write_csv("slegn96_disp.csv",   depths, res.disp,   f_disp);
        write_csv("slegn96_stress.csv", depths, res.stress,  f_stress);
        write_csv("slegn96_dc2db.csv",  depths, res.dc2db,   f_dc2db);
        write_csv("slegn96_dc2dh.csv",  depths, res.dc2dh,   f_dc2dh);
        write_csv("slegn96_dc2dr.csv",  depths, res.dc2dr,   f_dc2dr);

        /* --- Summary CSV row --- */
        {
            FILE *fp = fopen("slegn96_summary.csv", "w");
            if (fp) {
                fprintf(fp, "quantity,cpp_value,fortran_value\n");
                fprintf(fp, "cp,%.10e,%.10e\n", res.cp, f_cp);
                fprintf(fp, "cg,%.10e,%.10e\n", res.cg, f_cg);
                fclose(fp);
                printf("  wrote slegn96_summary.csv\n");
            }
        }
    }

    printf("\n");

    /* ----------------------------------------------------------------
     * slegnpu test at T=20 s, t1=19 s, t2=21 s
     * ---------------------------------------------------------------- */
    {
        double T  = 20.0;
        double t1 = 19.0;
        double t2 = 21.0;

        // Get phase velocities
        std::vector<double> pv = dispersion(MODEL, {T, t1, t2}, WaveType::Love, 0,
                                            VelocityType::Phase, EarthModel::Flat);
        double cp0 = pv[0], cp1_v = pv[1], cp2_v = pv[2];
        printf("slegnpu test: T=%.1f s, cp=%.6f km/s\n", T, cp0);
        printf("              t1=%.1f s, cp1=%.6f km/s\n", t1, cp1_v);
        printf("              t2=%.1f s, cp2=%.6f km/s\n", t2, cp2_v);

        /* --- Fortran call --- */
        std::vector<double> f_disp(nlayer), f_stress(nlayer);
        std::vector<double> f_dc2db(nlayer), f_dc2dh(nlayer), f_dc2dr(nlayer);
        std::vector<double> f_du2db(nlayer), f_du2dh(nlayer), f_du2dr(nlayer);
        double f_cp = cp0, f_cg = 0.0;
        double f_t1 = t1, f_cp1 = cp1_v;
        double f_t2 = t2, f_cp2 = cp2_v;

        slegnpu_(thk.data(), vs.data(), rho.data(),
                 nlayer,
                 &T, &f_cp, &f_cg,
                 f_disp.data(), f_stress.data(),
                 &f_t1, &f_cp1,
                 &f_t2, &f_cp2,
                 f_dc2db.data(), f_dc2dh.data(), f_dc2dr.data(),
                 f_du2db.data(), f_du2dh.data(), f_du2dr.data(),
                 iflsph);

        printf("  Fortran: cp=%.8f cg=%.8f\n", f_cp, f_cg);

        /* --- C++ call --- */
        auto res = slegn::slegnpu(MODEL, T, cp0, t1, cp1_v, t2, cp2_v,
                                  EarthModel::Flat);
        printf("  C++    : cp=%.8f cg=%.8f\n", res.cp, res.cg);

        double cp_err = rel_err(res.cp, f_cp);
        double cg_err = rel_err(res.cg, f_cg);
        printf("  cp_rel_err=%.2e  cg_rel_err=%.2e\n", cp_err, cg_err);

        double err_du2db = max_rel_err(res.du2db, f_du2db);
        double err_du2dh = max_rel_err(res.du2dh, f_du2dh);
        double err_du2dr = max_rel_err(res.du2dr, f_du2dr);

        auto thresh = 1.0e-4;
        printf("  du2db  max_rel=%.2e [%s]\n", err_du2db, err_du2db < thresh ? "PASS" : "WARN");
        printf("  du2dh  max_rel=%.2e [%s]\n", err_du2dh, err_du2dh < thresh ? "PASS" : "WARN");
        printf("  du2dr  max_rel=%.2e [%s]\n", err_du2dr, err_du2dr < thresh ? "PASS" : "WARN");

        /* --- CSV output --- */
        write_csv("slegnpu_du2db.csv", depths, res.du2db, f_du2db);
        write_csv("slegnpu_du2dh.csv", depths, res.du2dh, f_du2dh);
        write_csv("slegnpu_du2dr.csv", depths, res.du2dr, f_du2dr);

        /* --- Summary CSV row --- */
        {
            FILE *fp = fopen("slegnpu_summary.csv", "w");
            if (fp) {
                fprintf(fp, "quantity,cpp_value,fortran_value\n");
                fprintf(fp, "cp,%.10e,%.10e\n", res.cp, f_cp);
                fprintf(fp, "cg,%.10e,%.10e\n", res.cg, f_cg);
                fclose(fp);
                printf("  wrote slegnpu_summary.csv\n");
            }
        }
    }

    printf("\nDone.\n");
    return 0;
}

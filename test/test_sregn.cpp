/*
 * test_sregn.cpp
 * Compare sregn::sregn96 / sregn::sregn96_hti / sregn::sregnpu (C++)
 * against the original Fortran sregn96_ / sregn96_hti_ / sregnpu_.
 */

#include "sregn.h"
#include "surfdisp.h"
#include <cstdio>
#include <cmath>
#include <vector>

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

static double rel_err(double a, double b)
{
    double ref = std::fabs(b) > 1.0e-10 ? std::fabs(b) : std::fabs(a);
    if (ref < 1.0e-30) return 0.0;
    return std::fabs(a - b) / ref;
}

static void write_csv(const char* fname,
                      const std::vector<double>& cpp_vals,
                      const std::vector<double>& fort_vals)
{
    FILE* fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }
    fprintf(fp, "layer,cpp_value,fortran_value\n");
    for (size_t i = 0; i < cpp_vals.size(); i++)
        fprintf(fp, "%zu,%.10e,%.10e\n", i, cpp_vals[i], fort_vals[i]);
    fclose(fp);
    printf("  wrote %s\n", fname);
}

static double max_rel_err(const std::vector<double>& a, const std::vector<double>& b)
{
    double mx = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        double e = rel_err(a[i], b[i]);
        if (e > mx) mx = e;
    }
    return mx;
}

static void unpack_model(std::vector<float>& thk, std::vector<float>& vp,
                         std::vector<float>& vs,  std::vector<float>& rho)
{
    int n = static_cast<int>(MODEL.size());
    thk.resize(n); vp.resize(n); vs.resize(n); rho.resize(n);
    for (int i = 0; i < n; i++) {
        thk[i] = MODEL[i].thickness;
        vp[i]  = MODEL[i].vp;
        vs[i]  = MODEL[i].vs;
        rho[i] = MODEL[i].density;
    }
}

int main()
{
    using namespace surfdisp;

    int nlayer = static_cast<int>(MODEL.size());
    int iflsph = 0;

    std::vector<float> thk, vp, vs, rho;
    unpack_model(thk, vp, vs, rho);

    printf("=== sregn C++ vs Fortran comparison ===\n\n");

    const double thresh = 1.0e-4;

    /* ----------------------------------------------------------------
     * sregn96 test at T=10 s
     * ---------------------------------------------------------------- */
    {
        double T = 10.0;
        std::vector<double> pv = dispersion(MODEL, {T}, WaveType::Rayleigh, 0,
                                            VelocityType::Phase, EarthModel::Flat);
        double cp_init = pv[0];
        printf("sregn96 test: T=%.1f s, cp=%.6f km/s\n", T, cp_init);

        std::vector<double> f_dispu(nlayer), f_dispw(nlayer);
        std::vector<double> f_stressu(nlayer), f_stressw(nlayer);
        std::vector<double> f_dc2da(nlayer), f_dc2db(nlayer);
        std::vector<double> f_dc2dh(nlayer), f_dc2dr(nlayer);
        double f_cp = cp_init, f_cg = 0.0;

        sregn96_(thk.data(), vp.data(), vs.data(), rho.data(),
                 nlayer, &T, &f_cp, &f_cg,
                 f_dispu.data(), f_dispw.data(),
                 f_stressu.data(), f_stressw.data(),
                 f_dc2da.data(), f_dc2db.data(),
                 f_dc2dh.data(), f_dc2dr.data(),
                 iflsph);
        printf("  Fortran: cp=%.8f  cg=%.8f\n", f_cp, f_cg);

        auto res = sregn::sregn96(MODEL, T, cp_init, EarthModel::Flat);
        printf("  C++    : cp=%.8f  cg=%.8f\n", res.cp, res.cg);

        printf("  cp_rel_err=%.2e  cg_rel_err=%.2e\n",
               rel_err(res.cp, f_cp), rel_err(res.cg, f_cg));

        double e_dispu   = max_rel_err(res.dispu,   f_dispu);
        double e_dispw   = max_rel_err(res.dispw,   f_dispw);
        double e_stressu = max_rel_err(res.stressu, f_stressu);
        double e_stressw = max_rel_err(res.stressw, f_stressw);
        double e_dc2da   = max_rel_err(res.dc2da,   f_dc2da);
        double e_dc2db   = max_rel_err(res.dc2db,   f_dc2db);
        double e_dc2dh   = max_rel_err(res.dc2dh,   f_dc2dh);
        double e_dc2dr   = max_rel_err(res.dc2dr,   f_dc2dr);

        printf("  dispu   max_rel=%.2e [%s]\n", e_dispu,   e_dispu   < thresh ? "PASS" : "WARN");
        printf("  dispw   max_rel=%.2e [%s]\n", e_dispw,   e_dispw   < thresh ? "PASS" : "WARN");
        printf("  stressu max_rel=%.2e [%s]\n", e_stressu, e_stressu < thresh ? "PASS" : "WARN");
        printf("  stressw max_rel=%.2e [%s]\n", e_stressw, e_stressw < thresh ? "PASS" : "WARN");
        printf("  dc2da   max_rel=%.2e [%s]\n", e_dc2da,   e_dc2da   < thresh ? "PASS" : "WARN");
        printf("  dc2db   max_rel=%.2e [%s]\n", e_dc2db,   e_dc2db   < thresh ? "PASS" : "WARN");
        printf("  dc2dh   max_rel=%.2e [%s]\n", e_dc2dh,   e_dc2dh   < thresh ? "PASS" : "WARN");
        printf("  dc2dr   max_rel=%.2e [%s]\n", e_dc2dr,   e_dc2dr   < thresh ? "PASS" : "WARN");

        write_csv("sregn96_dispu.csv",   res.dispu,   f_dispu);
        write_csv("sregn96_dispw.csv",   res.dispw,   f_dispw);
        write_csv("sregn96_dc2da.csv",   res.dc2da,   f_dc2da);
        write_csv("sregn96_dc2db.csv",   res.dc2db,   f_dc2db);
        write_csv("sregn96_dc2dh.csv",   res.dc2dh,   f_dc2dh);
        write_csv("sregn96_dc2dr.csv",   res.dc2dr,   f_dc2dr);

        {
            FILE* fp = fopen("sregn96_summary.csv", "w");
            if (fp) {
                fprintf(fp, "quantity,cpp_value,fortran_value\n");
                fprintf(fp, "cp,%.10e,%.10e\n", res.cp, f_cp);
                fprintf(fp, "cg,%.10e,%.10e\n", res.cg, f_cg);
                fclose(fp);
                printf("  wrote sregn96_summary.csv\n");
            }
        }
    }

    printf("\n");

    /* ----------------------------------------------------------------
     * sregn96_hti test at T=10 s
     * ---------------------------------------------------------------- */
    {
        double T = 10.0;
        std::vector<double> pv = dispersion(MODEL, {T}, WaveType::Rayleigh, 0,
                                            VelocityType::Phase, EarthModel::Flat);
        double cp_init = pv[0];
        printf("sregn96_hti test: T=%.1f s, cp=%.6f km/s\n", T, cp_init);

        std::vector<double> f_dispu(nlayer), f_dispw(nlayer);
        std::vector<double> f_stressu(nlayer), f_stressw(nlayer);
        std::vector<double> f_dc2da(nlayer), f_dc2db(nlayer);
        std::vector<double> f_dc2dh(nlayer), f_dc2dr(nlayer);
        std::vector<double> f_dc2dgc(nlayer), f_dc2dgs(nlayer);
        double f_cp = cp_init, f_cg = 0.0;

        sregn96_hti_(thk.data(), vp.data(), vs.data(), rho.data(),
                     nlayer, &T, &f_cp, &f_cg,
                     f_dispu.data(), f_dispw.data(),
                     f_stressu.data(), f_stressw.data(),
                     f_dc2da.data(), f_dc2db.data(),
                     f_dc2dh.data(), f_dc2dr.data(),
                     f_dc2dgc.data(), f_dc2dgs.data(),
                     iflsph);
        printf("  Fortran: cp=%.8f  cg=%.8f\n", f_cp, f_cg);

        auto res = sregn::sregn96_hti(MODEL, T, cp_init, EarthModel::Flat);
        printf("  C++    : cp=%.8f  cg=%.8f\n", res.cp, res.cg);

        printf("  cp_rel_err=%.2e  cg_rel_err=%.2e\n",
               rel_err(res.cp, f_cp), rel_err(res.cg, f_cg));

        double e_dc2dgc = max_rel_err(res.dc2dgc, f_dc2dgc);
        double e_dc2dgs = max_rel_err(res.dc2dgs, f_dc2dgs);
        printf("  dc2dgc  max_rel=%.2e [%s]\n", e_dc2dgc, e_dc2dgc < thresh ? "PASS" : "WARN");
        printf("  dc2dgs  max_rel=%.2e [%s]\n", e_dc2dgs, e_dc2dgs < thresh ? "PASS" : "WARN");

        write_csv("sregn96_hti_dc2dgc.csv", res.dc2dgc, f_dc2dgc);
        write_csv("sregn96_hti_dc2dgs.csv", res.dc2dgs, f_dc2dgs);
    }

    printf("\n");

    /* ----------------------------------------------------------------
     * sregnpu test at T=20 s, t1=19 s, t2=21 s
     * ---------------------------------------------------------------- */
    {
        double T  = 20.0;
        double t1 = 19.0;
        double t2 = 21.0;

        std::vector<double> pv = dispersion(MODEL, {T, t1, t2}, WaveType::Rayleigh, 0,
                                            VelocityType::Phase, EarthModel::Flat);
        double cp0 = pv[0], cp1_v = pv[1], cp2_v = pv[2];
        printf("sregnpu test: T=%.1f s, cp=%.6f km/s\n", T, cp0);
        printf("              t1=%.1f s, cp1=%.6f km/s\n", t1, cp1_v);
        printf("              t2=%.1f s, cp2=%.6f km/s\n", t2, cp2_v);

        std::vector<double> f_dispu(nlayer), f_dispw(nlayer);
        std::vector<double> f_stressu(nlayer), f_stressw(nlayer);
        std::vector<double> f_dc2da(nlayer), f_dc2db(nlayer);
        std::vector<double> f_dc2dh(nlayer), f_dc2dr(nlayer);
        std::vector<double> f_du2da(nlayer), f_du2db(nlayer);
        std::vector<double> f_du2dh(nlayer), f_du2dr(nlayer);
        double f_cp = cp0, f_cg = 0.0;
        double f_t1 = t1, f_cp1 = cp1_v;
        double f_t2 = t2, f_cp2 = cp2_v;

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
        printf("  Fortran: cp=%.8f  cg=%.8f\n", f_cp, f_cg);

        auto res = sregn::sregnpu(MODEL, T, cp0, t1, cp1_v, t2, cp2_v,
                                  EarthModel::Flat);
        printf("  C++    : cp=%.8f  cg=%.8f\n", res.cp, res.cg);

        printf("  cp_rel_err=%.2e  cg_rel_err=%.2e\n",
               rel_err(res.cp, f_cp), rel_err(res.cg, f_cg));

        double e_du2da = max_rel_err(res.du2da, f_du2da);
        double e_du2db = max_rel_err(res.du2db, f_du2db);
        double e_du2dh = max_rel_err(res.du2dh, f_du2dh);
        double e_du2dr = max_rel_err(res.du2dr, f_du2dr);
        double e_dc2da = max_rel_err(res.dc2da, f_dc2da);
        double e_dc2db = max_rel_err(res.dc2db, f_dc2db);
        double e_dc2dh = max_rel_err(res.dc2dh, f_dc2dh);
        double e_dc2dr = max_rel_err(res.dc2dr, f_dc2dr);

        printf("  dc2da   max_rel=%.2e [%s]\n", e_dc2da, e_dc2da < thresh ? "PASS" : "WARN");
        printf("  dc2db   max_rel=%.2e [%s]\n", e_dc2db, e_dc2db < thresh ? "PASS" : "WARN");
        printf("  dc2dh   max_rel=%.2e [%s]\n", e_dc2dh, e_dc2dh < thresh ? "PASS" : "WARN");
        printf("  dc2dr   max_rel=%.2e [%s]\n", e_dc2dr, e_dc2dr < thresh ? "PASS" : "WARN");
        printf("  du2da   max_rel=%.2e [%s]\n", e_du2da, e_du2da < thresh ? "PASS" : "WARN");
        printf("  du2db   max_rel=%.2e [%s]\n", e_du2db, e_du2db < thresh ? "PASS" : "WARN");
        printf("  du2dh   max_rel=%.2e [%s]\n", e_du2dh, e_du2dh < thresh ? "PASS" : "WARN");
        printf("  du2dr   max_rel=%.2e [%s]\n", e_du2dr, e_du2dr < thresh ? "PASS" : "WARN");

        write_csv("sregnpu_du2da.csv", res.du2da, f_du2da);
        write_csv("sregnpu_du2db.csv", res.du2db, f_du2db);
        write_csv("sregnpu_du2dh.csv", res.du2dh, f_du2dh);
        write_csv("sregnpu_du2dr.csv", res.du2dr, f_du2dr);

        {
            FILE* fp = fopen("sregnpu_summary.csv", "w");
            if (fp) {
                fprintf(fp, "quantity,cpp_value,fortran_value\n");
                fprintf(fp, "cp,%.10e,%.10e\n", res.cp, f_cp);
                fprintf(fp, "cg,%.10e,%.10e\n", res.cg, f_cg);
                fclose(fp);
                printf("  wrote sregnpu_summary.csv\n");
            }
        }
    }

    printf("\nDone.\n");
    return 0;
}

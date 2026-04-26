# surfdisp96 — C++ Port

C++ port of the classic Fortran surface-wave dispersion code **surfdisp96**
(Herrmann & Wang 1971–2023, Saint Louis University).  
Computes **Love** and **Rayleigh** wave **phase / group velocities** for any
1-D layered earth model.

---

## Directory layout

```
surfdisp/
├── CMakeLists.txt
├── include/
│   └── surfdisp.h              C++ public API + legacy C interface
├── src/
│   ├── surfdisp96.f            Original Fortran (unmodified)
│   ├── surfdisp_fwrapper.f90   bind(C) wrapper so C++ can call Fortran
│   └── surfdisp.cpp            C++ port + surfdisp::dispersion()
├── test/
│   ├── test_compare.cpp        Numerical accuracy comparison
│   └── benchmark.cpp           Wall-clock timing comparison
└── scripts/
    └── plot_comparison.py      Python plot of C++ vs Fortran curves
```

---

## Build

Requirements: CMake ≥ 3.14, a C++14 compiler, gfortran.

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

---

## C++ API

```cpp
#include "surfdisp.h"

// Define the 1-D velocity model (last layer = halfspace, thickness ignored)
std::vector<surfdisp::Layer> model = {
    //  thickness(km)  vp(km/s)  vs(km/s)  density(g/cm³)
    {  3.0f,   4.0f,  2.3f,  2.3f },
    {  7.0f,   5.8f,  3.3f,  2.7f },
    { 12.0f,   6.3f,  3.6f,  2.8f },
    { 15.0f,   6.8f,  3.9f,  2.9f },
    { 10.0f,   7.2f,  4.1f,  3.1f },
    { 10.0f,   7.5f,  4.3f,  3.2f },
    { 20.0f,   7.8f,  4.4f,  3.3f },
    { 30.0f,   8.0f,  4.5f,  3.4f },
    { 60.0f,   8.1f,  4.6f,  3.4f },
    {  0.0f,   8.2f,  4.7f,  3.4f },  // halfspace
};

// 50 log-spaced periods from 1 to 100 s
std::vector<double> periods = { /* ... */ };

// Rayleigh fundamental-mode phase velocity, flat earth
auto cg = surfdisp::dispersion(
    model, periods,
    surfdisp::WaveType::Rayleigh,      // Love | Rayleigh
    0,                                  // mode: 0 = fundamental, 1 = 1st overtone, …
    surfdisp::VelocityType::Phase,     // Phase | Group
    surfdisp::EarthModel::Flat         // Flat  | Spherical
);
```

### API reference

| Type / function | Description |
|---|---|
| `surfdisp::Layer` | `{thickness, vp, vs, density}` — one model layer |
| `surfdisp::WaveType` | `Love`, `Rayleigh` |
| `surfdisp::VelocityType` | `Phase`, `Group` |
| `surfdisp::EarthModel` | `Flat`, `Spherical` (applies earth-flattening transform) |
| `surfdisp::dispersion(model, periods, wave, mode, vel, earth)` | Returns `std::vector<double>` of velocities (km/s); 0 where a mode was not found |

Throws `std::invalid_argument` for bad input.  
Throws `std::runtime_error` if the fundamental-mode search fails completely.

---

## Accuracy test

Model: 10-layer crust/mantle, 50 log-spaced periods 1–100 s, flat earth.  
Reference: original Fortran `surfdisp96`.

```
$ ./test_compare
```

```
Case                                        max_rel     mean_rel
Love  mode=0 (fundamental)  phase        0.00e+00    0.00e+00   PASS
Love  mode=0 (fundamental)  group        1.43e-05    4.68e-06   PASS
Love  mode=1 (1st overtone) phase        0.00e+00    0.00e+00   PASS
Rayleigh mode=0 (fundamental) phase      0.00e+00    0.00e+00   PASS
Rayleigh mode=0 (fundamental) group      1.19e-05    4.26e-06   PASS
Rayleigh mode=1 (1st overtone) phase     0.00e+00    0.00e+00   PASS
```

- **Phase velocity**: bit-identical to Fortran (`max_rel = 0`).
- **Group velocity**: relative error < 1.5 × 10⁻⁵. The tiny difference is
  intentional — the original Fortran uses `sngl()` (double→float→double
  round-trip) in the finite-difference group-velocity formula; the C++ port
  replicates this exactly, and the residual is floating-point rounding noise.

The test also writes six CSV files (`love_fund_phase.csv`, …) readable by
`scripts/plot_comparison.py`.

---

## Speed benchmark

200 periods × 2000 repetitions, Apple M-series (ARM64), `-O2`.

```
$ ./benchmark
```

```
Case                           C++ (ms)   F90 (ms)   ratio F/C
----------------------------------------------------------------
Love  mode=0 phase               392        452        1.15×
Love  mode=0 group               751        872        1.16×
Love  mode=1 phase               655        821        1.25×
Rayleigh mode=0 phase           1068       1830        1.71×
Rayleigh mode=0 group           2055       3514        1.71×
Rayleigh mode=1 phase           1844       2981        1.62×
```

Times are **total** for 2000 calls; lower is better.

### Optimization-level comparison (Rayleigh mode=0 phase, total 2000 calls)

| Flag | C++ (ms) | Fortran (ms) | ratio F/C |
|---|---|---|---|
| `-O2` | 1068 | 1830 | 1.71× |
| `-O3` | 1079 | 1825 | 1.69× |
| `-O3 -march=native -ffast-math` | 1061 | 1975 | 1.86× |

**Takeaways**

- **Love wave**: C++ is ~15–25% faster. The SH propagator (`dltar1`) uses a
  simple 2-element vector, leaving little headroom; both compilers do well.
- **Rayleigh wave**: C++ is ~60–90% faster. The P-SV propagator (`dltar4`)
  involves a 5×5 Dunkin compound matrix recomputed every layer; the C++ static
  functions are fully inlined by the compiler, eliminating call overhead and
  enabling more aggressive register allocation compared to the Fortran version
  which is constrained by `COMMON` block aliasing rules.
- **O3 vs O2**: negligible difference (±3%). The algorithm is
  branch-heavy scalar code with irregular control flow; auto-vectorisation
  does not apply, and extra loop unrolling from O3 adds marginal overhead.
- **`-march=native -ffast-math`**: benefits Fortran Love wave slightly, but
  hurts Fortran Rayleigh wave — `ffast-math` reorders floating-point
  operations in the Dunkin matrix, altering the numerical path.
- **Recommendation**: use `-O2`. It gives the same performance as `-O3` with
  shorter compile times and no risk from aggressive floating-point rewrites.

---

## Plot

```bash
cd build
python3 ../scripts/plot_comparison.py
```

Generates `comparison.png` — six panels (Love/Rayleigh × fundamental/overtone
× phase/group), each showing the C++ and Fortran dispersion curves overlaid,
with a relative-difference residual panel below.

![comparison](build/comparison.png)

---

## Fortran-to-C++ translation notes

| Topic | Approach |
|---|---|
| Array indexing | Fortran 1-based → C++ 0-based throughout |
| `COMMON` blocks (`getsol_module`, `sphere_module`) | Converted to reference parameters (`del1st`, `dhalf`); no global state |
| `goto`-based control flow | Converted to `while` + `break`/`continue` |
| `implicit double precision (a-h,o-z)` | All matching variables explicitly `double` |
| Dunkin 5×5 matrix `ca(j,i)` (Fortran column-major) | `ca[j][i]` (C++ row-major); both `dnka` and `dltar4` use the same logical indexing, so the product `ee[i] = Σ_j e[j]·ca[j][i]` is preserved |
| `sngl()` precision reduction in group-velocity formula | Replicated as `(float)` cast to match Fortran output exactly |
| Neville iteration index arithmetic | 1-based Fortran `x(j)`, `y(m+1)` mapped to 0-based `x[j-1]`, `y[m]` |

---

## References

- Herrmann R. B. (1987). *Computer Programs in Seismology*, Vol. IV, Saint Louis University.

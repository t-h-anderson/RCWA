# RCWA — OOP refactor

This directory holds an object-oriented restructuring of the original
`RCWA/*.m` scripts. The originals are left in place untouched so any
thesis run can still be reproduced; the OOP code lives entirely under
`+rcwa/`.

The refactor was driven by the request **"make it an OOP codebase,
write tests, and find bugs"**. Findings are catalogued in
[`BUGS.md`](BUGS.md). Per the agreed scope, the OOP port preserves the
originals' numerical behaviour so thesis comparisons stay valid —
bugs are documented and surfaced via regression tests, not silently
fixed.

## Layout

```
oop/
├── +rcwa/                       MATLAB package: rcwa.*
│   ├── Constants.m              physical constants (mu0, eps0, c, h, q, eta0)
│   ├── Units.m                  m↔nm, deg↔rad, eV↔nm helpers
│   ├── Config.m                 replaces the `loc` struct + DefaultLoc/MishaLoc/VerifyLoc/ApplyVarargin
│   ├── Geometry.m               replaces BuildMaterial (layer categorisation)
│   ├── BandgapProfile.m         replaces EgProfile (with `legacy` toggle for BUG-001)
│   ├── PermittivityBuilder.m    replaces BuildEps + BuildEpsF
│   ├── InverseFourier.m         replaces InverseFourier.m
│   ├── Solver.m                 replaces RCWA.m (the actual algorithm)
│   ├── SpectrumRunner.m         replaces RCWAsetup.m (the orchestration)
│   ├── Specular.m               replaces specular.m
│   ├── Spectrum.m               replaces W2inmim2AM15G.m
│   ├── +gratings/               profile shapes (Square, Sinusoidal, Pyramid, Spherical)
│   └── +materials/              dispersive permittivities (Glass, Silver, AZO, FTO, ASiH, Constant)
├── tests/                       MATLAB unittest suite
│   ├── run_all.m                runner that adds the right paths
│   └── *Test.m                  per-class tests
├── BUGS.md                      defects found while porting
└── README.md
```

## Running the tests

```matlab
cd /path/to/RCWA/oop/tests
run_all
```

`run_all.m` adds both the original RCWA directory (for the legacy
data tables and cross-check calls) and the OOP package directory to
the path before invoking `TestSuite.fromFolder(pwd)`.

## What's covered, what isn't

The tests cross-check the OOP classes against the originals where
practical (geometry, Fourier of the grating slice, EgProfile,
InverseFourier, the full RCWA solver on a tiny homogeneous problem).
There is also a smoke run of `SpectrumRunner` end-to-end on the
smallest config the original ever accepted.

The tests were written and read-verified against the MATLAB unittest
API but were **not executed** in this environment (no MATLAB
available). The first run on a real MATLAB install should be treated
as the initial validation pass.

## Quick start

```matlab
% Reproduce the default thesis-style sweep
cfg     = rcwa.Config();              % uses DefaultLoc-style defaults
runner  = rcwa.SpectrumRunner(cfg);
out     = runner.run();
out.mAicm2JOpt    % optical short-circuit current density (mA/cm^2)

% Override individual parameters
cfg = rcwa.Config('Nt', 7, 'pol', 0, 'nlambda', 30);

% Use a struct that came from an existing run
cfg = rcwa.Config(existing_loc_struct);
```

## Materials

Wavelength-dependent materials are first-class objects:

```matlab
m = rcwa.materials.Glass();
m.permittivity([400 500 600])     % column of eps_r at three wavelengths

m = rcwa.materials.ASiH('eVEg', [1.6; 1.7; 1.8]);
m.permittivity(550)               % column over Eg, one wavelength
```

The original behaviour is preserved (see BUGS.md, BUG-007 in
particular) but the `squared` flag on `rcwa.materials.AZO` lets a user
opt into the physically correct `(n+ik)^2` form.

## Gratings

```matlab
g = rcwa.gratings.Square('nmLx', 400, 'zeta', 0.5, 'nmLg', 50);
g.relief(linspace(-200, 200, 51))   % z = g(x)

% Or pick from a Config:
g = rcwa.gratings.Grating.fromConfig(cfg);
```

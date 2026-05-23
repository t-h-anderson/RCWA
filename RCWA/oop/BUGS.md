# Bugs found in the original RCWA codebase

This document lists defects discovered while porting the original
MATLAB scripts to the `+rcwa` OOP package. Each entry records the
file/line, the symptom, an analysis, and any impact on results.

> Per the agreed scope, **bugs are documented, not fixed in place.** The
> OOP port faithfully reproduces the originals' behaviour (with notes
> in the code) so thesis numbers stay comparable. Where the original is
> demonstrably crash-y rather than silently wrong (BUG-008, BUG-009),
> the OOP port chooses the obviously-intended behaviour and exposes a
> regression test.

---

## BUG-001 — `EgProfile.m:32`: i-layer perturbation offset is `nmLn`, should be `nmLp`

```matlab
Eg(ilayer == 1) = loc.Eg0 + loc.A .* ((0.5*(sin( ...
    2*pi * (nmz(ilayer == 1) - loc.nmLn) ...      % <-- should be loc.nmLp
    * loc.kappa/loc.nmLi + 2*pi*loc.phi )+1)).^loc.alpha);
```

**Why it's wrong.** The i-layer starts at `z = nmLp`, not `z = nmLn`. The
intent is "phase 0 at the start of the i-layer". The literal effect of
the bug is to shift the periodic perturbation's phase by
`2*pi*kappa*(nmLp - nmLn)/nmLi`.

**Impact.** Zero when `nmLp == nmLn` (which is the case for the default
config: both 15 nm). Non-zero, and silent, whenever they differ.

**Reproduced by.** `BandgapProfileTest.bugDocumentedRegression`.

The OOP `BandgapProfile` exposes a `legacy` flag (default `true`) that
reproduces the original behaviour for thesis comparability and a
`legacy = false` mode that uses `nmLp`.

---

## BUG-002 — `BuildEps.m:45`: air permittivity assigned as `nsa`, not `nsa^2`

```matlab
eps(:, mat_cat == 4) = loc.nsa;        % loc.nsa is the refractive index
```

**Why it's wrong.** Every other branch in `BuildEps.m` is assigning a
permittivity (`epsd`, `epsW`, `epsm`, `epsJ`). The air region should
follow suit: `eps_air = nsa^2`. With the default `nsa = 1 + 1e-6 i`
this barely matters (1 vs 1 + 2e-6 i), but it would matter if the user
deliberately used a higher-index ambient.

**Impact.** Negligible at default (`nsa ≈ 1`), but the bug means
results scale wrong as a function of the ambient. The OOP
`PermittivityBuilder` preserves the behaviour with a documentation
comment.

---

## BUG-003 — `BuildEps.m:21` vs `:47`: case mismatch on window permittivity field

```matlab
% line 21 -- writes lowercase:
loc.epsw = glasseps;
% line 47 -- reads capital:
eps(:, mat_cat == 3) = loc.epsW;
```

**Why it's wrong.** MATLAB struct fieldnames are case-sensitive. With
`wmat = 1` (glass window), `BuildEps` writes a brand-new field
`loc.epsw`, leaving `loc.epsW` untouched at whatever DefaultLoc set
(`4 + 1e-6i`). The dispersive glass refractive index data is silently
ignored for the window.

**Impact.** Whenever `wmat == 1` was used, the window had a constant
permittivity `4 + 1e-6 i` instead of glass(λ)². Likely affected any
thesis run with `wmat = 1`.

**Reproduced by.** Code-level review only — the OOP port reads
`cfgOut.epsW` consistently after assigning the resolved value.

---

## BUG-004 — `BuildEps.m:19`: window-material block keyed on `dmat` instead of `wmat`

```matlab
if(loc.wmat==1)
    ...
elseif(loc.dmat==2)         % <-- should be loc.wmat==2
    glasseps = AZO(nmlambda);
    loc.epsw = glasseps;
end
```

**Why it's wrong.** Inside the window block, the `elseif` keys on
`dmat` (the grating dielectric material), not `wmat`. Setting
`wmat = 2` (AZO window) silently does nothing if `dmat != 2`, and
conversely setting `dmat = 2` for the grating dielectric also sets the
window dielectric.

**Impact.** Window material assignment is unreliable. Combined with
BUG-003 (case mismatch) this branch never had any effect even when the
condition was met. The OOP port keys correctly on `wmat`.

---

## BUG-005 — `RCWA.m:79`: `diag(2*loc.Nt+1)` returns a scalar, not `eye(N)`

```matlab
epst_temp = epsz(2*loc.Nt+1) * diag(2*loc.Nt+1);
```

**Why it's wrong.** `diag(s)` with a scalar argument `s` simply wraps
it (returns the 1×1 matrix `s`). The author plainly meant
`eye(2*loc.Nt+1)`. The resulting `epst_temp` is a single scalar
`(2*Nt+1) * epsz(N)`, which then gets assigned to a (N,N) slice of
`epst` via broadcasting (filling every entry with the scalar).

**Impact.** None on final results: `epst` is stored but never read
after the homogeneous-branch block — only `invepst` (computed
correctly as `eye(N)/epsz(N)` on line 80) is used downstream. So the
bug just means a 3D buffer gets filled with a wrong but unused value.
Verified by searching for downstream reads of `epst` in `RCWA.m`.

---

## BUG-007 — `Materials/AZO.m`: returns refractive index `n+ik`, callers treat as permittivity

```matlab
% AZO.m -- documented as "refractive index for wavelength lam0"
val = interp1(nmlambdas, AZO_n + 1i.*AZO_k, nmlambda);
```

```matlab
% BuildEps.m (both branches):
glasseps = AZO(nmlambda);   % <-- this is (n+ik), not eps
loc.epsd = glasseps;
```

**Why it's wrong.** The other branches of `BuildEps` square refractive
indices before storing (e.g. silver: `(Silver(...)^2)`, glass:
`(glass(...) + 1e-3 i).^2`). The AZO branch does not. The values
stored as `epsd`/`epsw` are off by a factor of approximately n,
typically halving the real part and flipping the imaginary-part
ordering.

**Impact.** Whenever `dmat = 2` or `wmat = 2` (AZO), the dielectric
permittivity is wrong. This is the default `dmat = 2` in `DefaultLoc`.

The OOP `rcwa.materials.AZO` preserves the original behaviour by
default (returns `n+ik`) and exposes a `squared` flag for the corrected
form.

---

## BUG-008 — `VerifyLoc.m:14-16`: writes `Lp`/`Li`/`Ln`, the actual fields are `nmLp`/`nmLi`/`nmLn`

```matlab
if(loc.Nz == 0)
    loc.Lp = 0;
    loc.Li = 0;
    loc.Ln = 0;
end
```

**Why it's wrong.** The DefaultLoc fields are `nmLp`, `nmLi`, `nmLn`.
The above lines just create new fields and leave the real ones
unchanged. The "no junction" path in `VerifyLoc` is a no-op.

**Impact.** Configurations with `Nz = 0` still carry the junction
thicknesses and end up with junction slices having zero count but the
overall layout still subtracting nmLp+nmLi+nmLn. The OOP `Config.verify`
zeroes the correct fields.

**Reproduced by.** `ConfigTest.verifyZeroesOutCorrectFields`.

---

## BUG-009 — `RCWAsetup.m:109`: divide by zero when `loc.Nz == 0`

```matlab
nmdJ = (loc.nmLp+loc.nmLi+loc.nmLn)/loc.Nz;
```

**Why it's wrong.** Combined with BUG-008 (which doesn't actually zero
the layer thicknesses), a user setting `Nz = 0` ends up with
`nmdJ = Inf`, silently propagated into the absorption integrals.

**Impact.** Any `Nz = 0` run gives `mAicm2JOpt = NaN` (Inf · 0 = NaN).
OOP `SpectrumRunner.junctionSliceThickness` returns `0` in that case.

**Reproduced by.** `SpectrumRunnerTest.junctionSliceThicknessGuardsZero`.

---

## BUG-010 — `RCWAsetup.m:26`: `fliplr` of `EgProfile` output reverses the i-layer profile

```matlab
Eg(mat_cat == 2) = fliplr(EgProfile(nmz(mat_cat == 2)-nmgtop, loc));
```

**Why it's wrong.** `nmz(mat_cat == 2)` is already in top-down order
(`BuildMaterial` calls `fliplr` on `nmz` and `nmdz` near the end).
`EgProfile` evaluated on these top-down z values returns Eg values
in the same top-down ordering — naturally giving the n-layer first
and the p-layer last (correct: n faces the air side).

Applying `fliplr` then inverts the entire profile, placing the
p-layer's Eg at the top (air side) and the n-layer's at the bottom
(mirror side). Because p- and n-layers carry the same constant
`Eg = 1.95`, the only visible effect is that the i-layer's perturbation
is traversed in z-reverse — which manifests as an effective phase shift
in any `kappa > 0` periodic Eg modulation.

**Impact.** Cosmetic (re-labels left/right) for uniform i-layer; for a
modulated i-layer it shifts the perturbation phase relative to the
physical p-i-n stack.

The OOP `SpectrumRunner.junctionBandgap` preserves the flip and
documents it; correcting it would require coordinating with BUG-001.

---

## BUG-011 — `RCWAsetup.m:141-145`: accumulators allocated inside the `tt` loop

```matlab
for ll = 1 : loc.nlambda
    Q_temp = zeros(1, loc.Nx, Nz);
    for tt = 1 : loc.ntheta
        ...
        Et          = zeros(loc.ntheta, 3, loc.Nx, Nz);   % reset each tt!
        normE2t     = zeros(loc.ntheta,    loc.Nx, Nz);   % reset each tt!
        die_abst    = zeros(loc.ntheta, 1);               % reset each tt!
        metal_abst  = zeros(loc.ntheta, 1);               % reset each tt!
        ...
        Et(tt, :, :, :) = E_temp;
        ...
    end
    El(ll, :, :, :, :)    = Et;             % only the last-tt slice non-zero
    normE2l(ll, :, :, :)  = normE2t;        % same
    die_absl(ll, :)       = die_abst;       % same
    metal_absl(ll, :)     = metal_abst;     % same
end
```

**Why it's wrong.** `Et`, `normE2t`, `die_abst`, `metal_abst` are
allocated fresh inside the `tt` loop. Only their `tt`-th slot is
populated each iteration; the previous iterations' values are wiped
out. When the loop ends, only the very last theta's data is in
slot `ntheta`; slots `1..ntheta-1` are zero.

**Impact.** Zero when `ntheta == 1` (the default — all single-angle
runs are fine). For any multi-angle run (`degtheta0 ≠ degtheta1`,
`ntheta > 1`) the per-angle outputs `El`, `normE2l`, `die_absl`,
`metal_absl` are corrupted: only the last angle survives.
**The summed Q and JOpt are still correct**, because the `Q_temp`
accumulator (allocated correctly outside the `tt` loop) does
correctly sum across angles.

The OOP `SpectrumRunner.absorption` allocates each accumulator at the
right scope.

---

## BUG-012 — `Faryad.m`: leftover scratch with side-effecting calls and dead allocation

```matlab
%function eps = Faryad(nmz, nmlambda)        % <-- signature commented out
aSiHGC('eVEg', 1.95, 'nmlambda', 400, 'optimum',1)
temp = xlsread('epsaSiC195.xlsx');
temp(1,1)
temp(1,2)
temp(1,3)
...
eps = zeros(length(nmz),1);
```

**Why it's wrong.** Without the `function` line uncommented, this file
is a script: calling `Faryad(nmz, nmlambda)` would not behave as a
function call. Inside, the call signature `aSiHGC('eVEg', ...)` does
not match the `aSiHGC` actually defined in `Materials/aSiH/aSiHGC.m`
(which takes positional args, not name/value pairs), so it errors.
The `temp(1,1)` etc. are dead expressions whose results are dropped.
Finally `eps = zeros(length(nmz),1)` is a return value that never
makes it back as a function output.

**Impact.** Any path with `loc.material == 2` (the Faryad SPIE test
case branch in `BuildEps`) cannot run. Likely never exercised since
the thesis development moved on.

---

## BUG-013 — `gepsz.m:10`: dead expression with no assignment

```matlab
kron(ones(length(nmx),1), nmz);
```

**Why it's wrong.** Side-effect free MATLAB expression at statement
level (no `;` ostensibly, but MATLAB still discards). Pure dead code.

**Impact.** None on results, but wastes time allocating a `len(nmx) x
len(nmz)` matrix only to drop it.

---

## BUG-014 — `RCWAsetup.m:16-17`: `nmx` immediately overwritten

```matlab
nmx = linspace(-0.25*loc.nmLx,  0.75*loc.nmLx, loc.Nx);
nmx = linspace(-0.5 *loc.nmLx,  0.5 *loc.nmLx, loc.Nx);
```

**Why it's wrong.** The first line is a no-op (overwritten by the
second). Dead code; left over from a manual experiment.

**Impact.** None.

---

## BUG-015 — `ApplyVarargin.m:18-35`: unknown fields in the leading struct are silently dropped

```matlab
if round(nArgs/2) ~= nArgs/2
    inStruct = varargin{1};
    if isstruct(inStruct)
        ...
        for i = 1:length(inNames)
            if any(strcmp(inNames{i}, optionNames))
                loc.(inNames{i}) = inStruct.(inNames{i});
            end
        end
    else
        error(...)
    end
end
```

**Why it's wrong.** Mismatched field names in the struct are silently
ignored, while the later name/value-pair loop *does* error on
unknown names. So `ApplyVarargin(loc, 'foo', 1)` errors, but
`ApplyVarargin(loc, struct('foo', 1))` silently does nothing.

**Impact.** Easy to typo a struct field and never notice.
OOP `rcwa.Config` errors on unknown name/value pairs and only silently
ignores struct fields that don't exist on the class (still a small
asymmetry, but rationalised: struct→Config is "copy what you can"
while name/value implies "set this property explicitly").

---

## Minor / style observations (not bugs)

- `RCWAold.m` is a pre-stable-algorithm draft of `RCWA.m`. Keep one,
  delete the other.
- `Workspace.m` / `Workspace2.m` are scratch scripts referencing
  outer-scope variables (`loc`, `nmlambda`, `El`, `mat_cat`,
  `normE2l`) and a 10 000 × 10 000 timing experiment. They're useful
  as REPL snippets, not as committed code.
- `DefaultLoc` in `MishaLoc.m` shadows the real `DefaultLoc.m`
  function name — `MishaLoc.m`'s first line is
  `function [loc] = DefaultLoc()`. If a user accidentally calls
  `MishaLoc` expecting Misha's preset, MATLAB might dispatch the
  wrong function depending on path order. The OOP `Config.withMisha()`
  is a proper alternative.
- `Materials/FTOnk.m` references a `.mat` file (`FTOnk1200.mat`) that
  doesn't exist in the repo — only `FTOnk.mat` does. Loading would
  fail.
- `printloc.m`'s `catch` block calls `sprintf` and discards the
  result — should be `disp(sprintf(...))` or `warning(...)`.
- `epsn.m:79`'s commented-out `epsz(2:2:end,:) = conj(epsz(2:2:end,:));`
  is an unresolved sign-convention question between explicit and FFT
  Fourier methods.

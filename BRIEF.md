# BRIEF — FallingPlate2D.jl

Handover for a fresh session. Everything needed to build the code is here.
Written 11 September 2026.

---

## 1. Objective

A **fast-running 2D code for toy problems**: a freely falling plate in an
incompressible quiescent fluid, Reynolds number 10 to 1000, extended to
**permeable plates** via a Darcy term. Runs on the user's current laptop.
Published on GitHub.

This is an **exploration tool**, not a production solver. Optimise for speed of
iteration and ease of changing geometry, not for quantitative accuracy.

---

## 2. Decision already taken — do not relitigate

**Build a thin package on top of [WaterLily.jl](https://github.com/WaterLily-jl/WaterLily.jl).**
Do not write a Navier–Stokes solver from scratch.

Rejected, with reasons:

| Option | Why not |
|---|---|
| Bespoke conformal vorticity–streamfunction solver | 2–3× the cost, high risk, and the conformal map locks geometry to ellipses — wrong for toy problems |
| Basilisk | Faster per case via adaptivity, but macOS build friction and sandbox-dependent moving-body code |
| OpenLB | GPU-oriented and heavy; no benefit without a GPU |
| The user's existing OpenFOAM setup | Not fast, not toy-friendly. **Keep it as the quantitative reference.** |
| Anything GPU | Apple GPUs have no FP64, no CUDA, OpenGL capped at 4.1. Run Float64 on the CPU. |

---

## 3. Critical context — read before starting

**The user already has a validated 2D falling-plate solver.** OpenFOAM v2306,
overset mesh, `overPimpleDyMFoam`, modified `sixDoFRigidBodyMotion`,
non-inertial frame, partitioned weak coupling. Verified with a √2-refinement
study, validated against Andersen et al. (2005) at Re = 1147. Published as
Pandi, Gungor, Bose, Attili & Viola, *Effect of Transverse Gust on Free-Falling
Plates*, [arXiv:2304.04317](https://arxiv.org/abs/2304.04317).

**Do not propose rebuilding it.** This new code exists only to be *fast for
exploration*. OpenFOAM stays the reference for any quantitative claim.

Related group work: permeable plates (Bose et al. 2025), permeable disks
(Cummins et al. 2017). Permeability is the research direction that matters.

---

## 4. Hardware and toolchain

- **Apple M2, 4 performance + 4 efficiency cores, 24 GB.** Use `-t 4`
  (performance cores only). Do not assume 8.
- **Julia 1.13.0** (aarch64, Apple Darwin), installed via `juliaup` and set as
  the default `release` channel. Confirmed working. Launch with
  `julia -t 4 --project=.`
- `gh` authenticated as `ignaziomviola`; git and Homebrew present.
- **Repository at `~/Developer/FallingPlate2D.jl`**, git on `main`.
  Deliberately **outside OneDrive** — git indexes and OneDrive sync corrupt
  each other, and simulation output must never sync. `logs/` is ignored.
- Package environment: `Project.toml` has WaterLily v1.8.0 and StaticArrays.
  A second environment `examples/` adds ParametricBodies and Plots for the
  upstream example only; keep plotting out of the package dependencies.
- **Julia buffers stdout when redirected to a file.** Every long script must
  `flush(stdout)` after each progress line or the log stays empty until exit.
  macOS has no `timeout` command.

---

## 5. Physics and parameters

Nondimensionalise on chord `c`, descent speed `v`, fluid density.

| Symbol | Definition |
|---|---|
| `Re` | `c·v/ν` |
| `χ` | `c/f`, chord-to-thickness ratio (Andersen's `β = 1/χ`) |
| `I*` | `32·I/(π·ρ_f·c⁴)` |
| `ρ*` | `ρ_s/ρ_f` |
| `Da` | `K/c²`, Darcy number |

Descent speed estimate: `v ≈ sqrt(2·(ρ* − 1)·f·g)`.

**Regimes:** steady descent below `Re ≈ 100`; fluttering-to-tumbling transition
with increasing `I*` above it, transition `I*` around 0.2–0.4.

**Duration:** one tumbling revolution ≈ 8 convective times `c/v`. A useful run
is **100–200 `c/v`**. This is the dominant cost driver.

---

## 6. Validation target

Andersen, Pesavento & Wang, *JFM* **541** (2005) 65–90, table 1. Aluminium
plates in water, `ν = 0.0089 cm²/s`:

| `f` (cm) | `β = 1/χ` | `I*` | `Re` | `⟨vx⟩` (cm/s) | `⟨vy⟩` (cm/s) | `⟨ω⟩` (rad/s) |
|---|---|---|---|---|---|---|
| 0.162 | 1/5 | 0.48 | 737 | 14.2 | −8.1 | 5.5 |
| 0.162 | 1/6 | 0.39 | 1770 | 15.2 | −16.2 | 16.1 |
| 0.081 | 1/8 | 0.29 | 837 | 15.9 | −11.5 | 14.5 |
| 0.081 | 1/14 | 0.16 | 1147 | 19.8 | −9.0 | 6.8 |

**Use the `β = 1/14`, `Re = 1147` row as the primary target** — the user's own
OpenFOAM paper validated against it, so there is a published comparison curve.

Note Andersen et al. used an **elliptical** cross-section to avoid the corner
singularity; their velocity scale and `I*` definitions above assume that.

---

## 7. Resolution and expected runtime

WaterLily uses the boundary data immersion method, which smears the interface
over ca. 2 cells. **Need ca. 8 cells across the thickness**, so `Δ ≈ c/(8χ)`.

| Mode | Setup | Grid | Measured (Phase A, 4 threads, Float64) |
|---|---|---|---|
| **Toy** | `χ = 4`, `Δ = c/32`, domain 8c × 16c | 258 × 514 | 42–50 ms/step, `Δt ≈ 0.35` cells, ca. 92 steps per `c/v` → **6.5–8 min per 100 `c/v`** |
| **Check** | `χ = 10`, `Δ = c/80`, domain 8c × 16c | 642 × 1282 | 100 ms/step, `Δt ≈ 0.2–0.3` cells, ca. 320 steps per `c/v` → **ca. 55 min per 100 `c/v`** |

Measured on 11 September 2026 with `scripts/phaseA_timing.jl` for a freely
rotating ellipse with re-measurement every step; free translation adds two
more moment/force evaluations per step, so expect ca. 10 % more. Cost per
cell-step is 125 ns on the check grid and 375 ns on the toy grid (thread
efficiency is poor on the small grid). Poisson iterations were 1–3 per step
in stable runs. Thin plates are where BDIM costs you; a thicker section is
what buys the toy-mode turnaround.

Use `T = Float64`.

---

## 8. WaterLily API — verified against v1.8.0 (Phase A, 11 September 2026)

- `RigidMap(x₀, θ; xₚ, V, ω)` exists and is exported. Fields: centre `x₀`,
  angle `θ`, pivot offset `xₚ`, linear velocity `V`, angular velocity `ω`.
  Update with `sim.body = setmap(sim.body; θ, ω, V, x₀)`; then
  `sim_step!(sim; remeasure=true)`. Use `AutoBody(sdf, RigidMap(...))` with a
  signed-distance function in body coordinates (the `RigidMap` docstring
  pattern); no need for ParametricBodies.
- Forces and moments: `pressure_force`, `viscous_force`, `total_force`,
  `pressure_moment(x₀, sim)`, `viscous_moment`, `total_moment` — all present.
  The upstream ellipse example uses the pressure moment only; **use
  `total_moment` in Phase B**, viscous damping matters for a thin plate at
  `Re ≤ 1000`.
- Body force keyword is **`g(i, x, t)`** — three arguments, not two.
- Far-field velocity `uBC` may be a function `uBC(i, x, t)`; WaterLily then
  adds `dU/dt` to the momentum equation automatically by forward-mode
  differentiation of that function. A closure reading the body velocity from
  mutable state differentiates to zero, so **supply the frame acceleration
  explicitly through `g`** (the rigid-body integrator knows it exactly) and
  feed the far-field velocity through a `Ref`.
- `sim_step!(sim; udf=f, kwargs...)` calls `f` on the right-hand side in both
  predictor and corrector: the hook for an explicit Brinkman term (Phase C).
- The implicit-added-mass update `α = (M + Iₐα)/(Iₛ + Iₐ)` in the upstream
  `TwoD_FreeRotatingEllipse.jl` is the scheme Phase B proposes. Run
  unmodified on this laptop (`examples/run_upstream_ellipse.jl`, Float32,
  384 × 256): both cases complete, the chaotic one (`χ = 2`, `ρ* = 1`,
  pivot −0.16, 50 `c/v`) with `|ω| c/v ≤ 1.2`, the stable one (`χ ≈ 8`,
  pivot 0.25, `Re = 250`) inverting to 165° and stopping; 27 min in total
  including package compilation and 650 plot frames. Log:
  `logs/phaseA_upstream.log`.

### Two WaterLily pitfalls found in Phase B (both fixed in `src/simulation.jl`)

1. **Convective exit under a coupled body acceleration.** WaterLily's
   `exitBC=true` recomputes and flux-corrects the exit only in the predictor;
   the corrector and projection keep that value while the uniform body force
   `g` has changed the interior flux, so the pressure problem is slightly
   incompatible every step. With a prescribed smooth acceleration this is
   harmless; with the body acceleration fed back from the plate force the
   error grows and the run dies at the exit face within ca. 4 `c/U` (any
   density ratio, with or without periodic sides, faster with smaller `Δt`).
   Fix: `update_exit!` recomputes the first-order upwind exit profile once per
   step from the start-of-step velocity, flux-corrects it to the new far-field
   speed, and serves it through the `uBC` function, so every `BC!` call in the
   step sees the same consistent values. Use `exitBC=false` in `Simulation`.
2. **Interior body velocity.** `measure!` fills the imposed velocity `V` only
   within `2 + ϵ` cells of the surface and leaves it zero deeper inside. For a
   moving (translating or rotating) body thicker than ca. 6 cells the interior
   velocity is discontinuous and has divergence where the Poisson operator has
   no coefficients; the unremovable residual gives an interior pressure of
   O(10⁴) that leaks into the surface band. Consequence: force *along* the
   motion (3× too large) and a pressure moment that *drives* the rotation.
   **This, not a coupling instability, was the `ω c/U ≈ 5` blow-up seen in
   Phase A**, and it explains why the maintainers' thin-body examples work.
   Fix: `rigid_velocity_fill!` sets `V = v + ω × r` on every face with
   `μ₀ == 0` after `measure!`; interior pressure becomes smooth, drag opposes
   the motion, and the pressure moment damps the rotation.

Verified after the fixes (`scripts/phaseB_fillcheck.jl`): first-step
acceleration equals `(m − A)g/(m + m₂₂)` with Lamb's `m₂₂` to 1 %; plate
translating at `0.5 U` has `C_D = 2.5` against the motion; constant rotation
at `ω c/U = 0.5` gives `M c⁻²U⁻² = −0.009`; the 3-DOF fall of a `χ = 4`,
`ρ* = 2.7`, `Re = 500` ellipse flutters stably for 12 `c/U` at `<vY> = −0.64 U`.

Translational coupling: since the far fluid accelerates with the frame there
is no pressure gradient at infinity, so the update is exactly
`(m + Mₐ)·a_new = F + Mₐ·a_old + (m − A)·g` with `Mₐ` the lab-frame added-mass
tensor — no displaced-volume term beyond buoyancy. Convention: forces and
moments on the body are `−total_force(sim)` and `−total_moment(x₀, sim)[1]`.

Guards kept in `run!`: abort on non-finite state, `|ω| c/U > 10`, or
`Δt < 10⁻³` cells, printing where `|u|` is largest (exit face vs body).

---

## 9. Build plan

### Phase A — verify the API — **DONE 11 September 2026**
Environment created, API verified (section 8), runtimes measured (section
7), coupling stability probed (section 8). Outcome: the plan stands.
`scripts/phaseA_timing_defs.jl` already contains a Float64 ellipse SDF,
`make_sim`, and the rotation loop — the seed for `src/bodies.jl` and
`src/rigid.jl`.

### Phase B — 3-DOF free fall — **code done 12 September 2026; validation run below**
- Body signed distance function for an ellipse/plate parameterised by `χ`.
- **Translating frame:** keep the body near the domain centre; apply the
  negative frame acceleration as a uniform body force through `g(i,t)`; set the
  far-field inflow to minus the body velocity. This is what makes long
  trajectories affordable — the plate falls O(100) chords, and a lab-frame
  domain would be absurd.
- Rigid-body integration: solve `(M + M_a)·a = F_fluid + (M − ρ_f·V)·g`.
  Use `total_force`/`total_moment`, not pressure only. Add the abort guards
  from section 8.
- **Treat added mass implicitly.** The added-mass tensor for an ellipse with
  semi-axes `a, b` is analytic (Lamb): `m₁₁ = ρ_f·π·b²`, `m₂₂ = ρ_f·π·a²`,
  `I_a = (π/8)·ρ_f·(a² − b²)²`. Moving it to the left-hand side removes the
  sub-iteration that cost Andersen et al. a factor of 10–20, and removes the
  instability at low density ratio.
- Validate against the `β = 1/14`, `Re = 1147` row.

**Risks:** coupling stability at low `ρ*`; BDIM force noise at coarse
resolution; torque is a small residual of nearly cancelling surface integrals,
so it is the sensitive quantity.

### Phase C — permeability (~half a day)
Add a Brinkman resistance `−(ν/K)·χ_s·(u − u_s)` to the BDIM forcing. BDIM is
already a smoothed volume-penalisation method, so this is a modification of
existing forcing rather than a new capability — the single strongest reason
WaterLily was chosen.

- Treat the term **implicitly**: it is stiff (`ν/K` large) but local and linear.
- **Resolution floor:** must resolve the Brinkman screening length `c·√Da`. At
  `Da = 10⁻⁴` this needs `Δ ≲ 10⁻²·c`, coarser than the thickness already
  demands, so it is free. It only bites below ca. `Da = 10⁻⁶`, where cost scales
  as `Da^(−1/2)` per direction.

### Phase D — repository hygiene (~2–3 h)
Tests, GitHub Actions running the fast ones, README, examples.

**Total estimate: 1.5–2 days, 0.9–1.6 M tokens.**

---

## 10. Repository layout (as built)

```
FallingPlate2D.jl
├── BRIEF.md              # this file
├── Project.toml          # WaterLily 1.8, StaticArrays; Test in [extras]
├── src/FallingPlate2D.jl # module; frames and conventions in the docstring
├── src/bodies.jl         # Ellipse, RoundedPlate: SDF, area, inertia, added mass, I*
├── src/rigid.jl          # RigidState, PlateDynamics, implicit-added-mass step!
├── src/simulation.jl     # falling_plate, step!, run!, update_exit!, rigid_velocity_fill!
├── src/diagnostics.jl    # readrun, summarise (10-line comparison with Andersen)
├── src/darcy.jl          # Phase C (to do)
├── examples/free_fall.jl # driver: julia -t 4 --project=. examples/free_fall.jl χ ρ Re L tend name
├── examples/             # own Project.toml (Plots, ParametricBodies) for the upstream example
├── scripts/phaseA_*.jl, phaseB_*.jl   # diagnostics that produced the findings above
├── test/runtests.jl      # fast checks (< 1 min): properties, signs, interior pressure, short fall
├── out/  logs/           # ignored by git; simulation CSVs and session logs
```

Run a case: `julia -t 4 --project=. examples/free_fall.jl 14 2.7 2090 80 60 name`
(χ, ρ*, Re_g, L cells, tend in c/U, output name). Tests: `julia --project=. -e 'using Pkg; Pkg.test()'`.

---

## 11. Working practices

- **One phase per session.** Token cost is dominated by re-reading context, not
  by writing code.
- Keep modules to 200–400 lines so only the edited file is re-read.
- **Terse numeric validation output** — a script printing ten lines of errors
  and observed orders is far cheaper to read than a field dump.
- Run simulations in the background, logs to file, summaries only into context.
- Update this BRIEF at the end of each phase so the next session resumes cheaply.

---

## 12. Do not

- Rebuild the OpenFOAM work.
- Use the GPU, or Float32.
- Put simulation output anywhere inside OneDrive.
- Aim for publication-grade force decomposition — that is OpenFOAM's job.

---

## 13. References

- Andersen, Pesavento & Wang, *JFM* **541** (2005) 65–90 —
  <https://dragonfly.tam.cornell.edu/publications/S002211200500594Xa.pdf>
- Pandi, Gungor, Bose, Attili & Viola — <https://arxiv.org/abs/2304.04317>
- WaterLily.jl — <https://github.com/WaterLily-jl/WaterLily.jl>,
  paper <https://arxiv.org/abs/2407.16032>
- WaterLily-Examples — <https://github.com/WaterLily-jl/WaterLily-Examples>

Session memory lives in
`~/.claude/projects/.../CFD-free-falling-2D-plate/memory/`.

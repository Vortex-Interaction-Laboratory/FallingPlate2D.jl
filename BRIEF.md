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

Andersen et al. define `Re = c⟨vy⟩/ν` on the **mean descent speed** and
`I* = 32 I/(π ρ_f c⁴)` with the **rectangular** inertia `I = ρ_s c f (c² + f²)/12`,
so `I* = (8/3π) ρ* β (1 + β²)`; that reproduces table 1. For an ellipse the
same definition gives `I* = ρ* β (1 + β²)/2`, i.e. 59 % of the rectangle's
value at equal `ρ*` and `β`. Their DNS therefore used a **thicker ellipse,
`β = 1/8`, with the density ratio chosen to match `I*`**: `I* = 0.16` needs
`ρ* = 2.52`. That is the validation configuration for this code:

    julia -t 4 --project=. examples/free_fall.jl 8 2.52 2090 64 60 name

(`Re_g = 2090` on the gravitational speed gives `Re ≈ 1147` on `⟨vy⟩` if
`⟨vy⟩ ≈ 0.55 U` as in the experiment; report the realised value.)

**Thin-plate result (Phase B, 12 September 2026):** the literal plate,
`χ = 14`, `ρ* = 2.7` (`I* = 0.097` as an ellipse), at `L = 80` glides,
rises briefly at the end of the glide, drops edge-on at `φ ≈ −1.3` and
accelerates almost without drag to `6 U` over 15 chords while turning only
slowly (`ω c/U` 0.1–0.5); the energy budget of a frictionless edge-on drop
allows `9 U`, so this is coherent low-`I*` 2D dynamics, but it drives `Re` to
ten times the target and `Δt` to 0.03 cells. `run!` now aborts at `|v| > 4 U`.
Log: `logs/phaseB_andersen_thin_b14_L80_dive.log`.

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

1. **Flow direction and the convective exit.** WaterLily's `exitBC=true`
   is a first-order *upwind* convection of the exit profile at the mean inflow
   speed, valid only for flow towards `+x`. In the first Phase B runs the grid
   was oriented so that the plate fell towards `+x`, i.e. the relative flow
   ran from `xmax` to `0`: the scheme became downwind and the coupled run died
   at the exit face within ca. 4 `c/U`, while the wake hit a Dirichlet wall
   three chords behind the plate (a `χ = 14` plate then dived edge-on to
   `5 U`). With the grid oriented as documented in `src/FallingPlate2D.jl`
   (fall towards `−x`, relative flow `+x`) WaterLily's exit is stable in the
   coupled 1-DOF fall. The default `exit=:own` (`update_exit!`) reproduces the
   same scheme through the `uBC` function but freezes the profile when the
   far-field speed becomes negative, which happens briefly when a fluttering
   plate rises at the end of a glide; `exit=:waterlily` and `exit=:dirichlet`
   remain available for comparison. An earlier note here blamed a
   predictor–corrector flux inconsistency; that was wrong.
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

### Phase B — 3-DOF free fall — **DONE 12 September 2026**

**Validation result** (`examples/free_fall.jl 8 2.52 2090 64 60`, ellipse
`χ = 8`, `ρ* = 2.52`, `I* = 0.160`, 514 × 514 cells, 8 cells across the
thickness, 67 ms/step, aborted by the `|v| > 4U` guard at 42.6 `c/U` during a
fast glide; means over 12.6–42.6 `c/U`, log `logs/phaseB_andersen_b8_L64.log`,
data `out/andersen_dns_b8_L64.csv`, figure `.png` alongside):

| quantity | this code | Andersen et al. table 1, `β = 1/14` |
|---|---|---|
| regime | tumbling, 2.7 revolutions in 30 `c/U` | tumbling |
| `⟨ω⟩ c/|⟨vy⟩|` | 0.90 | 0.857 |
| `Re = c|⟨vy⟩|/ν` | 1290 | 1147 |
| `⟨vx⟩/|⟨vy⟩|` | 1.24 | 2.20 |
| `⟨vy⟩/U` | −0.62 | −0.55 |

The tumbling is more irregular than in the experiment: phases of steady
tumbling alternate with long glides reaching `2–3 U`, and the plate briefly
rises at the end of some glides. Andersen's trajectories are periodic with
`vy` between 0.3 and 0.9 `v_g` and no rising. Likely causes, to be tested in
that order: (i) resolution — the boundary layer is ca. 2 cells thick at
`L = 64`, `Re = 1300`; repeat at `L = 96` and `128` (ca. 1 h and 2.5 h per
60 `c/U`); (ii) time step — `Δt_max = 0.5` cells while the plate is slow;
(iii) `Re_g` — lower it so the realised `Re` is 1147; (iv) the periodic
lateral domain of 8 chords when the plate drifts 20 chords. The user's
OpenFOAM results (Pandi et al.) are the quantitative reference for the same
case. The `χ = 14` literal plate is discussed in section 6.
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

### Phase C — permeability — **implemented 12 September 2026** (`src/darcy.jl`)

The implicit Brinkman update `u(1 + Δtσ) = u* − Δt∇p + ΔtσV`, `σ = ν/K`, has
exactly the form of WaterLily's BDIM blend `u = μ₀(u* − Δt∇p) + (1 − μ₀)V`
with `μ₀ = 1/(1 + Δtσ)`, and the same coefficient is the pressure Poisson
mobility. So a permeable plate is the BDIM body with its interior mobility
raised from 0 to `μ_p = 1/(1 + Δt ν/(Da c²))`, blended through the kernel:
`μ₀ ← μ₀ + (1 − μ₀) μ_p`, `μ₁ = 0`. `Da → 0` recovers the impermeable body
exactly. Nothing else changes: no sub-iteration, one Poisson iteration per
step as before. Two details that matter: the domain-boundary faces must be
reset to zero mobility afterwards (otherwise the pressure solver diverges to
NaN), and the interior velocity fill must test the SDF, not just `μ₀ == 0`.

Held broadside plate, `χ = 4`, `L = 32`, `Re = 200`, 4 `c/U`
(`scripts/phaseC_held.jl`):

| `Da` | `C_D` | seepage `u/U` at the centre | screening length `c√Da` / Δ |
|---|---|---|---|
| 0 | 1.55 | 0 | 0 |
| 10⁻⁵ | 1.52 | 0.009 | 0.1 |
| 10⁻⁴ | 1.51 | 0.08 | 0.3 |
| 10⁻³ | 0.94 | 0.55 | 1.0 |
| 10⁻² | 0.14 | 0.95 | 3.2 |

The trend matches Cummins et al. (2017): the drag collapses once the
screening length exceeds the cell size, i.e. the interesting range
`Da = 10⁻⁴`–`10⁻²` is resolved at `Δ = c/32`, while `Da ≤ 10⁻⁵` is
numerically indistinguishable from impermeable at this resolution (as the
resolution-floor argument above predicts). Force on a permeable plate is still
the BDIM surface integral; for a thin plate the pressure jump across it is the
Darcy drag, so this is consistent, but a momentum-balance force
`−∑(1 − μ₀)(V − u* + Δt∇p)/Δt` would be cleaner and is the first Phase C+ task.

Toy free-fall sweep (`scripts/phaseC_sweep.jl`, `χ = 4`, `ρ* = 2.7`,
`Re_g = 500`, `L = 32`, means over 10–20 `c/U`):

| `Da` | `⟨vy⟩/U` | `⟨vx⟩/U` | `⟨ω⟩ c/U` | `Re` on `⟨vy⟩` | motion |
|---|---|---|---|---|---|
| 0 | −0.53 | −0.18 | −0.38 | 266 | fluttering |
| 10⁻⁴ | −0.68 | 0.02 | 0.04 | 340 | steady, broadside |
| 10⁻³ | −1.54 | 0.00 | 0.00 | 770 | steady, broadside |
| 10⁻² | −2.11 | −0.02 | 0.00 | 1054 | steady; run stopped at 6.4 `c/U` by the `|v| > 4U` guard, still accelerating (values are over 0–6.4 `c/U`) |

Permeability kills the flutter already at `Da = 10⁻⁴` and the descent speed
then rises steeply with `Da`, as for permeable disks (Cummins et al. 2017).
At `Da = 10⁻²` the plate is nearly transparent and has no terminal speed
within 20 `c/U` at this `Re_g`.

To do: free-fall sweeps in `Da` at fixed `I*`, `Re` (the research question);
a check against Bose et al. (2025) / Cummins et al. (2017) drag curves;
momentum-balance force.

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

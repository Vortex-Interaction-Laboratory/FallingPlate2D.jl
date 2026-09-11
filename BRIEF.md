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
- **Repository skeleton already exists at `~/Developer/FallingPlate2D.jl`**
  (git initialised on `main`, one commit: MIT licence + `.gitignore`).
  Deliberately **outside OneDrive** — git indexes and OneDrive sync corrupt
  each other, and simulation output must never sync.

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

| Mode | Setup | Grid | Expected |
|---|---|---|---|
| **Toy** | `χ = 4`, `Δ = c/32`, domain 8c × 16c | 256 × 512 | 2–4 min |
| **Check** | `χ = 10`, `Δ = c/80`, domain 8c × 16c | 640 × 1280 | 20–40 min |

These are estimates on 4 threads, **to be measured, not trusted**. Thin plates
are where BDIM costs you; a thicker section is what buys the toy-mode
turnaround.

Use `T = Float64`.

---

## 8. WaterLily API — verify before building

Confirmed present in
[WaterLily-Examples](https://github.com/WaterLily-jl/WaterLily-Examples):

- `TwoD_FreeRotatingEllipse.jl` — ellipse free to rotate about a pivot under
  fluid forces. **Closest starting point: your problem with translation locked.**
- `TwoD_CircleVIV.jl` — vortex-induced vibration, body ODEs via OrdinaryDiffEq.jl.
- `TwoD_SlowStartCircle.jl` — accelerating reference frame, time-dependent BCs.
- Force functions: `total_force(sim)`, `pressure_force(sim)`, `viscous_force(sim)`.
- `g(i,t)` keyword — spatially uniform, time-varying body force. This is both
  gravity and the non-inertial frame term.

**Uncertain:** `RigidMap` (storing linear/angular position and velocity, updated
via `update!`) appears in current documentation but **not** in the 2024 CPC
paper, so it is a recent addition. **Verify it exists in the installed version
before designing around it.**

---

## 9. Build plan

### Phase A — verify the API (~1–2 h)
Add WaterLily to the project (`julia -t 4 --project=.` then
`] add WaterLily`), run `TwoD_FreeRotatingEllipse.jl` unmodified, confirm the
force functions, the `g(i,t)` keyword, and whether `RigidMap` exists. Report
what the installed version actually provides. **This phase validates the whole
plan — do it first and report before proceeding.**

### Phase B — 3-DOF free fall (~1 day, highest risk)
- Body signed distance function for an ellipse/plate parameterised by `χ`.
- **Translating frame:** keep the body near the domain centre; apply the
  negative frame acceleration as a uniform body force through `g(i,t)`; set the
  far-field inflow to minus the body velocity. This is what makes long
  trajectories affordable — the plate falls O(100) chords, and a lab-frame
  domain would be absurd.
- Rigid-body integration: solve `(M + M_a)·a = F_fluid + (M − ρ_f·V)·g`.
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

## 10. Planned repository layout

```
FallingPlate2D.jl
├── BRIEF.md          # this file
├── src/bodies.jl     # SDFs: plate, ellipse, perforated plate, parameterised by χ
├── src/rigid.jl      # 3-DOF body (x, y, θ): gravity, buoyancy, implicit added mass
├── src/darcy.jl      # Brinkman resistance in the BDIM forcing, parameterised by Da
├── src/diagnostics.jl
├── examples/
└── test/             # validation against Andersen et al. (2005)
```

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

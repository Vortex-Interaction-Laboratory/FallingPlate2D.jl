# FallingPlate2D.jl

A fast 2D toy solver for a freely falling plate section in a quiescent fluid,
impermeable or permeable (Darcy–Brinkman), built as a thin layer over
[WaterLily.jl](https://github.com/WaterLily-jl/WaterLily.jl). Reynolds numbers
10–1000, chord-to-thickness ratios up to ca. 14, runs on a laptop in minutes.

It is an exploration tool. For quantitative claims the group's validated
OpenFOAM overset setup (Pandi, Gungor, Bose, Attili & Viola,
[arXiv:2304.04317](https://arxiv.org/abs/2304.04317)) is the reference.

## What it does

- Solves the incompressible Navier–Stokes equations with WaterLily's boundary
  data immersion method in a frame translating with the plate, so the plate
  can fall hundreds of chords in an 8 × 8 chord domain (periodic sides,
  convective exit).
- Integrates the three rigid-body degrees of freedom with the added mass and
  added moment of inertia treated implicitly, which removes the sub-iteration
  and the low-density-ratio instability of weakly coupled schemes.
- Makes the plate permeable by raising the immersed-boundary mobility to the
  implicit Brinkman value inside the body; the Darcy number is one parameter
  and `Da = 0` recovers the impermeable plate exactly.

## Install and run

Julia ≥ 1.10. From the repository root:

```bash
julia -t 4 --project=. -e 'using Pkg; Pkg.instantiate()'
julia -t 4 --project=. examples/free_fall.jl 8 2.52 2090 64 60 myrun 0
```

Arguments: chord-to-thickness ratio `χ`, density ratio `ρ*`, gravitational
Reynolds number `Re_g = c·sqrt(2(ρ*−1) f g)/ν`, cells per chord `L`, duration
in convective times `c/U`, output name, Darcy number `Da = K/c²`. The run
writes `out/myrun.csv` (nondimensional time series of position, angle,
velocities, loads) and prints a ten-line comparison with Andersen, Pesavento
& Wang (2005). Plot it with

```bash
julia --project=examples examples/plot_run.jl out/myrun.csv
```

Use `-t 4` on an Apple M2 (performance cores only); Float64 on the CPU.

## Conventions

Lengths in chords `c`, velocities by `U = sqrt(2(ρ*−1) f g)`, fluid density 1.
Lab frame `(X, Y)` with `Y` upwards; the chord line makes the angle `φ` with
the horizontal, positive anticlockwise, `φ = 0` broadside to the fall. Forces
are nondimensionalised by `U²c`, so coefficients are twice the listed values.
Frames and grid orientation are documented in `src/FallingPlate2D.jl`.

## Validation so far

Ellipse `χ = 8`, `ρ* = 2.52` (dimensionless moment of inertia `I* = 0.16`, the
configuration of Andersen et al.'s direct numerical simulations), 64 cells per
chord, means over 30 convective times:

| quantity | this code | Andersen et al. (2005), table 1, β = 1/14 |
|---|---|---|
| regime | tumbling | tumbling |
| ⟨ω⟩ c / \|⟨v_y⟩\| | 0.90 | 0.857 |
| Re on ⟨v_y⟩ | 1290 | 1147 |
| ⟨v_x⟩ / \|⟨v_y⟩\| | 1.24 | 2.20 |

The tumbling is more irregular than in the experiment; a resolution study is
in progress. Details, measured runtimes and the open questions are in
`BRIEF.md`.

## Two things to know about WaterLily for moving bodies

Both are handled in `src/simulation.jl` and documented in `BRIEF.md` §8:

1. `measure!` fills the imposed body velocity only within three cells of the
   surface. A moving body thicker than that gets a discontinuous interior
   velocity, a spurious interior pressure of O(10⁴), and forces of the wrong
   sign. `rigid_velocity_fill!` fixes it.
2. The convective exit is upwind for flow towards `+x` only. Orient the
   problem accordingly, or use the exit provided here, which also handles a
   momentarily reversed far field.

## Layout

```
src/bodies.jl       sections (ellipse, rounded plate): SDF, area, inertia, added mass, I*
src/rigid.jl        rigid-body state and implicit-added-mass update
src/simulation.jl   coupled simulation: frame, exit, interior fill, run!, CSV output
src/darcy.jl        Darcy–Brinkman mobility
src/diagnostics.jl  CSV reader and summary
examples/           free_fall.jl driver, plot_run.jl, upstream WaterLily example
scripts/            the diagnostic scripts behind the findings in BRIEF.md
test/               fast checks (< 1 min)
```

## References

- Andersen, Pesavento & Wang, *J. Fluid Mech.* **541** (2005) 65–90.
- Weymouth & Font, WaterLily.jl, [arXiv:2407.16032](https://arxiv.org/abs/2407.16032).
- Cummins et al., *J. Fluid Mech.* (2017), permeable disks; Bose et al. (2025), permeable plates.

## Licence

MIT.

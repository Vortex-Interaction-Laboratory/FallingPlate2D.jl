# Response of a free-falling permeable plate to a modulated gravitational force

Ignazio Maria Viola · 14 September 2026 · prepared with Claude Opus 5

## 1. Objective

To measure how a transient modulation of the gravitational force displaces a
freely falling permeable plate, and whether that displacement reverses when the
order of the modulation is reversed. The case is a two-dimensional plate at a
Galilei number $Ga = 150$ and a Darcy number $Da = 10^{-4}$.

## 2. Method

### 2.1 Scaling

Lengths are scaled on the plate length $c$, velocities on the gravitational
velocity $u_g \equiv (\hat m \hat g/\hat\rho \hat c^2)^{1/2}$, and the fluid
density is unity. The governing groups are the Galilei number
$Ga \equiv u_g c/\nu$, the nondimensional mass $m$, the chord-to-thickness
ratio $\chi \equiv c/f$, the nondimensional moment of inertia $I_O$ about the
plate centre, the centre-of-gravity offset $e$ along the body normal, and the
Darcy number $Da \equiv K/c^2$, where $K$ is the permeability. The solid volume
of a filamentous body is neglected, so buoyancy is null and the nondimensional
gravitational acceleration is $g = 1/m$.

Values are those of the established D2 case, the two-dimensional equivalent of
the three-dimensional permeable disk: $\chi = 10$, $m = 0.80$, $e = 1.16$,
$I_O = 1.15$, hence $g = 1.25$. Only $Ga$ and $Da$ were set by the present
study.

### 2.2 Flow solver

The incompressible Navier–Stokes equations are solved with WaterLily.jl, which
uses the boundary data immersion method on a uniform Cartesian grid. The
equations are solved in a frame translating with the plate centre, so that the
plate remains at a fixed grid location and only its angle changes; the frame
acceleration enters as a uniform body force, and the far-field velocity equals
minus the plate velocity. The lateral direction is periodic and the downstream
boundary is a convective exit.

### 2.3 Permeability

In a Brinkman medium the fluid within the solid obeys
$$\partial_t \bm{u} + \ldots = -\nabla p + \nu\nabla^2\bm{u} - \sigma(\bm{u} - \bm{V}), \quad \sigma = \nu/K, \tag{1}$$
where $\bm{V}$ is the solid velocity. Treated implicitly over a step, this is
exactly the blend of the boundary data immersion method with a mobility
$\mu = (1 + \Delta t\,\sigma)^{-1}$. A permeable plate is therefore obtained by
raising the interior mobility from zero to $\mu_p$, which also enters the
pressure equation, so that pressure acts on the pore fluid with the reduced
mobility. The limit $Da \to 0$ recovers the impermeable body exactly.

### 2.4 Rigid-body coupling

Three degrees of freedom are integrated, with the added-mass tensor treated
implicitly, so that
$$(m\bm{I} + \bm{M}_a)\,\bm{a} = \bm{F} + \bm{M}_a\bm{a}_\text{old} + m\,\bm{g}_\text{eff}. \tag{2}$$
The fixed point of this relaxation is exact and its error amplification is
$|\bm{M}_{a,\text{num}} - \bm{M}_a|/(m + \bm{M}_a)$. Implicit treatment removes
the sub-iteration that otherwise costs a factor of 10 to 20 for a freely
falling body.

### 2.5 Gravity modulation

The modulation is applied as a transient change of the gravitational
acceleration, $\bm{g}_\text{eff} = \bm{g} - \dot{\bm{u}}_G$, over a single
period $T_G = 2$ convective times beginning at $t_0$. With a semi-amplitude $A$
the effective gravity follows
$$g_\text{eff}(t)/g = 1 \mp A\sin\!\left[2\pi (t - t_0)/T_G\right], \tag{3}$$
so that over one period it passes through $(1-A)g$, $g$, $(1+A)g$ and back to
$g$ for the upper sign, and the reverse order for the lower sign. The net
impulse is null. Unlike a fluid gust, the fluid is not displaced.

For $A = 0.5$ this is a vertical isolated gust of ratio
$GR = \mp A g T_G/\pi = \mp 0.3979$, with a peak
$|\dot u_G| = 0.625 = 0.5\,g$.

### 2.6 Restart

Restart files hold the velocity and pressure fields, the time-step history, the
rigid-body state and the frame bookkeeping. The immersed-boundary fields are
functions of the body position and are rebuilt on load. This allows several
branches to be run from one common state.

### 2.7 Numerical parameters

The plate length is 64 cells, giving 6.4 cells across the thickness, in a
domain of $16c \times 16c$, i.e. a grid of $1026^2$ including ghost cells. The
kinematic viscosity is $\nu = 6.667\times 10^{-3}$ in grid units. Time
integration is semi-implicit Euler for the body and the solver's own adaptive
step for the flow, capped at 0.5 cells.

## 3. Verification

Three checks were made.

**Restart fidelity.** A run advanced to $t = 0.3$, checkpointed, and continued
to $t = 0.6$ was compared with a fresh simulation restarted from that
checkpoint and advanced to the same time. Position, velocity, angle, angular
velocity and elapsed time agreed to machine zero: $\max|\Delta X| = 0$. The
restart is therefore bitwise equivalent, and the branches differ only by the
imposed gravity.

**Modulation extremes.** Equation (3) was verified to reach $0.5g$, $g$, $1.5g$
and $g$ at quarter-period intervals to a tolerance of $10^{-12}$, with null
gust velocity at the end of the period.

**Terminal state.** Convergence was assessed from the cycle-mean descent speed
over successive pendulum periods (section 4.1).

## 4. Results

### 4.1 The terminal state is a limit cycle

The plate does not approach a constant descent speed. Because the centre of
gravity is offset, it develops a pendulum oscillation whose amplitude grows and
then saturates. The measured period is $5.8\,c/u_g$, against the estimate
$2\pi(mge/I_O)^{-1/2} = 6.26$.

| window $t\,u_g/c$ | mean $v_Y/u_g$ | $\varphi$ semi-amplitude [rad] |
|---|---|---|
| 8.0–13.8 | −1.1397 | 0.036 |
| 19.6–25.4 | −1.0451 | 0.055 |
| 31.2–37.0 | −1.0127 | 0.080 |
| 37.0–42.8 | −1.0071 | 0.086 |
| 42.8–48.6 | −1.0077 | 0.087 |

The last two cycle means differ by $5.7\times10^{-4}$, i.e. 0.06 %, and the
amplitude has flattened. The terminal state is therefore a saturated limit
cycle with a cycle-mean descent speed of $-1.008\,u_g$ (figure
`lap1_approach.png`).

An earlier branch attempted at $t = 20$, where $v_Y = -1.054$, was discarded:
the base state was still drifting by 6 % over five convective times, which
would have contaminated the response.

### 4.2 Response to the modulation

Three branches were run from the common state at $t_0 = 50$, where
$Y/c = -51.006$, $v_Y = -1.0125$, $\varphi = +0.0452$ and $\omega = -0.1011$,
each for ten convective times.

| branch | $Y/c$ at $t = 60$ | $\Delta Y/c$ | $v_Y/u_g$ at $t = 60$ |
|---|---|---|---|
| constant $g$ | −61.0745 | — | −1.01420 |
| $g$ decreasing first | −61.0331 | $+0.0414$ | −1.01482 |
| $g$ increasing first | −61.1185 | $-0.0439$ | −1.01285 |

The plate whose gravity is reduced first ends higher, and the plate whose
gravity is raised first ends lower (figure `gravity_three.png`). The separation
from the constant-gravity plate is monotonic in sign throughout: it grows
during the first half-period, peaks shortly after the modulation reverses, at
$+0.1124$ at $t = 51.37$ and $-0.1068$ at $t = 51.35$ respectively, then
relaxes to a residual that is constant thereafter. The peak separation is
therefore some 2.6 times the lasting one.

All three final velocities agree to within $2\times10^{-3}$, confirming that
the modulation carried no net impulse.

### 4.3 The response is close to antisymmetric

Reversing the order of the modulation reverses the sign of the displacement
almost exactly:
$$\Delta Y(\text{decreasing first}) + \Delta Y(\text{increasing first}) = -0.0025\,c, \tag{4}$$
against individual responses of some $0.042\,c$. The displacement is therefore
linear in the signed semi-amplitude to within 6 %, and only that residual is
genuinely nonlinear.

This excludes rectification by quadratic drag as the leading mechanism, which
would give the same sign for both orderings. It is also not the response of a
simple first-order relaxation: for $\dot v = -v/\tau + f$ with $\int f\,dt = 0$
one obtains $\int v\,dt = 0$ and hence no lasting displacement. A non-zero
lasting displacement therefore requires memory, either in the wake and the
added mass or through coupling to the pendulum mode. Which of these dominates
has not been established.

## 5. Limitations

**Pendulum phase.** The base state is a limit cycle, so the response depends on
the phase at which the modulation begins. All results above are at one phase,
$\varphi = +0.045$, $\omega = -0.101$. The reported displacements are one
sample of a phase-dependent quantity, and the 6 % residual in equation (4) may
be phase sensitivity rather than nonlinearity.

**Brinkman resolution.** The screening length is $c\sqrt{Da} = 0.64$ cells at
$Da = 10^{-4}$ and $L = 64$, so the Brinkman layer within the plate is not
resolved. For a thin plate the relevant physics is the lumped pressure drop
$\sigma f u_\text{seep}$, which does not require it, but no resolution study of
the permeable case was made here.

**Discretisation.** The solver is second order and the immersed boundary smears
the interface over some two cells, with 6.4 cells across the plate thickness.
No grid convergence study was made for the displacements reported in section
4.2. Their sign and relative magnitude are robust, since the branches are
deterministic from a common state; their absolute values are not established.

## 6. Files

Simulations were run on an Apple M2 with four performance cores. The lap cost
$829 + 2090 = 2919$ s of wall time; each branch cost 400 s on four threads, or
861 s when two shared the machine.

| file | content |
|---|---|
| `checkpoints/lap1_Da1e-4_t50.jld2` | terminal state, the common branch point |
| `results/lap_Da1e-4_t0-20.csv`, `…_t20-50.csv` | the lap |
| `results/end1.csv` | constant gravity |
| `results/end2.csv` | gravity decreasing first |
| `results/end3.csv` | gravity increasing first |
| `figures/lap1_approach.png` | approach to the limit cycle |
| `figures/gravity_three.png` | the three branches |
| `figures/gravity_lap.png` | earlier two-branch comparison |

Code: <https://github.com/ignaziomviola/FallingPlate2D.jl>, commit `b69d69b`.
Drivers `scripts/gravity_lap.jl` and `examples/plot_gravity_three.jl`.

## 7. References

Andersen, A., Pesavento, U. & Wang, Z. J. 2005 Unsteady aerodynamics of
fluttering and tumbling plates. *J. Fluid Mech.* **541**, 65–90.

Bose, C. & Viola, I. M. 2023 (gust model and permeable-plate convention).

Pandi, J. S. S., Gungor, A., Bose, C., Attili, A. & Viola, I. M. 2026 Effect of
transverse gust on free-falling plates. arXiv:2304.04317.

Weymouth, G. D. et al. 2025 WaterLily.jl: a differentiable and
backend-agnostic Julia solver. *Comput. Phys. Commun.*

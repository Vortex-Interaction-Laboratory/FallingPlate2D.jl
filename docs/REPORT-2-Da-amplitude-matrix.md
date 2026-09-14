# Permeability and the lasting displacement of a free-falling plate under a modulated gravitational force

Ignazio Maria Viola · 14 September 2026 · prepared with Claude Opus 5

## 1. Objective

To measure how the response of a freely falling permeable plate to a transient
modulation of the gravitational force depends on the Darcy number, the
semi-amplitude of the modulation, and the order in which the modulation is
applied. Twenty-one simulations were run at a Galilei number $Ga = 150$.

## 2. Method

The method is that of Report 1, to which the reader is referred for the
scaling, the flow solver, the Brinkman permeability model, the implicit
added-mass coupling and the restart procedure. The plate is the D2 case,
$\chi = 10$, $m = 0.80$, $e = 1.16$, $I_O = 1.15$, on a grid of $1026^2$ with a
plate length of 64 cells.

The matrix spans three Darcy numbers, $Da = 10^{-5}, 10^{-4}, 10^{-3}$; three
semi-amplitudes, $A = 0.5, 1.5, 2.5$; and two orderings, in which the
gravitational force either falls to $(1-A)g$ first and then rises to $(1+A)g$
("down first"), or the reverse ("up first"). The modulation occupies a single
period $T_G = 2$ convective times and carries no net impulse. For $A \ge 1.5$
the effective gravity reverses sign during part of the cycle, reaching $-0.5g$
and $-1.5g$ respectively. A constant-gravity reference was run for each Darcy
number, giving $3 \times (1 + 3 \times 2) = 21$ simulations.

**Each Darcy number required its own lap**, because the terminal state depends
on permeability. Two new laps of 50 convective times were run for
$Da = 10^{-5}$ and $10^{-3}$, and the existing state was reused for $10^{-4}$.
All branches start from the terminal state of their own Darcy number at
$t_0 = 50$ and run for ten convective times. No run triggered the stability
guards.

## 3. Results

### 3.1 Permeability sets the character of the terminal state

Measured from the constant-gravity runs over $50 \le t \le 60$:

| $Da$ | mean $v_Y/u_g$ | $\varphi$ semi-amplitude [rad] | $\max|\omega| c/u_g$ | character |
|---|---|---|---|---|
| $10^{-5}$ | −0.937 | 0.503 | 0.665 | fluttering |
| $10^{-4}$ | −1.007 | 0.088 | 0.115 | pendulum limit cycle |
| $10^{-3}$ | −2.110 | 0.0000 | 0.0000 | steady |

Permeability both **stabilises the descent and accelerates it**. Over two
decades of $Da$ the plate passes from large-amplitude fluttering, through a
small saturated pendulum oscillation, to a perfectly steady descent, while the
mean descent speed rises in magnitude from 0.937 to 2.110 $u_g$. This is
consistent with the stabilising role of permeability reported experimentally,
and with the fall in drag coefficient with $Da$ recorded previously for this
solver.

### 3.2 The peak excursion is linear in the semi-amplitude

The largest separation from the constant-gravity plate, reached shortly after
the modulation reverses:

| $Da$ | $A = 0.5$ | $A = 1.5$ | $A = 2.5$ | peak$/A$ |
|---|---|---|---|---|
| $10^{-4}$ | 0.112 | 0.354 | 0.613 | 0.22, 0.24, 0.25 |
| $10^{-3}$ | 0.254 | 0.762 | 1.274 | 0.51, 0.51, 0.51 |

For $Da = 10^{-3}$ the peak excursion is linear in $A$ to three significant
figures, despite the effective gravity reversing sign for $A \ge 1.5$. For
$Da = 10^{-4}$ it is linear to within 12 %. The more permeable, faster-falling
plate responds about twice as strongly during the modulation.

### 3.3 Permeability decides whether the displacement lasts

The transient excursion and the lasting displacement behave quite differently.
Residual $\Delta Y/c$ at $t = 60$, eight convective times after the modulation
ends:

| $Da$ | $A$ | peak | residual, down first | residual, up first | residual/peak |
|---|---|---|---|---|---|
| $10^{-3}$ | 0.5 | 0.254 | $+0.0015$ | $-0.0004$ | 0.006 |
| $10^{-3}$ | 1.5 | 0.762 | $+0.0090$ | $+0.0016$ | 0.012 |
| $10^{-3}$ | 2.5 | 1.274 | $+0.0264$ | $+0.0075$ | 0.021 |
| $10^{-4}$ | 0.5 | 0.112 | $+0.0414$ | $-0.0439$ | 0.37 |
| $10^{-4}$ | 1.5 | 0.354 | $+0.0972$ | $-0.0737$ | 0.27 |
| $10^{-4}$ | 2.5 | 0.613 | $+0.1198$ | $-0.1428$ | 0.20 |

At $Da = 10^{-3}$ the plate **forgets the perturbation almost completely**: the
excursion decays from 0.25 c at the end of the modulation to 0.0015 c eight
convective times later, and the two orderings converge on the same trajectory.
At $Da = 10^{-4}$ between a fifth and a third of the excursion **persists
indefinitely**, with the sign set by the ordering.

The distinction follows the character of the base state. A steady, strongly
damped descent relaxes back onto its original trajectory; a limit cycle retains
the perturbation in its phase, and that memory does not decay.

### 3.4 Antisymmetry and its limits

At $Da = 10^{-4}$ reversing the ordering reverses the sign of the residual. The
residual sum, which would vanish for a response linear in the signed
semi-amplitude, is

| $A$ | 0.5 | 1.5 | 2.5 |
|---|---|---|---|
| sum $/c$ | $-0.0025$ | $+0.0234$ | $-0.0231$ |
| as a fraction of the mean response | 6 % | 27 % | 18 % |

So the response is close to antisymmetric at the smallest amplitude and departs
from it as the amplitude grows, as expected once gravity reverses. The
sensitivity $\Delta Y/A$ **falls** with amplitude — 0.083, 0.065, 0.048 — so
the lasting displacement saturates rather than growing in proportion.

At $Da = 10^{-3}$ the transient is antisymmetric to better than 1 %
($+0.2089$ against $-0.2084$ at $A = 0.5$), but the residual is too small for
the comparison to be meaningful.

### 3.5 The $Da = 10^{-5}$ results are not interpretable

They are reported for completeness but should not be used.

| $A$ | residual, down first | residual, up first | sum |
|---|---|---|---|
| 0.5 | $-0.0743$ | $+0.0656$ | $-0.0087$ |
| 1.5 | $+0.5773$ | $+0.5664$ | $+1.1437$ |
| 2.5 | $+0.6171$ | $+0.9312$ | $+1.5483$ |

At $A \ge 1.5$ **both orderings displace the plate the same way**, by an amount
an order of magnitude larger than at any other Darcy number, and $\Delta Y(t)$
never settles. The base state flutters with a semi-amplitude of 0.50 rad, so
any perturbation shifts the phase of the flutter and the measured separation is
dominated by the divergence of two trajectories with different phase, not by
the response to the modulation. A deterministic response cannot be extracted
from a single realisation in this regime.

Note also that at $A = 0.5$ the residual has the **opposite sign** to that at
$Da = 10^{-4}$, which is a further indication that the quantity being measured
is not the same one.

## 4. Limitations

**Brinkman resolution.** The screening length $c\sqrt{Da}$ is 2.02, 0.64 and
0.20 cells at $Da = 10^{-3}, 10^{-4}, 10^{-5}$ for a plate of 64 cells. Only
the largest Darcy number resolves it. For a thin plate the governing physics is
the lumped pressure drop across the plate, which does not require the internal
layer to be resolved, but the two lower Darcy numbers should be regarded as
under-resolved. Resolving $Da = 10^{-5}$ would need a plate of at least 320
cells, some 25 times the cost.

**Single phase.** Every branch starts at one phase of its base state. For
$Da = 10^{-4}$ the residuals are therefore one sample of a phase-dependent
quantity, and the 6–27 % departures from antisymmetry in section 3.4 may be
phase sensitivity rather than nonlinearity. Separating the two requires a phase
sweep, which was not run.

**No grid convergence study.** The signs and the relative ordering of the
results are robust, since all branches are deterministic continuations of a
common state and differ only in the imposed gravity. The absolute values are
not established.

**Gravity reversal.** For $A \ge 1.5$ the effective gravity is negative for
part of the cycle. This is a mathematically well-posed perturbation of the
body force but has no obvious physical realisation.

## 5. Computational cost

All runs on an Apple M2 with four performance cores. The two new laps took 4828
and 5100 s, run concurrently on two threads each, finishing in 85 minutes. The
21 branches ran sequentially on four threads at some 500 s each, taking a
further 3 hours. Total wall time about 4.5 hours; no run required
intervention.

## 6. Files

`results/m_<Da>_const.csv` and `results/m_<Da>_A<A>_<down|up>.csv` for the 21
runs; `results/mlap_<Da>.csv` and `checkpoints/mlap_<Da>.jld2` for the laps;
`figures/matrix.png`.

Code: <https://github.com/ignaziomviola/FallingPlate2D.jl>. Drivers
`scripts/matrix.jl`, `scripts/chain_matrix.sh`, `examples/plot_matrix.jl`.

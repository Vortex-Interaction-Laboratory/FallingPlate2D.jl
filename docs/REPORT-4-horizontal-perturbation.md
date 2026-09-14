# Horizontal versus vertical gravity perturbation of a free-falling plate

Ignazio Maria Viola · 14 September 2026 · prepared with Claude Opus 5

## 1. Objective

Report 3 applied the perturbation to the weight itself. Here the perturbation
is applied **horizontally**, the vertical component of gravity being held
constant at $g$, and the amplitudes of the two responses are compared. Twelve
simulations were run with the angular position constrained, at Darcy numbers
$10^{-5}$ and $10^{-4}$, semi-amplitudes $A = 0.5, 1.5, 2.5$, and both signs.

## 2. Method

As Report 3, with the perturbation direction changed from $\bm{e}_Y$ to
$\bm{e}_X$, so that the effective gravity is
$$\bm{g}_\text{eff}(t) = \left(\pm A g \sin\!\left[2\pi(t-t_0)/T_G\right],\; -g\right), \tag{1}$$
over one period $T_G = 2$ convective times. The vertical component is constant
throughout and the net horizontal impulse is null. The sign convention is
"pos" for a force along $+X$ first, "neg" for the reverse.

The base state does not depend on the perturbation, so the terminal
checkpoints and the constant-gravity references of Report 3 were reused
unchanged. No new laps were needed and only twelve runs were required. All
twelve completed and end at $t = 40.00$.

Because the base state is exactly vertical, with $X/c$ and $v_X$ zero to five
decimal places, the horizontal response is measured directly rather than as a
small difference between two finite quantities.

## 3. Results

### 3.1 The horizontal response is much larger

| $Da$ | $A$ | peak, horizontal | peak, vertical | ratio | lasting, horizontal | lasting, vertical | ratio |
|---|---|---|---|---|---|---|---|
| $10^{-5}$ | 0.5 | 0.425 | 0.095 | 4.5 | 0.150 | 0.010 | 15 |
| $10^{-5}$ | 1.5 | 1.260 | 0.307 | 4.1 | 1.082 | 0.049 | 22 |
| $10^{-5}$ | 2.5 | 3.427 | 0.538 | 6.4 | 3.413 | 0.118 | 29 |
| $10^{-4}$ | 0.5 | 0.364 | 0.109 | 3.3 | 0.075 | 0.017 | 4.5 |
| $10^{-4}$ | 1.5 | 1.059 | 0.350 | 3.0 | 0.489 | 0.084 | 5.8 |
| $10^{-4}$ | 2.5 | 1.690 | 0.613 | 2.8 | 1.015 | 0.167 | 6.1 |

The same forcing produces a peak excursion **three to six times larger** when
applied horizontally, and a lasting displacement **four to twenty-nine times
larger**. At $Da = 10^{-4}$, $A = 2.5$ the plate is permanently displaced
sideways by more than one chord.

### 3.2 The character of the response also changes

Decomposing as in Report 3 into a part whose sign follows the ordering and a
part that does not:

| $Da$ | $A$ | antisymmetric | rectified |
|---|---|---|---|
| $10^{-4}$ | 0.5 | $+0.0750$ | $+0.0004$ |
| $10^{-4}$ | 1.5 | $+0.4892$ | $+0.0066$ |
| $10^{-4}$ | 2.5 | $+1.0150$ | $+0.0000$ |
| $10^{-5}$ | 0.5 | $+0.1496$ | $-0.0441$ |
| $10^{-5}$ | 1.5 | $+1.0820$ | $-0.1136$ |
| $10^{-5}$ | 2.5 | $+3.4126$ | $+0.0144$ |

The horizontal response is **antisymmetric**, at $Da = 10^{-4}$ to four decimal
places: reversing the sign of the perturbation reverses the displacement
exactly, and the rectified part vanishes. This is the opposite of the vertical
case with the same constraint, where the response was rectified and the
antisymmetric part was small.

### 3.3 The lasting fraction approaches unity

The ratio of the lasting displacement to the peak excursion:

| $Da$ | $A = 0.5$ | $A = 1.5$ | $A = 2.5$ |
|---|---|---|---|
| $10^{-4}$ | 0.21 | 0.46 | 0.60 |
| $10^{-5}$ | 0.35 | 0.86 | 1.00 |

At $Da = 10^{-5}$, $A = 2.5$ the lasting displacement equals the peak to three
significant figures: the plate is pushed sideways and stays there, with
essentially nothing recovered. The less permeable plate retains more, at
$Da = 10^{-5}$ some three times the displacement of the $10^{-4}$ plate.

### 3.4 Interpretation

The vertical and horizontal directions are not equivalent for this body, for
two reasons that both act in the same sense.

**There is no equilibrium horizontally.** Vertically the plate sits at a
terminal velocity, a stable fixed point set by the balance of weight and drag,
so a perturbation is opposed and decays back towards it. Horizontally the plate
has no preferred position and no restoring force at all; the direction is
neutrally stable. A perturbation with zero net impulse therefore leaves an
offset that only drag can erode, and sections 3.2 and 3.3 show how little of it
is eroded.

**The plate is slender in the direction of the forcing.** Held broadside, it
presents its thin edge to horizontal motion. For an ellipse of semi-axes $a$
and $b$ the added masses are $m_{11} = \rho\pi b^2$ along the chord against
$m_{22} = \rho\pi a^2$ normal to it, a ratio of $\chi^2 = 100$ at
$\chi = 10$; and the resisting drag is correspondingly smaller.

Together these account for both observations: the larger amplitude, because the
plate is easier to move sideways and nothing pushes it back; and the
antisymmetry, because a displacement about a neutrally stable direction with
little damping stays close to linear in the signed forcing, whereas the
vertical response is a perturbation about a stable equilibrium whose nonlinear
restoring law rectifies it.

**These statements are consistent with the measurements but are not
demonstrated by them.** No force decomposition, quasi-steady baseline or
flow-field data was collected, so the mechanism is inferred from the
displacement record alone.

## 4. Limitations

Those of Report 3 carry over: the Brinkman screening length is unresolved at
both Darcy numbers, and no grid-convergence study was made. Two further points
are specific to this campaign.

**The lateral direction is periodic** and the frame translates with the plate,
so a lateral displacement of more than one chord is not obstructed by the
domain. It does, however, mean that the plate at $Da = 10^{-5}$, $A = 2.5$ has
moved 3.4 chords across a domain 16 chords wide, and any interaction with its
own periodic image is untested.

**Superlinearity.** The antisymmetric response per unit amplitude grows with
$A$ — 0.150, 0.326, 0.406 at $Da = 10^{-4}$ — whereas the peak excursion per
unit amplitude is almost constant at 0.73, 0.71, 0.68. The lasting fraction,
not the excursion, is what grows with amplitude.

## 5. Computational cost

Twelve runs of ten convective times, sequential on four threads, 438 s each,
88 minutes in total. Reusing the laps and references of Report 3 saved a
further hour.

## 6. Files

`results/h_<Da>_A<A>_<pos|neg>.csv`; `figures/horiz.png`. Base states and
references are those of Report 3.

Code: <https://github.com/ignaziomviola/FallingPlate2D.jl>. Driver
`scripts/matrix_horiz.jl`, chain `scripts/chain_horiz.sh`, figure
`examples/plot_horiz.jl`. `gravity_modulation` now takes a `dir` argument.

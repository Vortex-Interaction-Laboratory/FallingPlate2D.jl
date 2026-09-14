# Where the lasting displacement lives: a free-falling plate with its rotation constrained

Ignazio Maria Viola · 14 September 2026 · prepared with Claude Opus 5

## 1. Objective

Reports 1 and 2 established that a gravity modulation of zero net impulse
leaves a freely falling plate permanently displaced, but left open whether that
memory resides in the rotational degree of freedom, as the phase of a pendulum
oscillation, or in the translational one, through the wake and the added mass.

The present campaign answers that by repeating the matrix with the angular
position constrained: the plate is held horizontal, $\varphi = 0$, and cannot
rotate. Fourteen simulations were run, at Darcy numbers $10^{-5}$ and $10^{-4}$
and semi-amplitudes $A = 0.5, 1.5, 2.5$, for both orderings, with a
constant-gravity reference for each Darcy number.

## 2. Method

The method is that of Reports 1 and 2, with `free = (true, true, false)` and
$\varphi_0 = 0$. Locking the rotation removes the pendulum mode at
$Da = 10^{-4}$ and the flutter at $10^{-5}$ entirely; both plates descend
exactly vertically, with $X/c$ and $v_X$ zero to five decimal places
throughout. The centre-of-gravity offset $e = 1.16$ is retained for consistency
but is inert when $\omega = \alpha = 0$.

Because the terminal state differs, two new laps were run, of 30 convective
times rather than 50: without the pendulum this is a simple drag relaxation and
saturates far sooner. Branches start at $t_0 = 30$ and run for ten convective
times. No run triggered the stability guards and all fourteen completed.

## 3. Results

### 3.1 Constraining the rotation accelerates the descent

| $Da$ | mean $v_Y/u_g$, locked | free-rotating | change |
|---|---|---|---|
| $10^{-5}$ | −1.303 | −0.937 | 39 % faster |
| $10^{-4}$ | −1.285 | −1.007 | 28 % faster |

Holding the plate broadside removes the drag penalty of fluttering and of the
pendulum swing. More striking, the two locked cases descend at very nearly the
same speed, −1.303 and −1.285, whereas the free-rotating ones differed by 7 %.
**Over this decade of $Da$ the direct effect of permeability on the descent
speed is therefore small; the difference seen in Report 2 came almost entirely
from the different oscillation amplitudes, not from the drag of the permeable
plate itself.**

### 3.2 The response becomes rectified

With the rotation locked, **both orderings displace the plate in the same
direction**, upwards, at every Darcy number and every amplitude:

| $Da$ | $A$ | peak | residual, down first | residual, up first |
|---|---|---|---|---|
| $10^{-5}$ | 0.5 | 0.095 | $+0.0037$ | $+0.0158$ |
| $10^{-5}$ | 1.5 | 0.307 | $+0.0425$ | $+0.0555$ |
| $10^{-5}$ | 2.5 | 0.538 | $+0.1233$ | $+0.1131$ |
| $10^{-4}$ | 0.5 | 0.109 | $+0.0220$ | $+0.0116$ |
| $10^{-4}$ | 1.5 | 0.350 | $+0.1076$ | $+0.0610$ |
| $10^{-4}$ | 2.5 | 0.613 | $+0.2111$ | $+0.1230$ |

This is the opposite of the free-rotating case, where at $Da = 10^{-4}$ the two
orderings gave residuals of opposite sign.

### 3.3 The two mechanisms separate cleanly

Decomposing the residual into a part whose sign follows the ordering and a part
that does not,
$$\Delta Y_\text{anti} \equiv \tfrac{1}{2}(\Delta Y_\downarrow - \Delta Y_\uparrow), \qquad
  \Delta Y_\text{rect} \equiv \tfrac{1}{2}(\Delta Y_\downarrow + \Delta Y_\uparrow), \tag{1}$$
gives, at $Da = 10^{-4}$:

| $A$ | free, antisym. | free, rectified | locked, antisym. | locked, rectified |
|---|---|---|---|---|
| 0.5 | $+0.0427$ | $-0.0013$ | $+0.0052$ | $+0.0168$ |
| 1.5 | $+0.0854$ | $+0.0117$ | $+0.0233$ | $+0.0843$ |
| 2.5 | $+0.1313$ | $-0.0115$ | $+0.0441$ | $+0.1670$ |

**The two components are almost perfectly complementary.** With the rotation
free, the response is antisymmetric and the rectified part is null to within
0.012 c. With the rotation locked, the antisymmetric part collapses by a factor
of three to eight and the rectified part becomes dominant, growing to 0.167 c.

The conclusion is that the lasting displacement has two distinct origins:

1. **An antisymmetric contribution carried by the rotational degree of
   freedom.** The modulation shifts the phase of the pendulum oscillation, and
   the sign of that shift follows the ordering. Lock the rotation and it
   disappears.
2. **A rectified contribution carried by the translation.** The plate ends
   higher whichever way the modulation runs, because during the low-gravity
   excursion it descends more slowly and its quadratic drag falls faster than
   linearly, so the high-gravity excursion does not recover the lost distance.
   This is the mechanism I proposed in Report 1 and then rejected on the
   evidence of the free-rotating runs; the present campaign shows it was
   present all along, merely masked by the larger rotational contribution.

### 3.4 Scaling of the rectified part

The rectified residual grows faster than linearly with amplitude — at
$Da = 10^{-4}$, 0.017, 0.084, 0.167 for $A = 0.5, 1.5, 2.5$, so roughly as
$A^{1.5}$ — consistent with a quadratic-drag origin, which cannot produce a
term linear in the perturbation. It also grows with permeability, from 0.118 to
0.167 at $A = 2.5$ between $Da = 10^{-5}$ and $10^{-4}$.

### 3.5 The $Da = 10^{-5}$ case is now interpretable

In Report 2 these runs were unusable because the base state fluttered at a
semi-amplitude of 0.50 rad and the measured separation was dominated by
flutter-phase divergence. With the rotation locked there is no flutter, the
residuals are smooth and monotonic in $A$, and the antisymmetric component is
null to within 0.007 c at every amplitude. This confirms the earlier diagnosis:
the anomalous $Da = 10^{-5}$ numbers in Report 2 measured trajectory
divergence, not a response.

## 4. Limitations

**Brinkman resolution** is unchanged from Report 2: the screening length is
0.64 and 0.20 cells at $Da = 10^{-4}$ and $10^{-5}$ for a 64-cell plate, so
neither resolves the internal layer.

**No grid convergence study.** The separation of mechanisms in section 3.3 rests
on a comparison between campaigns run at identical resolution, so it is robust
to a common discretisation error; the absolute values are not established.

**The rectified scaling in section 3.4 rests on three points** and an exponent
inferred from them should be treated as indicative.

**A single decade of $Da$.** No locked run was made at $Da = 10^{-3}$, where
the free-rotating plate descends at 2.11 $u_g$ and forgets the perturbation
almost entirely. Adding it would test whether the rectified component also
vanishes in the strongly permeable limit.

## 5. Computational cost

Two laps of 30 convective times, run concurrently on two threads each, and
fourteen branches run sequentially on four threads at some 400 s each. Total
about 2 h 50 min on an Apple M2, interrupted once by a machine sleep of some 35
minutes, which suspended and resumed the runs without loss.

A buffering defect was found and fixed during this campaign: `run!` opened its
output file but never flushed it, so a live run's CSV lagged behind its true
progress and an interrupted run would have been truncated at an arbitrary
buffer boundary. The handle is now flushed after every row.

## 6. Files

`results/f_<Da>_const.csv` and `results/f_<Da>_A<A>_<down|up>.csv`;
`results/flap_<Da>.csv` and `checkpoints/flap_<Da>.jld2`; `figures/fixed.png`.

Code: <https://github.com/ignaziomviola/FallingPlate2D.jl>. Drivers
`scripts/matrix_fixed.jl`, `scripts/chain_fixed.sh`, `examples/plot_fixed.jl`.

## 7. Next step

To apply the modulation **horizontally** rather than vertically and compare the
amplitudes of the two responses. The perturbation direction is `Gust.dir` in
`src/gust.jl`; `gravity_modulation` currently hard-codes the vertical direction
and needs a direction argument first.

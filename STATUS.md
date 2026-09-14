# STATUS — 14 September 2026, 19:20

Current state for a fresh session. Read this instead of the previous
conversation. Background and build history are in `BRIEF.md`; results are in
`docs/REPORT-1..4.md` and `report/report.tex` (Report 5, the narrative one).

---

## 1. What is running right now

`scripts/chain_ga.sh` (nohup, `caffeinate -ims` bound to its PID): Galilei-number
sweep at `Da=1e-5`, semi-amplitude `A=2.5`, rotation locked.

- 2 laps to `t=30` (`Ga=50`, `Ga=100`), then per `Ga`: a constant-gravity
  reference plus four branches (`v down`, `v up`, `h pos`, `h neg`), `t=30→40`.
- Output `out/g_Ga{50,100}_*.csv`, logs `logs/g_Ga*.log`, chain log
  `logs/chain_ga.log`. Finishes on `GA MATRIX DONE`.
- Launched 19:13, ETA ~21:15 plus any lid-close sleep.
- **Resume after any interruption:** `nohup ./scripts/chain_ga.sh >> logs/chain_ga.log 2>&1 &`
  (completed runs are skipped only if their CSV ends at `t=40`).

### FIRST ACTION for a new session: arm a waiter

The waiter armed in the previous session died with it, so nothing will report
when the campaign finishes. Arm one before doing anything else, as a
**background** Bash call so it notifies on completion:

```bash
cd ~/Developer/FallingPlate2D.jl && until grep -qE "GA MATRIX DONE|exit=[1-9]" logs/chain_ga.log 2>/dev/null; do sleep 60; done; echo "GA CAMPAIGN COMPLETE $(date +%T)"; grep -c "end   " logs/chain_ga.log
```

Note the pattern matches a non-zero exit as well as success: a watcher that
greps only for the success marker stays silent through a crash, and silence
then looks identical to still running. Do the same for any future campaign.
If the campaign has already finished by the time the new session starts, skip
this and go straight to the analysis.

Deliberately uses the same lap length and branch window as the `Ga=150`
campaign, for comparability rather than convergence. See section 3.

---

## 2. Campaigns completed

All at `χ=10, m=0.80, e=1.16, I_O=1.15` (the D2 plate), `L=64`, domain `16c×16c`.

| set | prefix | config | branch at |
|---|---|---|---|
| 1 | `m_*` | `Ga=150`, free rotation, `Da∈{1e-5,1e-4,1e-3}`, `A∈{0.5,1.5,2.5}`, both orders, vertical | `t=50` |
| 2 | `f_*` | as above but **rotation locked**, `Da∈{1e-5,1e-4}` | `t=30` |
| 3 | `h_*` | as set 2 but **horizontal** perturbation | `t=30` |

Headline results, all in the reports:
1. Free to rotate: the plate gains **or** loses altitude, sign set by the order
   of the fluctuation (antisymmetric response).
2. Rotation locked: **every** case gains altitude (rectified response).
3. Horizontal fluctuation: **always** a gain, 13–90× larger than vertical at
   `t=40`.

---

## 3. IMPORTANT — the base states were never converged

**This is not yet written in any report and it qualifies Reports 3, 4 and 5.**

- The **locked** laps were run to `t=30` on the argument that removing the
  pendulum leaves a simple drag relaxation. They had **not** converged: `vY`
  was still rising in magnitude monotonically (−1.147, −1.195, −1.222, −1.243,
  −1.263 at `t=10,15,20,25,30` for `Da=1e-4`), and the constant-gravity
  reference still drifts by −0.022 over the last four convective times of its
  own window.
- Consequence: "terminal state" is the wrong word throughout Reports 3–5. The
  branches are compared at a fixed time after release, not at a converged
  condition.
- What survives: every branch is an exact continuation of its own reference
  from an identical checkpoint, so the common drift largely cancels and the
  `ΔY` comparisons stand.
- What does not: the deficit analysis. Statements that the horizontally
  perturbed plate is "still descending 7–38 % more slowly" compare two
  quantities that are both still evolving; exponential fits to the deficit came
  out inconsistent for this reason.
- The user's physical argument, which is correct: a transverse perturbation of
  zero net impulse **cannot** change a unique terminal state, so `vY` must
  recover. Demonstrating that needs a converged reference, which we do not have.
- A lap extension to `t=120` was launched and then **stopped** on the user's
  instruction. `ckpt/flap_{1e-5,1e-4}.jld2` (t=30) are intact and resumable:
  `julia -t 2 --project=. scripts/matrix_fixed.jl extend <Da> <tend>`.

For the free-rotation set the base state **was** verified converged (cycle-mean
`vY` changing by 0.06 % between the last two pendulum periods at `t≈43–49`).

---

## 4. Defects found and fixed today

- `run!` opened its CSV but never flushed it, so a live run's file lagged and an
  interrupted run would have been truncated at an arbitrary buffer boundary.
  Now flushed after every row (`src/simulation.jl`).
- Chain scripts tested only that an output file was non-empty, so a partial file
  would be skipped rather than redone. `chain_ga.sh` now requires the CSV to end
  at `t=40`. **`chain_matrix.sh` and `chain_fixed.sh` still have the weak test.**
- `matrix_horiz.jl` saved no end checkpoint, so the `h_*` runs cannot be
  extended and would have to be re-run from `t=30`. Fixed for `matrix_horiz.jl`
  and `matrix_ga.jl`; `matrix_fixed.jl` branches still do not checkpoint.

---

## 5. Known limitations to carry into any report

- **Brinkman layer unresolved**: `c√Da` is 2.02, 0.64 and 0.20 cells at
  `Da=1e-3, 1e-4, 1e-5` for `L=64`. Only `1e-3` resolves it.
- **No grid-convergence study** anywhere in this work.
- **No mechanism evidence.** Only trajectories and integrated loads were
  recorded — no pressure/viscous force split, no quasi-steady baseline, no
  flow fields, no circulation. Any statement about wakes, vortices, separation
  or quasi-steady behaviour is inference, not measurement. Reports 1 and 3
  contain a quadratic-drag attribution that is an argument from scaling only.
- **`Da=1e-5` with free rotation is uninterpretable**: the base state flutters
  at 0.50 rad semi-amplitude, so the measured separation is flutter-phase
  divergence, not a response.
- **Locked campaigns span one decade of `Da`** (1e-5 to 1e-4), not the full
  range explored with free rotation.
- In the largest horizontal case the plate moves 3.4 chords laterally in a
  16-chord domain with periodic lateral boundaries; image interaction untested.

---

## 6. Where things live

- Repo `~/Developer/FallingPlate2D.jl`, remote `github.com/ignaziomviola/FallingPlate2D.jl`,
  `main` up to date (`bade23f`).
- OneDrive `Edinburgh/AI/CFD free falling 2D plate/`: `REPORT-1..4.md`,
  `REPORT-5-altitude-response.pdf`, `reports/` (LaTeX source + figure PDFs),
  `results/`, `checkpoints/`, `figures/`, `logs/`.
- `out/` and `ckpt/` are git-ignored; OneDrive holds the copies.

---

## 7. Next steps, in the order they were agreed

1. Analyse the `Ga` sweep when it lands; add the two points to the narrative.
2. Decide whether to state the non-convergence in Reports 3–5 or re-run the
   locked campaign from a converged base state (~8 h for the horizontal set).
3. Outstanding from earlier: six locked runs at `Da=1e-3` (~90 min) would
   extend claim 2 beyond one decade of `Da`.

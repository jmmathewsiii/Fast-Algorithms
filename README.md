# 2D N-Body Simulation

A 2D gravitational N-body simulator with three runnable variants:

| Target            | Algorithm    | Backend |
|-------------------|--------------|---------|
| `run-naive`       | Direct O(N²) | Serial CPU |
| `run-serial-BH`   | Barnes-Hut O(N log N) | Serial CPU |
| `run-CUDA-direct` | Direct O(N²) | CUDA GPU |

> **Note:** A CUDA Barnes-Hut variant is **not implemented yet**.

All variants use a 2nd-order kick-drift-kick (KDK) symplectic Verlet integrator with Plummer softening.

> **All `make` and run commands below assume you are inside the `n_body/` subdirectory.** From the repository root:
>
> ```bash
> cd n_body
> ```

---

## Requirements

- **g++** with C++14 support (any reasonably recent GCC works)
- **GNU make**
- **gnuplot** — for viewing the resulting plots/animations
- **CUDA toolkit (nvcc)** — *only* required for `run-CUDA-direct`. The two serial targets build without CUDA.
- A CUDA-capable GPU at runtime for `run-CUDA-direct`. The code was tested on an A100 (Alpine cluster) but should run on any compute-capability ≥ 6.0 device.

If you don't have CUDA, `make run-naive` and `make run-serial-BH` still work fine on their own.

---

## Building & running

From inside `n_body/`:

```bash
make run-naive          # builds + runs serial direct
make run-serial-BH      # builds + runs serial Barnes-Hut
make run-CUDA-direct    # builds + runs CUDA direct
```

Each target prompts you for three values at startup. Just hit `enter` to accept the default:

```
Number of particles [100]:
Number of iterations [1000]:
Mode (speed / plot) [plot]:
```

- **Number of particles** — `N`. Sensible range: 100 to ~10⁵ for direct, up to ~10⁶ for BH on CPU, more for CUDA.
- **Number of iterations** — total integrator steps. With `dt = 0.001`, 1000 steps = 1 time unit.
- **Mode** — `speed` or `plot`:
  - **speed**: prints periodic energy values during the run and total/force/tree-build timings at the end. No `.plt` or animation files are emitted.
  - **plot**: dumps per-frame `.plt` files and gnuplot scripts to `output/` and `scripts/`. No timings are printed. (This is the default.)

---

## Plotting results (plot mode)

`plot` mode writes a particle animation script and (for BH) a quadtree animation script. To view:

```bash
cd scripts
gnuplot
gnuplot> load 'direct.gp'        # particle animation, direct run
gnuplot> load 'bh.gp'            # particle animation, BH run
gnuplot> load 'bh_tree.gp'       # quadtree animation, BH run
```

Hit `q` to close a window. The scripts read their data from `../output/*.plt`.

---

## Lyapunov / BH-error analysis

When you run **both** `run-naive` and `run-serial-BH` in plot mode, the BH run produces an error-vs-time CSV plus a gnuplot script that fits the Lyapunov exponent.

Workflow:

```bash
make run-naive          # mode = plot, accept defaults (or pick larger N/n_iter)
make run-serial-BH      # mode = plot, with the SAME N and n_iter
```

Both runs must use the same `N`, `n_iter`, and seed (the seed is hard-coded). At the end of the second run you'll see something like:

```
[validate] wrote 101 rows to ./output/lyapunov/error_bh_vs_direct_serial.csv
[validate] wrote gnuplot script: ./scripts/lyapunov/error_bh_vs_direct_serial.gp
```

To view the error-growth plot and Lyapunov fit:

```bash
cd scripts/lyapunov
gnuplot
gnuplot> load 'error_bh_vs_direct_serial.gp'
```

The script plots all three error metrics on a semilog-y axis and fits the exponential-growth region. The fitted slope is the largest Lyapunov exponent λ. **Tune the fit window** — open the `.gp` file and edit `fit_step_lo` and `fit_step_hi` to bracket the visually-linear (in log-y) region after re-running gnuplot.

The Lyapunov files (~one `.dat` per snapshot stride per run, plus the CSV) live in `output/lyapunov/` so they don't clutter the main `output/` directory.

---

## Speed-mode validation

In `speed` mode, each non-baseline run compares its final particle positions against the `direct_serial` baseline saved on disk. Run the baseline first:

```bash
make run-naive          # mode = speed
make run-serial-BH      # mode = speed   → prints abs/rel error vs direct_serial
make run-CUDA-direct    # mode = speed   → prints abs/rel error vs direct_serial
```

If the baseline file is missing, the comparison is skipped with a warning.

---

## Running on the CU Alpine cluster

The repo ships with two SLURM batch files (`slurm_cpu.bat`, `slurm_gpu.bat`) for running on Alpine.

**1. Compile on a compile node.** SSH into Alpine, then from a compile node:

```bash
cd n_body
module load gcc/11.2.0           # required for the serial targets
module load cuda                 # additionally required for cuda_direct

make serial_direct serial_bh     # CPU targets
make cuda_direct                 # GPU target
```

(You can also just run `make` to build all three.)

**2. Submit to a compute node** with `sbatch` (still from inside `n_body/`):

```bash
sbatch slurm_cpu.bat             # runs serial_direct then serial_bh on amilan-style node
sbatch slurm_gpu.bat             # runs cuda_direct on an A100 node (atesting_a100)
```

Both batch files are wired for **speed mode** and feed `N`, `n_iter`, and `speed` to the binaries via a heredoc. To change problem size, edit the `N=` and `NITER=` lines at the top of the relevant `.bat` file. Output goes to `logs/<job>-<jobid>.out`.

Each batch file's `#SBATCH` block (partition, time, GPU request, etc.) is preserved as-is — adjust those lines if your account/QOS needs differ.

---

## Cleaning

```bash
make clean       # removes obj/ and bin/
make cleanall    # also clears output/ and scripts/
```

---

## Directory layout

```
Fast-Algorithms/
├── README.md             (this file)
└── n_body/               <-- everything below lives here
    ├── bin/              built executables
    ├── include/          headers
    ├── obj/              build artifacts
    ├── output/           final-position dumps + plot data
    │   └── lyapunov/     Lyapunov snapshot dumps + error CSV
    ├── scripts/          generated gnuplot scripts
    │   └── lyapunov/     Lyapunov fit script
    ├── src/              C++ / CUDA sources
    │   └── kernels/      force-evaluation kernels
    ├── Makefile
    ├── slurm_cpu.bat
    └── slurm_gpu.bat
```

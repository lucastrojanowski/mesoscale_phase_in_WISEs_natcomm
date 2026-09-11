This repository includes all of the code/figures used in our numerical simulations of phase separation via the Cahn-Hilliard equation detailed in [insert_doi_here]. Please see the publication for a full discussion of our findings associated with this repository.

Notably, the contents of this repository can be run as-is to reproduce Figures 3 f,g, and h as well as Fig S22 from the associated publication. A minimal usage entails cloning the repository, creating a conda environment (details below), installing the prerequisite python packages into the environment, and simply running the respective cahn-hilliard.py scripts to generate the data for Fig 3 f, g, and h and/or Fig S22. It is not necessary to run create-ic.py to create initial conditions for the simulations since the initial conditions already exist, but for completeness one can re-generate the initial conditions used in the publication. No edits to the existing scripts are necessary. Please find a thorough discussion of the code and other contents of this repository below.

================================================================================
Cahn-Hilliard mesoscale aggregation simulations
================================================================================

This directory reproduces Fig. 3 f, g, and h as well as Fig. S22.

Layout:

  initial-conditions/create-ic.py           - generates initial conditions
                                               (.npy arrays + quick-look .pdf
                                               heatmaps), saved in this same
                                               directory
  fig_3/cahn-hilliard_fig3.py               - runs the Fig. 3 simulation
                                               (order_param = 0.2)
  fig_s22/cahn-hilliard_fig_s22.py          - runs the Fig. S22 simulation
                                               (order_param = 0.6)

Both simulation scripts load their initial condition from
../initial-conditions/ (relative to their own location) and write their
outputs (.mp4/.pdf/.png) next to themselves, regardless of the directory
you run them from — i.e. `python fig_3/cahn-hilliard_fig3.py` from anywhere
lands outputs in fig_3/, and equivalently for fig_s22/.

--------------------------------------------------------------------------------
Environment
--------------------------------------------------------------------------------

These scripts were run with the following environment:

  Python        3.14.7
  numpy         2.5.3
  matplotlib    3.11.1
  imageio       2.37.4
  imageio-ffmpeg 0.6.0
  pillow        12.3.0

imageio-ffmpeg bundles its own ffmpeg binary, so no separate system-wide
ffmpeg install is required as long as imageio-ffmpeg is installed.

To reproduce this environment:

  conda create -n <env_name> python=3.14
  conda activate <env_name>
  pip install -r requirements.txt

(requirements.txt pins the exact versions listed above.)

Both scripts should also run unmodified on a wider range of reasonably
recent versions (Python >= 3.9, numpy >= 1.24, matplotlib >= 3.6), but the
versions above are what was actually used to generate the figures/movies
currently in this directory. The only hard requirement beyond numpy/
matplotlib is that `imageio` (with an ffmpeg-capable backend, e.g.
imageio-ffmpeg) be installed, since it is used to write the .mp4 output.

--------------------------------------------------------------------------------
1. Generating initial conditions: initial-conditions/create-ic.py
--------------------------------------------------------------------------------

create-ic.py builds 250x250 initial condition arrays for the Cahn-Hilliard
order parameter field u (which takes values in [-1, 1], where
-1 = pure salt, +1 = pure solvent). The solvent volume fraction c (reported in the publication) is related
to u by c = (u+1)/2, which maps [-1, 1] to [0, 1].

The __main__ block sweeps a list of `order_params` (uniform mean value of
u) and `deltas` (amplitude of uniform random noise added on top), and for
each combination calls:

  create_uniform_mix(path, order_param, delta, dim=250, seed=12041997)

which saves both:
  - {order_param}-{delta}.npy  (the raw array, u0)
  - {order_param}-{delta}.pdf  (a quick-look heatmap)
into the same directory as the script (initial-conditions/).

As shipped, the __main__ block regenerates exactly the two initial
conditions actually consumed by the simulation scripts in this directory
(order_param 0.2 for fig_3, 0.6 for fig_s22; delta 0.01 for both). Because
the random seed is fixed, re-running it reproduces these two .npy files
byte-for-byte. However, other initial conditions could easily be generated.

Run it (from anywhere):

  python initial-conditions/create-ic.py

Edit the `order_params` / `deltas` lists at the bottom of the file to
generate additional initial conditions.

--------------------------------------------------------------------------------
2. Running simulations: fig_3/cahn-hilliard_fig3.py, fig_s22/cahn-hilliard_fig_s22.py
--------------------------------------------------------------------------------

Both scripts integrate:

  du/dt = D * [ Laplacian(u^3 - u) - gamma * Laplacian(Laplacian(u)) ]

on a periodic grid using an explicit finite-difference scheme, starting
from an initial condition u0, via the shared core function:

  solve_cahn_hilliard(diff_coeff, dt, N, grid_size, save_path, gamma,
                       snapshot_steps, u0=None, mp4_freq=100)

  - diff_coeff     : diffusion coefficient D
  - dt             : integration time step
  - N              : number of time steps to integrate
  - grid_size      : side length of the square grid (must match u0's shape)
  - save_path      : output path for the .mp4 movie
  - gamma          : interfacial energy penalty
  - snapshot_steps : which integer time-step indices (t, not t*dt) to pull
                      out and save as labeled panels
  - u0             : initial condition array (2D numpy array, grid_size x
                      grid_size); if None, a uniform random field is used
  - mp4_freq       : write one movie frame every this many steps

The function raises a ValueError if the chosen (diff_coeff, dt, gamma,
grid_size) combination violates the explicit-scheme stability conditions
(Dγdt/dx^2 <= 0.01 and Ddt/dx^2 <= 0.1, checked at the top of the
function), so start from the example values below and adjust dt/grid_size
together if you change diff_coeff or gamma.

Run them (from anywhere):

  python fig_3/cahn-hilliard_fig3.py
  python fig_s22/cahn-hilliard_fig_s22.py

Parameters actually used, and outputs produced, per script:

  fig_3/cahn-hilliard_fig3.py
    dt=0.01, N=100000, diff_coeff=1, gamma=1, grid_size=250, mp4_freq=10
    order_param=0.2 -> loads initial-conditions/0.2-0.01.npy
    snapshot_steps=(0, 3990, 28980) -> titles "t=0.00", "t=39.90", "t=289.80"
    writes, next to the script:
      fig3_0.6.mp4                  (0.6 = solvent volume fraction (0.2+1)/2)
      fig3_0.6_t0.pdf/.png, fig3_0.6_t3990.pdf/.png, fig3_0.6_t28980.pdf/.png
      fig3_0.6_colorbar.pdf/.png

  fig_s22/cahn-hilliard_fig_s22.py
    dt=0.01, N=100001, diff_coeff=1, gamma=1, grid_size=250, mp4_freq=100
    order_param=0.6 -> loads initial-conditions/0.6-0.01.npy
    snapshot_steps=(0, 7100, 100000) -> titles "t=0.00", "t=71.00", "t=1000.00"
    writes, next to the script:
      figs22_0.8.mp4                (0.8 = solvent volume fraction (0.6+1)/2)
      figs22_0.8_t0.pdf/.png, figs22_0.8_t7100.pdf/.png, figs22_0.8_t100000.pdf/.png
      figs22_0.8_colorbar.pdf/.png
      figs22_0.8_combined.pdf/.png  (fully labeled side-by-side panel figure
                                     with all three snapshots and a shared
                                     colorbar; fig_3's script does not
                                     produce this extra combined figure)

Notes:
  - grid_size must match the initial condition's array shape (250 for both
    initial conditions used here).
  - snapshot_steps are given as step indices, not simulation time; the
    corresponding simulation time is t_index * dt, which is what gets
    printed in the per-frame panel titles.

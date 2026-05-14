# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System

Copy a config template from `config/` to `config.mk` before building:

```bash
# Serial build
cp config/unix-gcc-default.mk config.mk
make serial && make

# Parallel (MPI) build
cp config/unix-mpicc-default.mk config.mk
make build

# Current config: CUDA + PETSc (RTX 4070, sm_89)
make build          # builds src/Ludwig.exe and libludwig.a
make test           # unit tests + d3q19-short regression tests
make unit           # unit tests only
make clean
```

The active `config.mk` currently targets NVCC with PETSc (`HAVE_PETSC=true`), OpenMPI, and CUDA sm_89. Debug build (`-O0 -g -G`). For production add `-DNDEBUG`.

Run unit tests directly:
```bash
LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH ./tests/unit/a.out
```

Run regression tests by model:
```bash
make -C tests d3q19-short   # fastest suite
make -C tests d2q9
make -C tests d3q27
```

## Architecture

Ludwig is a lattice Boltzmann solver for complex fluids (mixtures, colloids, liquid crystals, electrokinetics). The main simulation loop is in `src/ludwig.c`.

### Core flow (`ludwig_run`)
1. **Colloid update** — rebuild lattice links if needed
2. **Order parameter gradients** — halo exchange + gradient stencils
3. **Electrokinetics** (if `psi` active):
   - Scatter charges from particles/fluid onto lattice (subgrid kernels)
   - Solve Poisson equation (FFT via cuFFT, or PETSc, or SOR)
   - Compute electrostatic forces on fluid (`psi_force_gradmu`) and particles (`subgrid_update_Esub`, `subgrid_update_forces_electrokinetics`)
   - Nernst-Planck transport (multistep)
4. **Free energy forces** — stress divergence or grad-mu
5. **LB collision + propagation**
6. **Bounce-back** on colloid links

### Key source modules

| File | Purpose |
|------|---------|
| `src/ludwig.c` | Main driver and time-stepping loop |
| `src/collision.c` | LB collision operators (largest file) |
| `src/psi.c` / `psi.h` | Electrostatic potential and charge density |
| `src/psi_fft.c` / `psi_fft_pn.c` | FFT-based Poisson solver (cuFFT) |
| `src/psi_force.c` | Force on fluid from electrostatics |
| `src/subgrid.c` / `subgrid.h` | Subgrid particle scatter/gather kernels for charge and force |
| `src/ewald_charge.c` | Ewald summation for charges |
| `src/nernst_planck.c` | Nernst-Planck ion transport |
| `src/colloids.c` | Colloid data structures |
| `src/build.c` | Build lattice-colloid links |
| `src/bbl.c` | Bounce-back on links |

### Subgrid kernels
Particle-to-mesh and mesh-to-particle operations use regularized delta functions selected by `subgrid_kernel_t`:
- `SUBGRID_KERNEL_PESKIN4` — Peskin 4-point kernel
- `SUBGRID_KERNEL_BSPLINE4` / `BSPLINE6` — B-spline kernels
- `SUBGRID_KERNEL_KB4` — Kaiser-Bessel kernel
- `SUBGRID_KERNEL_PESKIN6` — Peskin 6-point kernel

The active kernel is set by `kernel_g` in `ludwig.c`.

### Interlacing mode (reciprocal-space PME)
Controlled by `interlacing_mode_` (static in `ludwig.c`):
- `0` — standard path (no interlacing)
- `1` — single shifted grid
- `2` — average of two grids (full interlacing, equivalent to P3M/PME deconvolution)

### FFT Poisson solver
`psi_fft_pn_t` (`psi_fft_pn.c`) solves ∇²ψ = -ρ/ε in reciprocal space using cuFFT. Deconvolution to correct for the smearing kernel is controlled by `psi_fft_deconv_t` (`deconv_g`).

### Code modification conventions
All local modifications are wrapped with `/*CHANGE INIT - <date/description> */` ... `/*CHANGE END - ... */` markers. Do not delete commented-out code inside these blocks; it records the original behavior.

### Free energy models
Selected at runtime via input file. Key types: `FE_ELECTRO`, `FE_SYMMETRIC`, `FE_LC_DROPLET`, `FE_TERNARY`, `FE_BRAZOVSKII`. Polymorphic interface via `fe_t*` with function pointer table.

### GPU/CPU abstraction
TargetDP (`target/`) provides `tdp*` wrappers for CUDA/HIP/OpenMP. Data lives on both host and device; use `field_memcpy`, `hydro_memcpy`, `lb_memcpy`, `map_memcpy`, `colloids_memcpy` to synchronize.

### Parallelism
MPI domain decomposition via `cs_t` (coordinate system). Halo exchanges: `psi_halo_psi`, `psi_halo_rho`, `psi_halo_psijump`, `hydro_u_halo`, `field_halo`.

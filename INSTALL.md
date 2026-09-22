# Instalación de Ludwig — dependencias y entorno

Esta guía describe, en orden, todo lo que hay que instalar para compilar y
ejecutar Ludwig con la configuración actual del repo (`config.mk`):
**CUDA + PETSc (con hypre) + OpenMPI**, más el entorno Python para los
scripts de análisis/gráficos (`microgel/`, `util/`, `docs/figures/`).

Referencia de la máquina en la que este `config.mk` fue probado:
- GPU: NVIDIA RTX 4070 (sm_89), driver 595.x
- CUDA Toolkit 12.8 (`nvcc`)
- OpenMPI compilado en `/usr/local/ompi`
- PETSc compilado en `/usr/local/petsc-cuda-hypre` (variante con hypre/BoomerAMG en GPU)

Si se usa otra GPU, cambiar `sm_89` por la arquitectura correspondiente en
`config.mk` (`-arch=sm_XX` en `CFLAGS` y `LDFLAGS`).

## Orden de instalación

1. Toolchain base (gcc, make, git, etc.)
2. Driver NVIDIA + CUDA Toolkit
3. OpenMPI (compilado con soporte CUDA/UCX)
4. PETSc (compilado contra ese OpenMPI y CUDA, con hypre)
5. Ludwig (compilar con `config.mk`)
6. Entorno Python (para scripts de post-proceso, independiente del build de Ludwig)

El orden 2 → 3 → 4 importa: OpenMPI debe compilarse con soporte CUDA para que
PETSc pueda usar GPU-aware MPI, y PETSc debe compilarse contra ese mismo
OpenMPI y contra el mismo CUDA Toolkit.

---

## 1. Toolchain base

```bash
sudo apt update
sudo apt install -y build-essential gfortran git cmake pkg-config \
    python3 python3-venv python3-pip wget
```

## 2. Driver NVIDIA + CUDA Toolkit

- Instalar el driver NVIDIA correspondiente a la GPU (via `apt` o el
  instalador `.run` de NVIDIA). Verificar con:
  ```bash
  nvidia-smi
  ```
- Instalar el CUDA Toolkit (12.x). Opción recomendada: paquete oficial de
  NVIDIA (`cuda-toolkit-12-8` o similar) o el instalador `.run` desde
  https://developer.nvidia.com/cuda-downloads.
- Verificar:
  ```bash
  nvcc --version
  ```
- Anotar la arquitectura de la GPU (`sm_XX`) — se necesita para `config.mk`.
  Para RTX 4070 es `sm_89`. Se puede consultar con:
  ```bash
  nvidia-smi --query-gpu=name,compute_cap --format=csv
  ```

## 3. OpenMPI (con soporte CUDA)

Ludwig usa MPI para el paralelismo de dominio y PETSc lo necesita también.
Se recomienda compilar OpenMPI desde fuente con soporte CUDA para permitir
transferencias GPU-aware (UCX):

```bash
# Dependencias de UCX (opcional pero usado por LAUNCH_MPIRUN_CMD en config.mk)
sudo apt install -y libucx-dev ucx-utils

wget https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.x.tar.gz
tar xzf openmpi-5.0.x.tar.gz
cd openmpi-5.0.x
./configure --prefix=/usr/local/ompi \
    --with-cuda=/usr/local/cuda \
    --with-ucx
make -j$(nproc)
sudo make install
```

Agregar al `PATH`/`LD_LIBRARY_PATH` (por ejemplo en `~/.bashrc`):

```bash
export PATH=/usr/local/ompi/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/ompi/lib:$LD_LIBRARY_PATH
```

Verificar:
```bash
/usr/local/ompi/bin/mpicc --version
```

Nota: si no se necesita GPU-aware MPI, un OpenMPI de `apt` alcanza, pero
entonces hay que ajustar `MPI_INC_PATH`/`MPI_LIB_PATH` en `config.mk` y quitar
`-mca pml ucx` de `LAUNCH_MPIRUN_CMD`.

## 4. PETSc (con soporte CUDA + hypre)

PETSc resuelve la ecuación de Poisson cuando `HAVE_PETSC = true`. La
configuración actual usa la variante con hypre (BoomerAMG en GPU).

```bash
git clone -b release https://gitlab.com/petsc/petsc.git petsc
cd petsc

./configure \
    --prefix=/usr/local/petsc-cuda-hypre \
    --with-mpi-dir=/usr/local/ompi \
    --with-cuda=1 \
    --with-cuda-arch=89 \
    --download-hypre \
    --with-debugging=0 \
    COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3'

make PETSC_DIR=$(pwd) PETSC_ARCH=<arch-que-imprime-configure> all
sudo make PETSC_DIR=$(pwd) PETSC_ARCH=<arch-que-imprime-configure> install
```

Ajustar `--with-cuda-arch` a la arquitectura real de la GPU (89 para sm_89).
El script de `configure` imprime al final el `PETSC_ARCH` a usar en `make`.

Si no se necesita hypre (variante `petsc-cuda`, más liviana), omitir
`--download-hypre` e instalar en `/usr/local/petsc-cuda`.

Verificar que quedaron los archivos que espera `config.mk`:
```bash
ls /usr/local/petsc-cuda-hypre/include/petsc.h
ls /usr/local/petsc-cuda-hypre/lib/libpetsc.*
```

## 5. Compilar Ludwig

```bash
cd /home/bater/Sim/ludwig
cp config/unix-mpicc-default.mk config.mk   # o editar el config.mk existente
make build
```

Esto construye `src/Ludwig.exe` y `libludwig.a`. Para tests:

```bash
make test    # tests unitarios + regresión d3q19-short
make unit    # solo tests unitarios
```

Para correr los tests unitarios directamente:
```bash
LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH \
    ./tests/unit/a.out
```

Si `config.mk` apunta a rutas distintas de PETSc/OpenMPI/CUDA (por ejemplo
otra máquina), editar `PETSC_INC`, `PETSC_LIB`, `MPI_INC_PATH`,
`MPI_LIB_PATH` y `-arch=sm_XX` en consecuencia.

## 6. Entorno Python (scripts de análisis y gráficos)

Los scripts de post-proceso (`util/plot_PM_self_force.py`,
`microgel/cluster/*.py`, `docs/figures/*.py`) son independientes del build de
Ludwig — solo necesitan Python + unas pocas librerías científicas. Se
recomienda un entorno virtual dedicado, sin tocar el Python del sistema:

```bash
cd /home/bater/Sim/ludwig
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install numpy scipy matplotlib
```

Paquetes usados actualmente por los scripts del repo:
- `numpy`
- `scipy` (`scipy.interpolate.Rbf`, `scipy.spatial.ConvexHull`)
- `matplotlib` (incluye `mpl_toolkits.mplot3d`)

Activar el entorno antes de correr cualquier script:
```bash
source .venv/bin/activate
python microgel/cluster/plot.py ...
```

No existe actualmente un `requirements.txt` en el repo; si se quiere fijar
versiones, generarlo con:
```bash
pip freeze > requirements.txt
```

---

## Resumen de variables de entorno a persistir

Agregar a `~/.bashrc` (o script de activación del proyecto):

```bash
export PATH=/usr/local/ompi/bin:/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:/usr/local/cuda/lib64:$LD_LIBRARY_PATH
```

## Checklist de verificación rápida

```bash
nvidia-smi                              # driver + GPU visible
nvcc --version                          # CUDA toolkit
/usr/local/ompi/bin/mpicc --version     # OpenMPI
ls /usr/local/petsc-cuda-hypre/lib/libpetsc.*   # PETSc instalado
cd /home/bater/Sim/ludwig && make build # compila Ludwig.exe
source .venv/bin/activate && python -c "import numpy, scipy, matplotlib"
```

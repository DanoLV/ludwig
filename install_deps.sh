#!/usr/bin/env bash
#
# Instala todas las dependencias de INSTALL.md para compilar y correr Ludwig
# en esta maquina (RTX 4070, sm_89): toolchain, driver NVIDIA, CUDA Toolkit,
# OpenMPI (CUDA/UCX), PETSc (CUDA+hypre), compila Ludwig, y crea el venv
# de Python.
#
# Idempotente: se puede correr mas de una vez, salta lo que ya este hecho.
# IMPORTANTE: el paso del driver NVIDIA requiere reiniciar la maquina para
# que el modulo de kernel nuevo quede cargado. Si este script instala el
# driver, va a parar ahi y pedirte que reinicies y lo vuelvas a correr;
# la segunda corrida sigue sola con el resto.
#
# Uso:
#   chmod +x install_deps.sh
#   ./install_deps.sh
#
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PETSC_PREFIX=/usr/local/petsc-cuda-hypre
OMPI_PREFIX=/usr/local/ompi
CUDA_ARCH=89
SM_ARCH=sm_89
BUILD_DIR="${HOME}/src-ludwig-deps"

log()  { printf '\n\033[1;32m==>\033[0m %s\n' "$*"; }
warn() { printf '\n\033[1;33m!!\033[0m %s\n' "$*"; }
die()  { printf '\n\033[1;31mERROR:\033[0m %s\n' "$*" >&2; exit 1; }

mkdir -p "$BUILD_DIR"

##############################################################################
# 1. Toolchain base
##############################################################################
log "1/6 Toolchain base"
if ! dpkg -l build-essential >/dev/null 2>&1 || ! command -v pkg-config >/dev/null 2>&1; then
    sudo apt update
    sudo apt install -y build-essential gfortran git cmake pkg-config \
        python3 python3-venv python3-pip wget curl
else
    echo "ya instalado, salteando"
fi

##############################################################################
# 2. Driver NVIDIA
##############################################################################
log "2/6 Driver NVIDIA"
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
    echo "driver activo:"
    nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv
else
    if dpkg -l 2>/dev/null | grep -q "nvidia-driver-"; then
        die "El paquete del driver NVIDIA ya esta instalado pero nvidia-smi no responde.
Esto normalmente significa que falta reiniciar para cargar el modulo de kernel.
Reiniciá la máquina y volvé a correr este script: sudo reboot"
    fi

    log "instalando nvidia-driver-595-open (recomendado por ubuntu-drivers para esta GPU)"
    sudo apt update
    sudo apt install -y nvidia-driver-595-open

    warn "Driver instalado. Hace falta REINICIAR para que el modulo de kernel cargue."
    warn "Reiniciá la máquina (sudo reboot) y volvé a correr este mismo script:"
    warn "  cd '$REPO_DIR' && ./install_deps.sh"
    exit 0
fi

##############################################################################
# 3. CUDA Toolkit
##############################################################################
log "3/6 CUDA Toolkit"
if command -v nvcc >/dev/null 2>&1; then
    echo "ya instalado:"
    nvcc --version
else
    CODENAME="$(lsb_release -rs | tr -d '.')"   # ej: 26.04 -> 2604
    CANDIDATES=("ubuntu${CODENAME}" ubuntu2404 ubuntu2204)
    KEYRING_OK=""
    for C in "${CANDIDATES[@]}"; do
        URL="https://developer.download.nvidia.com/compute/cuda/repos/${C}/x86_64/cuda-keyring_1.1-1_all.deb"
        if curl -fsIL "$URL" >/dev/null 2>&1; then
            KEYRING_OK="$URL"
            echo "usando repo NVIDIA CUDA para $C"
            break
        fi
    done
    [ -n "$KEYRING_OK" ] || die "No encontre un repo de CUDA de NVIDIA compatible (probe: ${CANDIDATES[*]}). Instalar CUDA a mano con el .run desde https://developer.nvidia.com/cuda-downloads"

    curl -fsSL -o "$BUILD_DIR/cuda-keyring.deb" "$KEYRING_OK"
    sudo dpkg -i "$BUILD_DIR/cuda-keyring.deb"
    sudo apt update

    CUDA_PKG=cuda-toolkit-12-8
    if ! apt-cache show "$CUDA_PKG" >/dev/null 2>&1; then
        warn "$CUDA_PKG no esta disponible en el repo (comun en distros nuevas donde NVIDIA todavia no publico CUDA 12.8)."
        CUDA_PKG=$(apt-cache search '^cuda-toolkit-[0-9]+-[0-9]+$' \
            | grep -oE '^cuda-toolkit-[0-9]+-[0-9]+' | sort -t- -k3,3n -k4,4n | tail -1)
        [ -n "$CUDA_PKG" ] || die "no encontre ningun paquete cuda-toolkit-N-M disponible en el repo"
        warn "instalando en su lugar: $CUDA_PKG (arch sm_89 es compatible con CUDA 12 y 13)"
    fi
    sudo apt install -y "$CUDA_PKG"

    if ! grep -q '/usr/local/cuda/bin' ~/.bashrc 2>/dev/null; then
        {
            echo ''
            echo '# CUDA (agregado por install_deps.sh de Ludwig)'
            echo 'export PATH=/usr/local/cuda/bin:$PATH'
            echo 'export LD_LIBRARY_PATH=/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}'
        } >> ~/.bashrc
    fi
    export PATH=/usr/local/cuda/bin:$PATH
    export LD_LIBRARY_PATH=/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}
    nvcc --version
fi

##############################################################################
# 4. OpenMPI (con soporte CUDA/UCX)
##############################################################################
log "4/6 OpenMPI"
if [ -x "$OMPI_PREFIX/bin/mpicc" ]; then
    echo "ya instalado:"
    "$OMPI_PREFIX/bin/mpicc" --version
else
    sudo apt install -y libucx-dev ucx-utils 2>/dev/null || warn "libucx-dev/ucx-utils no disponibles via apt, sigo sin UCX explicito (openmpi puede traer su propio soporte)"

    log "resolviendo ultima version 5.0.x de Open MPI"
    # OJO: download.open-mpi.org es un bucket S3/CloudFront con listado de
    # directorio deshabilitado (403 Forbidden). La pagina human-friendly en
    # www.open-mpi.org si lista las versiones; la descarga del tarball en si
    # funciona bien contra download.open-mpi.org (solo el listado falla).
    OMPI_TARBALL=$(curl -fsSL https://www.open-mpi.org/software/ompi/v5.0/ \
        | grep -oE 'openmpi-5\.0\.[0-9]+\.tar\.gz' | sort -V -u | tail -1)
    [ -n "$OMPI_TARBALL" ] || die "no pude determinar la ultima version de OpenMPI 5.0.x"
    OMPI_VER="${OMPI_TARBALL%.tar.gz}"
    echo "version elegida: $OMPI_VER"

    cd "$BUILD_DIR"
    [ -f "$OMPI_TARBALL" ] || wget "https://download.open-mpi.org/release/open-mpi/v5.0/${OMPI_TARBALL}"
    rm -rf "$OMPI_VER"
    tar xzf "$OMPI_TARBALL"
    cd "$OMPI_VER"
    ./configure --prefix="$OMPI_PREFIX" --with-cuda=/usr/local/cuda --with-ucx
    make -j"$(nproc)"
    sudo make install

    if ! grep -q "$OMPI_PREFIX/bin" ~/.bashrc 2>/dev/null; then
        {
            echo ''
            echo '# OpenMPI (agregado por install_deps.sh de Ludwig)'
            echo "export PATH=$OMPI_PREFIX/bin:\$PATH"
            echo "export LD_LIBRARY_PATH=$OMPI_PREFIX/lib:\${LD_LIBRARY_PATH:-}"
        } >> ~/.bashrc
    fi
    export PATH="$OMPI_PREFIX/bin:$PATH"
    export LD_LIBRARY_PATH="$OMPI_PREFIX/lib:${LD_LIBRARY_PATH:-}"
fi

##############################################################################
# 5. PETSc (con soporte CUDA + hypre)
##############################################################################
log "5/6 PETSc"
if ls "$PETSC_PREFIX"/lib/libpetsc.* >/dev/null 2>&1; then
    echo "ya instalado en $PETSC_PREFIX"
else
    # /usr/local es root:root, asi que ./configure (corriendo como usuario
    # normal) no puede crear/escribir en el prefix para --download-f2cblaslapack.
    # Creamos el directorio de antemano y lo hacemos nuestro; make install
    # entonces no necesita ni sudo.
    if [ ! -d "$PETSC_PREFIX" ]; then
        sudo mkdir -p "$PETSC_PREFIX"
        sudo chown "$(id -u):$(id -g)" "$PETSC_PREFIX"
    fi

    cd "$BUILD_DIR"
    if [ ! -d petsc ]; then
        git clone -b release https://gitlab.com/petsc/petsc.git petsc
    fi
    cd petsc
    git fetch --quiet origin release 2>/dev/null || true

    PETSC_DIR="$(pwd)"
    ./configure \
        --prefix="$PETSC_PREFIX" \
        --with-mpi-dir="$OMPI_PREFIX" \
        --with-cuda=1 \
        --with-cuda-arch="$CUDA_ARCH" \
        --download-hypre=1 \
        --with-hypre-cuda=1 \
        --download-f2cblaslapack=1 \
        --with-debugging=0 \
        COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3' \
        | tee "$BUILD_DIR/petsc-configure.log"

    PETSC_ARCH=$(grep -oE 'PETSC_ARCH=[^ ]+' "$BUILD_DIR/petsc-configure.log" | tail -1 | cut -d= -f2)
    [ -n "$PETSC_ARCH" ] || die "no pude determinar PETSC_ARCH de la salida de configure; revisar $BUILD_DIR/petsc-configure.log"
    echo "PETSC_ARCH=$PETSC_ARCH"

    make PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH" all
    make PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH" install

    ls "$PETSC_PREFIX/include/petsc.h"
    ls "$PETSC_PREFIX"/lib/libpetsc.*
fi

##############################################################################
# 6. Compilar Ludwig
##############################################################################
log "6/6 Compilando Ludwig"
cd "$REPO_DIR"
[ -f config.mk ] || die "no existe config.mk en $REPO_DIR (se esperaba uno ya presente)"

export PATH=/usr/local/cuda/bin:$OMPI_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$PETSC_PREFIX/lib:$OMPI_PREFIX/lib:/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}

make build

log "entorno Python (venv)"
if [ ! -d .venv ]; then
    python3 -m venv .venv
    # shellcheck disable=SC1091
    source .venv/bin/activate
    pip install --upgrade pip
    pip install numpy scipy matplotlib
    deactivate
else
    echo "ya existe .venv, salteando"
fi

log "listo. Verificacion rapida:"
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv
nvcc --version | tail -1
"$OMPI_PREFIX/bin/mpicc" --version | head -1
ls "$PETSC_PREFIX"/lib/libpetsc.* >/dev/null && echo "PETSc OK"
ls -la "$REPO_DIR/src/Ludwig.exe"
echo "Ludwig instalado y compilado correctamente."

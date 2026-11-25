# Scripts de Ludwig: qué hace cada uno y cómo se usa

Generado por `util/documentar_scripts.py` el 2026-11-17 09:47 a partir de la documentación que trae cada archivo (docstring o comentarios iniciales). No editar a mano: documentar el script en su propio archivo y volver a correr el generador.

Alcance: `util/`, `code_changes/`, scripts de la raíz y `<proyecto>/scripts/` (scripts de cada proyecto; en `local/` quedan enlaces a ellos, ver `<proyecto>/scripts/README_ubicacion.md`). No incluye las copias de scripts dentro de cada corrida.

| Sección | Scripts | Sin documentación |
|---|---|---|
| Herramientas generales (Python / shell) | 10 | 1 |
| Programas en C de util/ (se compilan con util/Makefile) | 19 | 0 |
| Proyecto capacitor (capacitor/scripts/) | 76 | 0 |
| Proyecto microgel (microgel/scripts/) | 170 | 30 |

## Herramientas generales (Python / shell)


### `code_changes/`

#### `code_changes/extraer_cambios.py`
> Genera code_changes/CAMBIOS_CODIGO.md a partir de las marcas de cambio local
> del código de Ludwig (convención del proyecto, ver CLAUDE.md):
>
>     /*CHANGE INIT - <AAAAMMDD> <descripción> */  ...  /*CHANGE END - ... */
>
> Para cada marca: fecha, descripción, archivo y rango de líneas (INIT..END).
> Agrupa por fecha (más reciente primero) y, dentro de cada fecha, por descripción.
> Las entradas escritas a mano en code_changes/NOTAS_MANUALES.md (si existe) se
> copian al principio del documento: sirven para explicar el porqué de un cambio
> con más detalle del que cabe en la marca.
>
> Uso (desde la raíz de Ludwig o desde cualquier lugar):
> … (2 líneas más en el archivo)

Uso (del archivo):
    python3 code_changes/extraer_cambios.py


### `./`

#### `install_deps.sh`
> Instala todas las dependencias de INSTALL.md para compilar y correr Ludwig
> en esta maquina (RTX 4070, sm_89): toolchain, driver NVIDIA, CUDA Toolkit,
> OpenMPI (CUDA/UCX), PETSc (CUDA+hypre), compila Ludwig, y crea el venv
> de Python.
> Idempotente: se puede correr mas de una vez, salta lo que ya este hecho.
> IMPORTANTE: el paso del driver NVIDIA requiere reiniciar la maquina para
> que el modulo de kernel nuevo quede cargado. Si este script instala el
> driver, va a parar ahi y pedirte que reinicies y lo vuelvas a correr;
> la segunda corrida sigue sola con el resto.
> Uso:
> chmod +x install_deps.sh
> ./install_deps.sh


### `util/`

#### `util/archive_runs.py`
> Para cada proyecto de la forma <proyecto>/local/ en la raiz de Ludwig (p.ej.
> capacitor/, microgel/), crea/actualiza <proyecto>/runs/ como un espejo
> liviano de <proyecto>/local/: copia la estructura completa de carpetas, pero
> excluyendo los volcados de datos pesados de Ludwig (dist-/psi-/qsi-/vel-/
> rho-/efield-/config.cds-<paso>, sus "-metadata", binarios, logs) y cualquier
> archivo >2MB que no sea .png. Al espejar el arbol completo (no solo las
> carpetas con "input"), las carpetas compartidas de graficos/scripts que
> viven al lado de varias corridas (como en capacitor/local/WALL-DIRICHLET/
> graficos/) se copian solas, en su lugar natural.
>
> Copia liviana LOCAL de las corridas: lo minimo necesario para volver a
> correr cada simulacion (input, .petscrc, condiciones iniciales chicas,
> … (13 líneas más en el archivo)

Uso (del archivo):
    Uso (desde cualquier directorio):
    python3 util/archive_runs.py --dry-run   # solo reporta, no copia nada
    python3 util/archive_runs.py             # copia de verdad
    python3 util/archive_runs.py microgel    # solo ese proyecto (atajo)

#### `util/documentar_scripts.py`
> Genera SCRIPTS.md en la raíz de Ludwig: qué hace y cómo se usa cada script
> (2026-11-13, pedido del usuario).
>
> Fuente: la documentación que trae cada archivo (no se ejecuta nada):
>   * Python: docstring del módulo (o comentarios iniciales); "Uso"/"Usage" y argumentos
>     de argparse si los hay; sys.argv[...] usados.
>   * Shell: bloque de comentarios inicial (tras el shebang).
>   * C (util/*.c): bloque de comentario inicial.
> Alcance:
>   * herramientas generales: util/, code_changes/, scripts de la raíz;
>   * scripts de cada proyecto: <proyecto>/scripts/ (desde 2026-11-13; en local/ quedan
>     enlaces a estos archivos).
> … (3 líneas más en el archivo)

Uso (del archivo):
    Uso:  python3 util/documentar_scripts.py

#### `util/extract.sh`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

#### `util/guardar_bases.py`
> Guarda en <proyecto>/saved_data/ lo necesario para rehacer cada corrida (2026-11-13,
> pedido del usuario), sin guardar un input ni un archivo de coloides por corrida:
>
>   saved_data/<carpeta-serie>/input.base          input real (con comentarios) de la corrida
>                                                  que menos difiere del resto de la carpeta
>   saved_data/<carpeta-serie>/coloides_<k>.cds    archivo(s) de coloides base de la carpeta
>   saved_data/<carpeta-serie>/deltas.json         por corrida: claves del input que cambian /
>                                                  se agregan / se quitan respecto de input.base,
>                                                  base de coloides que usa y líneas que cambian
>
> Coloides: dentro de una carpeta-serie el contenido más frecuente es la primera base; otro
> contenido con el mismo número de líneas se guarda como cambios de línea si difiere en
> … (8 líneas más en el archivo)

Uso (del archivo):
    Uso:  python3 util/guardar_bases.py [capacitor] [microgel]

#### `util/plot_PM_self_force.py`
> plot_PM_self_force.py
>
> Plots PM self-force components and magnitude from PM_self_force_diag.csv.
>
> Usage examples
> --------------
> # 1D: fix dx=0, dy=0, vary dz
> python plot_PM_self_force.py --fix dx=0.0 dy=0.0 --vary dz
>
> # 1D: fix dx=0.2, dz=0.4, vary dy
> python plot_PM_self_force.py --fix dx=0.2 dz=0.4 --vary dy
>
> … (5 líneas más en el archivo)

Uso (del archivo):
    python plot_PM_self_force.py --fix dx=0.0 --vary dy dz --plot3d Fmod
    python plot_PM_self_force.py --fix dz=0.0 --vary dx dy

Argumentos (argparse): `--csv`; `--fix`: Fixed axes, e.g. --fix dx=0.0 dy=0.2; `--vary`: Varying axis/axes: one name for 1D, two for 3D; `--plot3d`: Single column for 3D surface (omit to plot all four); `--tol`; `--no-show`: Save PNG without opening a window

#### `util/registro_propositos.py`
> Idea / propósito de cada carpeta-serie, para util/registro_simulaciones.py.
> (2026-11-13). Editable a mano: es la fuente de la columna "Descripción".
>
> Cada entrada: (proyecto, regex sobre el nombre del grupo, idea, fuente).
> Se usa la PRIMERA que coincida (poner las más específicas antes).
> Fuentes:
>   README          = README de la carpeta (o README de WALL-DIRICHLET, tabla "Contenido")
>   cabecera input  = comentarios al principio del input de las corridas
>   sesiones        = conversaciones previas con Claude donde se armaron las corridas
>   INFERIDO        = deducido del nombre de la carpeta y de los parámetros: VERIFICAR

#### `util/registro_simulaciones.py`
> Registro de simulaciones por proyecto (2026-11-13, pedido del usuario).
>
> Genera <proyecto>/Registro_simulaciones.xlsx para capacitor/ y microgel/:
>   * una pestaña por SERIE (agrupación científica de carpetas, ver SERIES_*),
>   * dentro de cada pestaña, una fila-título por carpeta-serie y una fila por
>     simulación (= cada carpeta que contiene un archivo "input"),
>   * columnas: ruta, descripción, fecha, ubicación (local / runs / borrada),
>     ID del archivo de coloides y los parámetros del input que varían en la
>     pestaña; los que no varían se listan en una fila "Comunes" arriba,
>   * columna "Scripts" (2026-11-16): scripts de crear/correr/analizar asociados a la
>     corrida por una regla heurística (ver run_scripts y la hoja Leyenda),
>   * hoja "Coloides": un ID por contenido distinto de config.cds.init.001-001,
> … (10 líneas más en el archivo)

Uso (del archivo):
    python3 util/registro_simulaciones.py            # ambos proyectos
    python3 util/registro_simulaciones.py microgel   # uno solo

#### `util/rehacer_corrida.py`
> Reconstruye el input y el config.cds.init.001-001 de una corrida a partir de
> <proyecto>/saved_data/<carpeta-serie>/ (input.base, coloides_*.cds, deltas.json),
> generados por util/guardar_bases.py (2026-11-13).
>
>   python3 util/rehacer_corrida.py <proyecto> <ruta de la corrida> [carpeta destino]
>
> <ruta de la corrida> es la columna "Ruta" del Registro_simulaciones.xlsx
> (p. ej. SLIDES-LD8-hann08-corr_yes/pos_12.00_24.00_24.00). Sin carpeta destino solo
> muestra los cambios respecto de la base. Con destino escribe ahí "input" y
> "config.cds.init.001-001" (no sobrescribe archivos existentes).
>
> El input reconstruido tiene las mismas claves y valores que el original; los
> … (4 líneas más en el archivo)

Lee argumentos de línea de comandos: sys.argv[1, 2, 3]


## Programas en C de util/ (se compilan con util/Makefile)


### `util/`

#### `util/PM_self_force.c`
> PM_self_force.c
> Compute the particle-mesh (PM) self-force on a single charged subgrid
> particle using the Peskin4 spread/gather kernel and an FFT Poisson solver.
> For each position on a grid covering one voxel (between lattice nodes),
> the program:
> 1. Places a particle with charge q at position r.
> 2. Spreads charge to the mesh with the Peskin4 kernel.
> 3. Solves Poisson's equation via FFT.
> 4. Computes the electric field E = -grad(psi).
> 5. Gathers E back to the particle with the same Peskin4 kernel.
> 6. Computes the self-force F = q * E_gathered.
> Output: tab-separated table to stdout:
> … (18 líneas más en el archivo)

#### `util/bath.c`
> bath.c
> Compile and link with -lm for the maths library.
> This utility produces an output file suitable for initialising
> a capillary structure within a bath in Ludwig.
> See also capillary.c
> $Id: bath.c,v 1.2 2010-08-11 08:21:39 kevin Exp $
> Edinburgh Soft Matter and Statistical Physics Group and
> Edinburgh Parallel Computing Centre
> Kevin Stratford (kevin@epcc.ed.ac.uk)
> (c) 2008 The University of Edinburgh

#### `util/capillary.c`
> capillary.c
> This utility produces an output file suitable for initialising
> a capillary structure in Ludwig.
> Some examples of uniform, and less uniform capillary initialisations.
> Options are selected at compile time (below) at the moment.
> The output should be capillary.dat      [for human consumption]
> capillary.001-001  [for initial input to run]
> Edinburgh Soft Matter and Statistcal Physics Group and
> Edinburgh Parallel Computing Centre
> (c) 2008-2024 The University of Edinburgh
> Contributing authors:
> Kevin Stratford (kevin@epcc.ed.ac.uk)

#### `util/coll_squ_subgrid_init.c`
> coll_squ_subgrid_init.c
> Produce a file of colloid/squirmer/subgrid particle information suitable
> for reading into the main code.
> Particle positions are initialised at random to provide a requested
> volume fraction. Flat boundary walls may be included. A 'grace'
> distance dh may be specified to prevent the initial particle
> positions being too close together, or too close to the wall.
> For compilation instructions see the Makefile.
> $ make
> $ ./coll_squ_subgrid_init
> should produce a file config.cds.init.001-001 in the specified format.
> Edinburgh Soft Matter and Statistical Physics Group and
> … (4 líneas más en el archivo)

#### `util/colloid_file.c`
> colloid_file.c
> Produce a file of colloid information suitable for reading into
> the main code in serial or in parallel.
> First decide how many colloids are required and with what
> properties. The format of the output should be one 4-byte
> integer which is the number of colloids in the file,
> followed by this number of colloid_state_t structs. This
> is set out below.
> The output is always in a single file, and can be ascii or
> binary. Use the standard functions in colloid.c to do the
> I/O as below.
> Compile with, e.g., cc -I../src ../src/colloid.c colloid_file.c
> … (7 líneas más en el archivo)

#### `util/colloid_init.c`
> colloid_init.c
> Produce a file of colloid information suitable for reading into
> the main code.
> Colloid positions are initialised at random to provide a requested
> volume fraction. Flat boundary walls may be included. A 'grace'
> distance dh may be specified to prevent the initial colloid
> positions being too close together, or too close to the wall.
> For compilation instructions see the Makefile.
> $ make colloid_init
> $ ./colloid_init [-a a0] [-h ah] [-v volume_fraction]
> should produce a file config.cds.init.001-001 in the specified format.
> Options:
> … (7 líneas más en el archivo)

#### `util/extract.c`
> extract.c
> This is under revision as part of the move to new metadata.
> For most users, there will be no need to consolidate
> parallel output. Only post-processing is relevant (unrolling
> Lees Edwards plane, diagonalisation of liquid crystal order parameter).
> ./extract data-file
> The relevant metadata file is located from the file stub.
> Most options are still available (see below).
> Older version ...
> This program deals with the consolidation of parallel I/O
> (and in serial with Lees-Edwards planes). This program is
> serial.
> … (31 líneas más en el archivo)

#### `util/extract_colloids.c`
> extract_colloids.c
> Convert an output file to a csv file suitable for Paraview.
> The csv file uses three extra particles at (xmax, 0, 0)
> (0, ymax, 0) and (0, 0, zmax) to define the extent of the
> system.
> If you want a different set of colloid properties, you need
> to arrange the header, and the output appropriately. The read
> can be ascii or binary and is set by the switch below.
> CHANGE -------------------
> Add bonds to the csv file. Nbonds and ids.
> This adds the info needed to calculate bond lengths later
> CHANGE -------------------
> … (12 líneas más en el archivo)

#### `util/extract_squirmer_polymer.c`
> extract_colloids.c
> Convert an output file to a csv file suitable for Paraview.
> The csv file uses three extra particles at (xmax, 0, 0)
> (0, ymax, 0) and (0, 0, zmax) to define the extent of the
> system.
> If you want a different set of colloid properties, you need
> to arrange the header, and the output appropriately. The read
> can be ascii or binary and is set by the switch below.
> For compilation instructions see the Makefile.
> $ make extract_colloids
> $ ./a.out <colloid file name stub> <nfile> <csv file name>
> where the
> … (10 líneas más en el archivo)

#### `util/initialise.c`
> initialise.c
> Toolkit to create parallel input files of the
> tensor order parameter phi and LB-distributions in
> row-major order.
> The details of the I/0-topology have to be given below,
> whereas the individual configuration of order parameter
> and the distribution are set are hard-coded in separtate
> routines.
> ./initialise
> WARNING: Backup of the input file is required as it
> might be replaced when the command is issued.
> Compile with $(CC) -o initialise initialise.c -lm
> … (2 líneas más en el archivo)

#### `util/length_from_sk.c`
> length_from_sk.c
> This program computes a domain length scale L from the spherically
> averaged structure factor S(k). The latter is, in turn, the
> Fourier transform of the order parameter distribution phi(r).
> This is intended for use with binary fluid calculations
> on a cubic lattice size NXxNYxNZ, which may be reduced
> to size NXRxNYRxNZR for the purposes of the Fourier
> tranform.
> This program requires FFTW3 for the Fourier transforms:
> details will vary locally.
> Quadrature proceeds by binning the structure factor in k-space
> on a regular interval.
> … (16 líneas más en el archivo)

#### `util/microgel_poly_cross_density.c`
> microgel_poly_cross_density.c
> Produces a file containing microgel polymer information suitable for input
> into the main code. The microgel is a network of monomers distributed within
> a spherical volume. Monomer positions are initialized using a Fibonacci
> approach to distribute them evenly over the sphere. Then, linear chains are
> grown inward from the surface monomers, with an average chain length adjusted
> to match the target density. Afterward, crosslinks are created between chains
> according to the specified crosslinking density.
> For compilation instructions see the Makefile.
> $ make microgel_poly_cross_density
> $ ./microgel_poly_cross_density
> should produce a file config.cds.init.001-001 in the specified format.
> … (10 líneas más en el archivo)

#### `util/multi_poly_init.c`
> multi_poly_init.c
> Produce a file of polymers information suitable for reading into
> the main code.
> Polymer positions are initialised at random. Multiple polymers can
> be generated. This code is suitable for dilute or intermediate polymer
> solutions.
> It should not be used if the system is dense. Boundary walls can not
> be included.
> A 'grace' distance dh may be specified to prevent the initial monomer
> positions being too close together.
> For compilation instructions see the Makefile.
> $ make multi_poly_init
> … (9 líneas más en el archivo)

#### `util/multi_squ_poly_init.c`
> multi_squ_poly_init.c
> Produce a file of squirmer and polymers information suitable for reading into
> the main code.
> Polymer/squirmer positions are initialised at random. Multiple polymers and
> squirmers can be generated simultaneously. This code is suitable for dilute
> or intermediate polymer solution. It should not be used if the system is dense.
> Boundary walls can not be included.
> A 'grace' distance dh may be specified to prevent the initial particle
> positions being too close together.
> For compilation instructions see the Makefile.
> $ make multi_squ_poly_init
> $ ./a.out
> … (7 líneas más en el archivo)

#### `util/noise_phi_sq.c`
> noise_phi_sk.c
> This program requires FFTW3 for the Fourier transforms:
> details will vary locally.
> Quadrature proceeds by binning the structure factor in k-space
> on a regular interval.
> It should be compiled with something like
> $(CC) -I $(FFTW3)/include length.c -L$(FFTW3)/lib -lfftw3 -lm
> where CC is the C compiler and FFTW3 is the location of the
> FFTW3 library
> A single file, name given on the command line, is assumed to
> provide the order parameter data. See the routine read_phi()
> below for details.
> … (8 líneas más en el archivo)

#### `util/polarizer.c`
> polarizer.c
> This program simulates the appearance of a liquid-crystalline sample
> under polarised light. Director and scalar order parameter data in
> vtk-format is required. The output consists of a file 'polar-...'
> that contains the light intensity in vtk-format.
> The file names should be (in the current directory) e.g.,:
> ./polarizer lcd-00100000.vtk lcs-00100000.vtk x|y|z
> Cartesian_direction is x, y, or z and defines the direction of the
> incident light.
> COMMAND LINE INPUT
> (1) the director filename (must be "lcd-*.vtk")
> (2) the scalar order parameter filename
> … (13 líneas más en el archivo)

#### `util/porous_solid_cube.c`
> porous_solid_cube.c
> Produces a file containing a porous solid network suitable for input
> into the main code, using the same colloid configuration file format
> as microgel_poly_cross_density.c.
> The cube side length, the mean bond distance, its tolerance ("delta")
> and the point density are specified on the command line. The cubic
> simulation volume is then filled with points at random subject to the
> constraint that no two points may lie closer than (lbond - delta).
> Bonds are then created between any pair of points whose separation
> falls within [lbond - delta, lbond + delta], up to a maximum number
> of bonds per point.
> Periodicity: each of the x, y, z directions can independently be
> … (45 líneas más en el archivo)

#### `util/process.c`
> process.c
> For a porous media data set (e.g., X-ray tomography) stored as
> ASCII (0 = fluid, 1 = solid), we generate a file suitable for
> input to Ludwig.
> There are two main concerns:
> 1) the structure must be periodic in the flow direction, meaning
> a 'reflection' is required. The flow direction is assumed to
> be the x-direction.
> 2) Similarly, the y-direction and z-direction must be periodic.
> To satify this, we pad the edges with solid sites.
> The input file is assumed to have one data value per line, with
> the z-direction running fastest.
> … (13 líneas más en el archivo)

Uso (del archivo):
    Usage:

#### `util/psi_init.c`
> psi_init.c
> Initialise the electrokinetic quantities and write to file
> suitable for reading at time t = 0.
> The user must ensure that the initial charges are conistent
> with the values for unit charge, valency, and so on, set in
> the input file at run time.
> Overall electroneutrality (solid plus fluid) must also be
> ensured by the user.
> Compilation:
> Compile the main code for serial execuation.
> Compile this file via
> $(CC) -I../mpi_s -I../targetDP psi_init.c ../src/libludwig.a \
> … (7 líneas más en el archivo)


## Proyecto capacitor (capacitor/scripts/)


### `capacitor/scripts/H-microgel-ruido/`

#### `capacitor/scripts/H-microgel-ruido/analizar_largo.py`
> Análisis de las corridas largas H-microgel-ruido/{neutro,superficie}-largo (2026-11-13).
>
> Posiciones cada 100 pasos (no hay pasos consecutivos: el desplazamiento por paso de
> analizar_ruido.py no se puede medir aquí). Para cada corrida imprime:
>   * serie por bloques de 5000 pasos (heurístico: 50 volcados por bloque) de
>       Rg (radio de giro, sin pesos), <U_enlaces>/kT (ref. equipartición armónica 35/2 = 17.5),
>       z del centro de masa, d mínimo a los planos z = 1.5 / 33.5 (h = d - 3.5);
>     sirve para decidir a ojo desde cuándo la corrida está estacionaria: NO elige la ventana.
>   * fluido: <u_a^2>/(kT/rho) en cada vel-* (equipartición LB, Dünweg & Ladd 2009).
>   * máximo desplazamiento de un monómero entre volcados (100 pasos), solo como control de saltos.
> Reutiliza read_cds y las constantes de analizar_ruido.py. Solo lee archivos.
>
> … (1 líneas más en el archivo)

Uso (del archivo):
    Uso:  python3 capacitor/scripts/H-microgel-ruido/analizar_largo.py [carpetas...]

#### `capacitor/scripts/H-microgel-ruido/analizar_pared.py`
> Análisis de la prueba de pared con ruido (H-microgel-ruido/pared-*, 2026-11-13).
>
> Para cada corrida (un monómero libre que empieza a d = 4.0 del plano z = 1.5):
>   * si terminó o abortó (mensaje de Ludwig);
>   * h = d - 3.5 (d = distancia del centro al plano sólido/fluido más cercano; 3.5 = radio del
>     kernel Hann-8: h < 0 significa que el kernel toca nodos de pared): mínimo, fracción h < 0;
>   * desplazamiento por paso: rms por componente y máximo frente a drmax = 0.8;
>   * V(h)/kT de la Steele (Hoang Ngoc Minh et al. 2023, ecs. 7-8, cortada en z* = 0.98698 sigma_w):
>     máximo alcanzado y, dentro del alcance de la pared (0 < h < z*), comparación del promedio
>     <V>/kT de la corrida con el de Boltzmann, int V e^{-V/kT} dh / int e^{-V/kT} dh en (0, z*).
>     Esa comparación supone equilibrio local junto a la pared (densidad uniforme sin pared);
>     con 2000 pasos el muestreo es corto: es un indicador, no un test estadístico.
> … (2 líneas más en el archivo)

Uso (del archivo):
    Uso:  python3 capacitor/scripts/H-microgel-ruido/analizar_pared.py [carpetas...]

#### `capacitor/scripts/H-microgel-ruido/analizar_ruido.py`
> Análisis de las pruebas cortas H-microgel-ruido (2026-11-13): ¿el ruido térmico es demasiado fuerte?
>
> Para cada corrida (argumentos = carpetas; por defecto las 4 de H-microgel-ruido):
>  1. Desplazamiento por paso de cada monómero (posiciones de config.cds en cada paso; imagen
>     mínima en x, y; z no periódico). rms por componente y máximo |dr| frente a drmax = 0.8
>     (límite de Ludwig). Referencias: sqrt(2 D_desnudo) con D_desnudo = kT/(6 pi eta)(1/a_h - 1/al)
>     (empujón aleatorio propio de la partícula) y sqrt(2 kT/(6 pi eta a_h)) (difusión total de Stokes).
>  2. Equipartición en los enlaces: U = sum_enlaces k/2 (r - r0)^2. Para fluctuaciones chicas la red
>     es armónica y <U> = (rango de la matriz de rigidez)/2 kT = 35/2 kT (rango 35 según el README de
>     G-microgel-relax). Teorema de equipartición; la aproximación armónica requiere sqrt(kT/k) << r0.
>  3. Fluido: <u_a^2> = kT/rho por componente (equipartición de las fluctuaciones del LB; Dünweg &
>     Ladd, Adv. Polym. Sci. 221, 89 (2009)), sobre nodos de fluido z = 2..33 de los vel-*.
> … (5 líneas más en el archivo)

Uso (del archivo):
    Uso:  python3 capacitor/scripts/H-microgel-ruido/analizar_ruido.py [carpetas...]

#### `capacitor/scripts/H-microgel-ruido/correr.sh`
> H-microgel-ruido: corre las 4 pruebas cortas en paralelo (reanudable: salta las terminadas).
> Uso: capacitor/local/WALL-DIRICHLET/H-microgel-ruido/correr.sh

#### `capacitor/scripts/H-microgel-ruido/correr_largo.sh`
> Corrida larga H-microgel-ruido: neutro-largo y superficie-largo en paralelo (reanudable).

#### `capacitor/scripts/H-microgel-ruido/correr_pared.sh`
> Prueba de pared (4 corridas de 1 monómero) en paralelo.

#### `capacitor/scripts/H-microgel-ruido/crear.py`
> H-microgel-ruido (2026-11-13, pedido del usuario): microgel de G-microgel-relax entre placas
> a psi = 0, ahora CON ruido térmico (LB fluctuante + empujón aleatorio de las partículas
> subgrid), en la escala de unidades de Ingo (escala "C"):
>
>   kT = 1e-4, eta = 0.01, rho = 1 (como INGO-1PART-LBTERM), a_h = 0.3078881 (el de Ingo;
>   al = 0.4743549 se deja como en G).  Para conservar el problema físico de G se re-escalan
>   con kT:  epsilon = 1/(4 pi lB kT) con lB = 0.795775 (igual que G) -> 1000;
>            bond_harmonic_k = 10 kT / LU^2 = 1e-3 (G: 1e-4 con kT 1e-5);
>            soft_sphere_epsilon = kT = 1e-4 (G: 1e-5 = kT).
>   Cocientes que se conservan: lB/r0 = 0.199, sqrt(kT/k)/r0 = 7.9 %, eps_rep/kT = 1, caja.
>
> Prueba corta: 5000 pasos, posiciones de los coloides en CADA paso (para el desplazamiento
> … (6 líneas más en el archivo)

Uso (del archivo):
    Uso:  python3 capacitor/scripts/H-microgel-ruido/crear.py

#### `capacitor/scripts/H-microgel-ruido/crear_largo.py`
> Corrida larga H-microgel-ruido (2026-11-13, pedido del usuario): neutro y superficie, pared
> Steele eps_w = 1e-5 (la que resultó estable con ruido en la prueba de pared), 50000 pasos,
> posiciones cada 100 pasos, volcados completos cada 10000 (permiten reiniciar / extender).
> Copia de <gel>-steele1e-5 (misma configuración inicial, misma semilla) cambiando solo la salida
> y la duración.
>
> Uso:  python3 capacitor/scripts/H-microgel-ruido/crear_largo.py
> Crea H-microgel-ruido/<gel>-largo/ (no sobrescribe).

#### `capacitor/scripts/H-microgel-ruido/crear_pared.py`
> Prueba de pared con ruido térmico (H-microgel-ruido, 2026-11-13, pedido del usuario).
>
> Un monómero libre, empezando a d = 4.0 del plano sólido/fluido inferior (z = 5.5; x = y = 16.5),
> mismas unidades y parámetros que las corridas H (escala de Ingo: kT 1e-4, eta 0.01, rho 1,
> a_h 0.308, al 0.474; epsilon 1e3), sin sal, placas a psi = 0. 2000 pasos, posición en cada paso.
> Cuatro corridas: carga q = 0 y q = +1 (la imagen lo atrae a la pared) x las dos paredes Steele
> (eps_w = 6.2572e-19 "steeleG" y 1e-5 "steele1e-5", sigma_w = 1).
> Base: el input de H-microgel-ruido/superficie-<pared> y el bloque del monómero 1 de su
> config.cds.init (q = 1, a_h ya cambiado), sin enlaces.
>
> Uso:  python3 capacitor/scripts/H-microgel-ruido/crear_pared.py [z0]   (z0 = 5.5 por defecto)
> Nota: con eps_w = 1e-5, z0 = 5.5 (h = 0.5) cae dentro de la zona repulsiva y aborta en el paso 1:
> … (2 líneas más en el archivo)

Lee argumentos de línea de comandos: sys.argv[1]


### `capacitor/scripts/SOLIDO-L48-lb3.5-relax/`

#### `capacitor/scripts/SOLIDO-L48-lb3.5-relax/convergencia.py`
> Seguimiento de la relajacion: extrae de run.log, cada freq_statistics pasos,
> la longitud minima/maxima de enlace, la energia armonica y la velocidad maxima
> de los coloides. Uso: python3 convergencia.py [run.log]  (guarda convergencia.png/.csv)

Lee argumentos de línea de comandos: sys.argv[1]

#### `capacitor/scripts/SOLIDO-L48-lb3.5-relax/exportar_csv.py`
> Exporta una configuracion de coloides de Ludwig (config.cds*.001-001) a ParaView.
> Uso: python3 exportar_csv.py [config] [prefijo_salida] [L]
> Genera <prefijo>_puntos.csv (id,x,y,z,nbonds,isfixed), <prefijo>_enlaces.csv
> (id1,id2,longitud con imagen minima,periodico) y <prefijo>.vtk (puntos + enlaces
> no periodicos como lineas; los que cruzan la caja se omiten para no dibujar
> lineas de lado a lado).

Lee argumentos de línea de comandos: sys.argv[1, 2, 3]

#### `capacitor/scripts/SOLIDO-L48-lb3.5-relax/run.sh`
> Lanza la relajacion en segundo plano:  ./run.sh
> Seguimiento:  tail -f run.log   |   python3 convergencia.py


### `capacitor/scripts/TEST-NOFLUX-ADV/`

#### `capacitor/scripts/TEST-NOFLUX-ADV/analizar.py`
> TEST-NOFLUX-ADV (2026-11-16): conservación de los iones con advección centred + LB térmico + paredes.
>
> Para cada corrida <ver>/<caso> y cada volcado qsi (cada 100 pasos):
>   N_a(t) = sum de n_a sobre los nodos FLUIDOS (con pared: z = 2..L-1; sin pared: todos), a = +, -
>   dN_a(t)/N_a(0) = cambio relativo respecto del paso 1 (la masa de iones de NP debería conservarse a redondeo)
> También: n en los nodos de pared (deben quedar en el valor inicial) y RMS de u_z en la primera capa de fluido
> (velocidad normal a la pared, la que según la lectura del código produce la pérdida).
> Al final, si existen antes/ y despues/ del caso sinpared, compara qsi byte a byte (control de regresión).
> Orden de los datos: x lento, z rápido (como los demás scripts del proyecto).
> Solo lee archivos. Uso: python3 capacitor/scripts/TEST-NOFLUX-ADV/analizar.py

#### `capacitor/scripts/TEST-NOFLUX-ADV/correr.sh`
> TEST-NOFLUX-ADV (2026-11-16): corre <ver>/{pared,sinpared} en serie (corridas cortas, 16^3).
> Uso: bash capacitor/scripts/TEST-NOFLUX-ADV/correr.sh <ver> [caso ...]

#### `capacitor/scripts/TEST-NOFLUX-ADV/crear.py`
> TEST-NOFLUX-ADV (2026-11-16, pedido del usuario): conservación de los iones con advección centred de NP,
> LB térmico, paredes y electrolito, ANTES y DESPUÉS de agregar la condición de no-flujo (np_no_flux_boundary)
> al driver de GPU de Nernst-Planck.
>
> Motivo: en el driver de GPU, np_advective_fluxes_kernel no mira el mapa; un enlace fluido-pared lleva el
> flujo W_p (u_i.c_p)/2 (n_i + n_pared)/2, que sale del nodo fluido y no entra en ninguno (la pared no se
> actualiza). Se espera que el total de iones en el fluido NO se conserve antes y SÍ después.
>
> Corridas (cada una con su Ludwig.exe copiado):
>   <ver>/pared     caja 16x16x18 con paredes sólidas en z = 1 y z = 18 (porous_media_init wall_z), placas a psi = 0
>   <ver>/sinpared  la misma sin paredes (periódica): control de regresión, antes y después deben ser idénticas
> con <ver> = antes (exe 221f3229) o despues (exe nuevo). Sin partículas.
> … (3 líneas más en el archivo)

Uso (del archivo):
    Uso:  python3 capacitor/scripts/TEST-NOFLUX-ADV/crear.py <ver> <exe> [caso ...]

Lee argumentos de línea de comandos: sys.argv[1, 2]


### `capacitor/scripts/WALL-DIRICHLET/`

#### `capacitor/scripts/WALL-DIRICHLET/analizar_B.py`
> Analiza la prueba B (capacitor con electrolito, lambda_D = 4, estacionario):
>   1) residuo de Poisson nodo a nodo con la regla de pared (ec. 15)
>   2) equilibrio de Boltzmann de los iones: ln(rho+/rho-) = -2 psi
>   3) perfil vs Debye-Hückel discreta (misma red) y continua (ec. 17)
>   4) campo en el nodo junto a la pared vs derivada de DH continua
> Uso: analizar_B.py [directorio de B] [paso]

Lee argumentos de línea de comandos: sys.argv[1, 2]

#### `capacitor/scripts/WALL-DIRICHLET/analizar_C.py`
> Prueba C: una particula subgrid cargada (q = +1, Hann-8) entre dos paredes a
> tierra, sin sal, en funcion de su distancia a la pared inferior.
>
> Geometria: caja 32 x 32 x 34, paredes MAP_BOUNDARY en z = 1 y z = 34, planos
> medios en z = 1.5 y 33.5 -> L = 32; d = z0 - 1.5. Periodica en x, y (A = 32^2).
>
> Teoria (unidades de Ludwig: psi reducido, 4 pi l_B = beta/epsilon = 10):
>
>  * Carga inducida (reciprocidad de Green, Jackson Prob. 1.13): una carga Q a distancia d de la placa
>    inferior induce -Q (L-d)/L en la inferior y -Q d/L en la superior. Se suma
>    la particula y el fondo de contraiones que agrega psi_electroneutral
>    (uniforme en los nodos fluidos, salvo n_near nodos pegados a la particula,
> … (15 líneas más en el archivo)

Uso (del archivo):
    Uso: analizar_C.py [directorio de la bateria] [orden Hann n]

Lee argumentos de línea de comandos: sys.argv[1, 2]

#### `capacitor/scripts/WALL-DIRICHLET/analizar_linealidad.py`
> Test de linealidad (Hann-8, L = 32), q = 0.25, 1, 2 en la misma configuracion.
> La fuerza medida es F = q E, y el campo de imagen es proporcional a q, asi que F/q^2 y Q_inf/q son
> independientes de q si Debye-Huckel linealizado es exacto; la teoria lineal da un unico valor de F/q^2.

#### `capacitor/scripts/WALL-DIRICHLET/cola_hannp4.sh`
> Cola Hann cos^4: espera a las corridas chicas (dx = 2, 1) y corre las dos dx = 1/2 de a una, en paralelo con la cola
> del Peskin-6 (dos corridas grandes juntas entran en los 8 GB; cinco arrancando juntas no).

#### `capacitor/scripts/WALL-DIRICHLET/cola_peskin6.sh`
> Cola Peskin-6: espera a las corridas chicas (dx = 2, 1) y corre las dos dx = 1/2 de a una (8 GB de GPU no alcanzan
> para ambas junto con otras).

#### `capacitor/scripts/WALL-DIRICHLET/comparar_C.py`
> Figura comparativa de las baterias C: error relativo del campo de imagen en
> funcion de la distancia a la pared, con/sin fondo de contraiones, con/sin
> correccion pm_sr, Hann-8 y Hann-6. Escribe graficos/C_comparacion.png.

#### `capacitor/scripts/WALL-DIRICHLET/comparar_pared_paper.py`
> Wall potential of Ludwig (wall_ss_cut, soft sphere with energy+force shift) vs the ion-wall potential of
> Hoang Ngoc Minh et al., J. Chem. Phys. 158, 104103 (2023), Eqs. (6)-(8):
>    V_ads(z) = 2 pi eps_w [ (2/5)(s/z)^10 - (s/z)^4 - sqrt(2) / (3 (z/s + 0.61/sqrt(2))^3) ]     (10-4-3 Steele, fcc (100) face, lattice parameter s = sigma_w)
>    V_rep(z) = V_ads(z) - V_ads(z*)  for z <= z* = 0.987 s (minimum of V_ads);  0 for z > z*       (purely repulsive wall)
> Ludwig: v(h) = eps (s/h)^nu - eps (s/hc)^nu + (eps nu/s)(s/hc)^(nu+1) (h - hc), h < hc;  eps/kT = 1.3105e-16, s = 1, nu = 12, hc = 3.
> Usage: python3 comparar_pared_paper.py -> ../graficos/comparacion_potencial_pared_paper.png

#### `capacitor/scripts/WALL-DIRICHLET/convergencia_campo.py`
> Orden de convergencia del campo en el nodo junto a la pared, sobre la solucion
> de Debye-Hückel discreta (la misma que Ludwig reproduce a 1e-8), con kappa*L=8.
> Compara: diferencia central actual, factor 2 en diferencias (implementado) e
> interpolacion cuadratica (descartada).

#### `capacitor/scripts/WALL-DIRICHLET/correr_dospot_electrolito.sh`
> Runs the four two-potential electrolyte batteries (lambda_D = 8, 10, 12, 14) in parallel.

#### `capacitor/scripts/WALL-DIRICHLET/correr_electrolito.sh`
> Corre las 6 baterias C-electrolito-* una tras otra, 3 posiciones en paralelo. Reanudable.

#### `capacitor/scripts/WALL-DIRICHLET/correr_lineal.sh`
> Test de linealidad: 4 baterias de 3 o 1 posiciones, 4 corridas en paralelo.

#### `capacitor/scripts/WALL-DIRICHLET/correr_tests.sh`
> Reproduce todas las pruebas de las paredes a potencial constante y regenera
> los graficos. Usa src/Ludwig.exe (enlazado en cada carpeta) y el .petscrc de
> GPU (CG + BoomerAMG). Duracion: A, segundos; regresion, ~3 min;
> B, ~10 min con la GPU libre; C, ~5 min por kernel.
> Uso: ./correr_tests.sh            (todo)
> ./correr_tests.sh A B        (solo esas)

#### `capacitor/scripts/WALL-DIRICHLET/crear_baterias_electrolito.py`
> Baterias C-electrolito-hann{6,8}-LD{4,6,8}: una particula subgrid q=+1 fija en
> x=y=16, a distancia d de la pared inferior, entre paredes a tierra, CON sal 1:1
> (lambda_D = 4, 6, 8). Una corrida por posicion hasta el estacionario ionico.
> Se parte de C-sin-fondo-hann08/ (misma caja 32x32x34, L=32, mismo bloque de particula).
> Uso: python3 crear_baterias_electrolito.py

#### `capacitor/scripts/WALL-DIRICHLET/crear_continuacion.py`
> Continuacion de C-electrolito-hann8-LD6-desplazada (crear_desplazada.py) desde el paso 5000 hasta el 10000, para
> verificar el estacionario. Las corridas originales no se tocan: cada z_* se copia entero (con las salidas del paso 5000:
> dist, rho, vel, psi, qsi, config.cds) a C-electrolito-hann8-LD6-desplazada-cont/, y en la copia el input pasa a
>   N_start 5000, N_cycles 5000 (pasos adicionales), freq_config 1000 (psi/efield cada 1000 pasos).
> Ludwig reanuda leyendo esos archivos (ludwig.c: "Re-starting simulation at step 5000"; colloids_rt.c lee
> config.cds00005000). Probado el 2026-11-05 con 1 paso: max|psi(5001) - psi(5000)| = 2e-15, y los archivos del paso 5000
> no se modifican.
> correr_cont.sh NO borra salidas previas (correr.sh si: borraria los archivos del paso 5000 que hacen falta para
> reanudar). Escribe run_cont.log.
> Las posiciones ya continuadas se saltean (se pueden agregar posiciones nuevas despues).
> Uso: python3 crear_continuacion.py   y despues   cd C-electrolito-hann8-LD6-desplazada-cont && bash correr_cont.sh 3 > bateria.log

#### `capacitor/scripts/WALL-DIRICHLET/crear_desplazada.py`
> Bateria C-electrolito-hann8-LD6-desplazada: el mismo sistema que C-electrolito-hann8-LD6/z_17.50 (Hann-8, lambda_D = 6,
> dx = 1, placas a tierra, q = +1 fija en x = y = 16), con la carga en el centro del canal y corrida en z hacia la placa
> superior por fracciones de celda:
>   z = 17.50 + delta,  delta = 0, 0.10, 0.25, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90
> (z = 18.50, delta = 1, seria de nuevo un punto medio entre nodos, simetrico; no se corre). Las corridas que ya
> existen no se recrean ni se modifican. Las de delta > 0.5 (agregadas despues) guardan psi/efield cada 1000 pasos
> (freq_config 1000, solo cambia la salida) para verificar el estacionario en el paso 5000 con psi(5000) - psi(4000).
> Acercamiento a la placa (2026-11-05): 10 posiciones de delta = 1 a 12, equiespaciadas y redondeadas a 0.1 (= dx/10):
>   delta = 1.0, 2.2, 3.4, 4.7, 5.9, 7.1, 8.3, 9.6, 10.8, 12.0
> delta = 12 es la mas cercana sin que el soporte del kernel, [z0 - 4, z0 + 4], cruce la placa superior (z = +16,
> a mitad del enlace): el borde del soporte queda justo en la placa, donde W = 0.
> Puntos medios entre esas posiciones (2026-11-05), multiplos de 0.1 (4.05 -> 4.0 y 8.95 -> 9.0, eleccion del usuario):
> … (11 líneas más en el archivo)

Uso (del archivo):
    Uso: python3 crear_desplazada.py   y despues   cd C-electrolito-hann8-LD6-desplazada && bash correr.sh 3 > bateria.log
    python3 crear_desplazada.py N=10000   las mismas corridas desde cero con N_cycles = 10000 y psi cada 1000 pasos

#### `capacitor/scripts/WALL-DIRICHLET/crear_desplazada_p4.py`
> Carga desplazada con Hann cos^4 de sigma (casi) igual al Hann-8 cos^2 (2026-11-09, pedido del usuario): prueba en
> delta = 1 con dx = 1 y dx = 2, para la convergencia y para comparar con las corridas cos^2 de la misma posicion.
>
> Ancho: sigma^2 = n^2 (1/12 - 1/(2 pi^2)) para cos^2 y n^2 (1/12 - 5/(8 pi^2)) para cos^4 (comentario de d_hann, subgrid.c).
>   cos^2, n = 8:  sigma = 1.44607 (fisico)
>   cos^4: n = 10.223 daria el mismo sigma, pero n debe ser entero (suma 1 exacta). n = 10 (eleccion del usuario):
>          sigma = 1.41449 fisico (2.2 % menor); n = 11 daria 1.55594 (7.6 % mayor).
> Escalado con el refinamiento igual que cos^2 (ancho fisico fijo 10): dx = 1 -> n = 10, dx = 2 -> n = 5 (impar: permitido,
> la suma 1 sigue exacta para n entero >= p/2 + 1 = 3; solo cambia entre 5 y 4 el numero de nodos activos con el offset).
>
> Cada corrida es una copia de la corrida cos^2 de la misma delta y la misma red; solo cambian subgrid_hann_order,
> subgrid_hann_power (nuevo, 4), pm_sr_sigma (no se usa: pm_sr_correction no) y la cabecera. Mismo Ludwig.exe (se copia
> … (5 líneas más en el archivo)

Uso (del archivo):
    Uso: python3 crear_desplazada_p4.py   y despues, en cada bateria,   bash correr.sh 1 > bateria.log

#### `capacitor/scripts/WALL-DIRICHLET/crear_desplazada_ref.py`
> Carga desplazada (crear_desplazada.py) en la red refinada dx = 1/2: mismo sistema fisico con Hann-16 (ancho fisico 8),
> con la misma receta de escalado que crear_refinado.py para r = 2:
>   - caja 64 x 64 x 66 (32r x 32r x (32r + 2)), carga en x = y = 16r = 32
>   - electrokinetics_epsilon 1e4/r, electrokinetics_init_rho_el 1.388889e-3/r^3 (lambda_D = 6r = 12 nodos)
>   - D por paso igual -> N = 5000 r^2 = 20000
> Posiciones (eleccion del usuario, 2026-11-05): extremos de la serie de acercamiento a la placa y dos intermedias lo mas
> equiespaciadas posible entre las ya corridas con dx = 1 (ideales 4.67 y 8.33):
>   delta = 1.0, 4.7, 8.3, 12.0   (z0 = delta fisico desde el plano medio; Ludwig z = 1.5 + r (16 + delta))
> Para verificar el estacionario psi/efield se guardan cada 4000 pasos (freq_config 4000): mismo tiempo fisico que 1000
> pasos con dx = 1 (el paso de tiempo fisico escala como 1/r^2); criterio sobre psi(20000) - psi(16000).
> Se parte del input de C-electrolito-hann8-LD6/z_17.50 (el mismo de las corridas desplazadas con dx = 1).
> Uso: python3 crear_desplazada_ref.py   y despues   cd C-electrolito-hann16-LD12-ref2-desplazada && bash correr.sh 2 > bateria.log
> … (4 líneas más en el archivo)

Uso (del archivo):
    python3 crear_desplazada_ref.py r=0.5   red dx = 2 (2026-11-06): Hann-4, caja 16 x 16 x 18, lambda_D = 3 nodos,

#### `capacitor/scripts/WALL-DIRICHLET/crear_refinado.py`
> Refinamiento de la red para C-electrolito-hann8-LD6/z_09.50 (carga a d = 8 de la placa inferior, placas a tierra).
>
> Mismo sistema fisico con paso de red dx = 1/r (r = n/8), de modo que Hann-n tenga el mismo ancho fisico (8) que Hann-8:
>   - longitudes x r: caja 32 -> 32r (x, y), canal L = 32 -> 32r (size_z = 32r + 2), carga en (16r, 16r, 1.5 + 8r)
>   - epsilon -> epsilon/r : q/(eps kT) en unidades de red = 10 r, asi psi(x_fisico) es el mismo (l_B -> r l_B)
>   - rho_el -> rho_el/r^3 : densidad numerica por nodo; con eps/r da lambda_D -> r lambda_D
>   - potenciales de pared, q, kT: sin cambios (psi = e phi / kT es adimensional)
>   - D (por paso) igual -> el tiempo de difusion en pasos crece como r^2: N = 5000 r^2 redondeado hacia arriba
> Solo se compara el estacionario, que no depende de D ni de la viscosidad (particula fija, u ~ 0).
>
> Uso: python3 crear_refinado.py [n ...]   (por defecto n = 10 12; p. ej. "python3 crear_refinado.py 16" crea
>      C-electrolito-hann16-LD12-ref2)
> … (14 líneas más en el archivo)

Uso (del archivo):
    python3 crear_refinado.py gauss [n ...]   kernel gaussiano en lugar de Hann, con r = n/8 (n = 8 incluido):
    python3 crear_refinado.py [gauss] q=0.25 [n ...]   carga q0 = 0.25 en lugar de 1 (sufijo -q0.25 en el nombre):
    python3 crear_refinado.py peskin6 [q=...] [n ...]   kernel de Peskin de 6 puntos (Bao et al. 2016, C^3) estirado por
    python3 crear_refinado.py p=4 [q=...] [n ...]   Hann generalizado w ~ cos^p (subgrid_hann_power p, C^(p-1)), mismo

#### `capacitor/scripts/WALL-DIRICHLET/formula_imagenes_kernel.py`
> Prueba C: campo sobre la particula entre placas a tierra, con el kernel Hann-n, por SUMA DE IMAGENES
> puntuales + correccion de hexadecapolo (cumulante kappa4 del kernel). Formula analitica para comparar con
> la simulacion (ver README, 'Formula analitica para el kernel Hann').
>
>   phi(r) = l_B q [ 1/r + (kappa4/24) sum_i d_i^4 (1/r) + O(r^-7) ]        (fuera del soporte)
>   kappa4 = m4 - 3 m2^2 = n^4 (3/(4 pi^4) - 1/120),  m2 = n^2 (1/12 - 1/(2 pi^2))
>   E_z(d) = -l_B q sum_{(i,j,n)} s d_Z [ 1/R + kappa4 (35 (X^4+Y^4+Z^4) - 21 R^4)/(4 R^9) ]
> con imagenes +q en (i l, j l, d+2nL) (salvo la propia) y -q en (i l, j l, -d+2nL); Z = d - z_imagen.
> Uso: python3 formula_imagenes_kernel.py [n=8]

Lee argumentos de línea de comandos: sys.argv[1]

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_campo_dospotenciales.py`
> Prueba C (sin fondo): potencial sobre la linea de la particula y carga
> inducida, medidos, contra una carga GAUSSIANA entre placas a tierra infinitas.
> Gaussiana isotropa con la misma varianza que el kernel Hann-n:
>   sigma^2 = n^2 (1/12 - 1/(2 pi^2)) -> sigma = 1.446 para n = 8.
> Potencial de una gaussiana en el vacio: psi = l_B q erf(r/(sqrt2 sigma))/r
> (valor central l_B q sqrt(2/pi)/sigma). Imagenes: la misma gaussiana, +q en
> d+2nL y -q en -d+2nL. Referencia adicional: caja periodica (Fourier, cada modo
> lateral multiplicado por exp(-k^2 sigma^2/2), convolucion con la gaussiana en z).
> Paredes a DOS potenciales (V1 abajo, V2 arriba): por linealidad, psi = V1 + (V2-V1) z/L +
> psi_carga(paredes a tierra). Corridas reales en C-dospotenciales-hann08/.
> Uso: python3 grafico_C_campo_dospotenciales.py -> ../graficos/C_campo_dospotenciales.png

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_campo_kernel.py`
> Prueba C, paredes a DOS potenciales (C-dospotenciales-hann08): campo medido
> contra la teoria cuyo objetivo es la FORMA REAL DEL KERNEL, no una puntual ni
> una gaussiana: carga q w(x) w(y) w(z), w(r) = (1/n)(1 + cos(2 pi r/n)), |r| <= n/2
> (Hann-n, src/subgrid.c: d_hann(x) d_hann(y) d_hann(z)), n = 8.
>
> Teoria (Fourier lateral, Dirichlet en z, constante de Green de la sec. C del README):
>   psi(r) = (4 pi l_B q / A) sum_k e^{ik.r} wh(kx) wh(ky) int G_k(z,z') w(z') dz'
>   wh(k) = int w(x) cos(kx) dx = sinc(u)/(1-u^2), u = k n/(2 pi)
>   G_k(z,z') = sinh(k z<) sinh(k(L-z>)) / (k sinh kL),  G_0 = z<(L-z>)/L
> Campo sobre la distribucion (fuerza/q; el autocampo se cancela por simetria):
>   E = E0 - (4 pi l_B/A) sum_k wh(kx)^2 wh(ky)^2 II w(z) w(z') d_z G_k(z,z') dz dz'
> Caja de lado ell = 32 (la de la simulacion) o ell = 200 (placas infinitas en la
> … (2 líneas más en el archivo)

Uso (del archivo):
    Uso: python3 grafico_C_campo_kernel.py -> ../graficos/C_campo_dospotenciales_kernel.png

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_carga_desplazada_puntual.py`
> Potential and field at the centre of the displaced charge, all positions and the three lattices (2026-11-08, user's
> request; same runs as grafico_C_convergencia_desplazada.py todos): dx = 2 (Hann-4), 1 (Hann-8), 1/2 (Hann-16).
>
> Ludwig: psi and E_z (efield dump) of the last step on the vertical line x = y = x_charge (the charge is on a node in x, y
> in the three lattices), linearly interpolated in z at z0 = delta. Trilinear interpolation at r0 reduces to this. The
> interpolation adds its own error (user's choice, 2026-11-08: pointwise value at r0, interpolated).
> Reference: linear Debye-Hueckel solution with the continuous Hann source (grafico_C_hoja_hann.py) at r0:
>   phi_DH(r0) = HH.phi(0, delta, 0, 0, mu),   E_DH,z(r0) = -d phi_DH/dz (closed form, HH.dphi via grafico_C_desplazada.ez_an).
> Not separated: the nonlinear Poisson-Boltzmann part (Ludwig is nonlinear, phi_DH linear).
> Figure (one column per quantity): value of Ludwig over the analytic curve; (Ludwig - DH)/DH in % with sign;
> observed order p = log2 |e(coarse)/e(fine)| at the same delta. delta = 0: E = 0, relative error of E undefined.
> dx = 1/2 runs that are not finished are left out.
> … (1 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_carga_desplazada_puntual.py   -> graficos/C_carga_desplazada_puntual.png

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_convergencia_desplazada.py`
> Convergence of the displaced charge between dx = 1 (Hann-8, C-electrolito-hann8-LD6-desplazada, crear_desplazada.py)
> and dx = 1/2 (Hann-16, C-electrolito-hann16-LD12-ref2-desplazada, crear_desplazada_ref.py): same physical system and
> the same physical kernel (continuous Hann of width 8), delta = 1.0, 4.7, 8.3, 12.0 (z0 = delta from the mid-plane).
>
> Errors (all against the analytic linear Debye-Hueckel solution with the continuous Hann charge, as in
> grafico_C_desplazada.py / grafico_C_fuerza_desplazada.py):
>   psi, E_z  on the z axis through the charge: d = Ludwig - analytic at each node. The z nodes of the two lattices do
>             not coincide (dx = 1: z = +-0.5, +-1.5, ...; dx = 1/2: z = +-0.25, +-0.75, ...), so the comparison is by
>             bands (user's choice, 2026-11-06): RMS of d on bands of width 2 with ends at z + 16 = 2j (2 nodes with
>             dx = 1, 4 nodes with dx = 1/2; grafico_C_refinamiento_anidado.band_rms).
>   F         F_L - F_kernel, F_L = Esub_z * R (Ludwig), F_kernel = int W E_DH (continuous limit, the target force).
> Observed order between the two lattices: p = log2(e(dx = 1) / e(dx = 1/2)); with two lattices it is a single number
> … (17 líneas más en el archivo)

Uso (del archivo):
    grafico_C_desplazada.py. Steady state of every dx = 2 run checked here as above. Usage: python3
    Usage: python3 grafico_C_convergencia_desplazada.py

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_desplazada.py`
> Charge at the centre of the channel, displaced towards the upper plate by fractions of a cell (Hann-8, dx = 1,
> lambda_D = 6, plates grounded; runs C-electrolito-hann8-LD6-desplazada built by crear_desplazada.py):
>   z0 = delta (physical, from the mid-plane),  delta = 0, 0.10, 0.25, 0.40, 0.50, 0.60, 0.75, 0.90
> delta = 0 is the mid-plane (half way between two nodes) and delta = 0.5 is a node; 0.10, 0.25, 0.40 are not lattice
> symmetry positions in z.
>
> On the z axis through the charge (x = y = 16, on nodes) and for each delta:
>   psi residual   d_psi(z) = psi_Ludwig(z) - phi_DH(z)
>   field residual d_E(z)   = E_z,Ludwig(z) - E_z,DH(z),   E_z,DH = -d phi_DH/dz
> phi_DH: analytic linear Debye-Hueckel solution with the continuous Hann-8 charge centred at z0 = delta
> (grafico_C_hoja_hann.py, its globals Z0, ZA, ZB set to the displaced charge). E_z,DH: closed-form derivative of phi_DH
> (grafico_C_hoja_hann.dphi, since 2026-11-06; before, central difference with h = 1e-4, agreement ~5e-11). E_z,Ludwig: Ludwig's own efield dump (27-point stencil), read with
> … (22 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_desplazada.py [cont | celda]

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_dospot_vs_LD.py`
> Effect of the lateral images as the screening decreases: two-potential electrolyte capacitor (plates at psi = 0.1 and 0.3),
> point charge at d = 8 from the lower plate, lambda_D = 6, 8, 10, 12, 14.  Ludwig potential vs the closed-system formula
>
>   phi = phi0 + [kL/(2 tanh(kL/2))] <phi0> [1 - cosh(kz)/cosh(kL/2)],   phi0 = phi_CAP + phi_q^p,   k = 1/lambda_D
> Dashed: phi0 of the single cell (infinite plates). Solid: phi0 summed over lateral images (10 per side).
> Lines rho = 6 and 10 (outside the kernel support).   Usage: python3 grafico_C_dospot_vs_LD.py

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_electrolito.py`
> Prueba C con electrolito (C-electrolito-hann{6,8}-LD{4,6,8}): particula subgrid q = +1 entre paredes a tierra,
> sal 1:1, N iones por especie conservado, estacionario (5000 pasos). Teoria de Debye-Huckel con numero fijo de iones.
>
> Linealizando rho_+- = rho_s [1 -+ (psi - <psi>)] (el promedio <psi> sobre el volumen fluido sale de conservar N_+-):
>    nabla^2 psi - kappa^2 (psi - <psi>) = -4 pi l_B rho_p,   kappa = 1/lambda_D = sqrt(8 pi l_B rho_s)
>  * modos laterales k != 0:   (k^2 + kappa^2) en lugar de k^2 en G_k(z,z') (misma G, q = sqrt(k^2+kappa^2)).
>    Es la G de Zhou, Bruch, Wang (arXiv 2405.13261) Eq. (19) con kappa uniforme y G = 0 en la placa metalica.
>  * modo k = 0: psi0 = phi_a + <psi> [1 - cosh(kappa(z-L/2))/cosh(kappa L/2)], phi_a con G_kappa y
>    <psi> = <phi_a> (kappa l)/tanh(kappa l), l = L/2 (autoconsistencia del promedio).
> Cargas inducidas: Q_inf = -(A/4 pi l_B) psi0'(0), Q_sup = +(A/4 pi l_B) psi0'(L).
> Fuente y pesos de fuerza: kernel Hann-n muestreado en nodos fluidos (tratamiento T0, ver README).
> Uso: python3 grafico_C_electrolito.py -> ../graficos/C_electrolito.png + tabla

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_electrolito_gauss.py`
> Check of the derived solution (electrolyte, closed box) with a GAUSSIAN charge, against
>  (1) an independent numerical solver (lateral Fourier modes + finite differences in z, <phi - phibar> = 0 imposed numerically)
>  (2) the Ludwig potential (C-electrolito-hann8-LD6, Hann-8 kernel; same variance as the gaussian).
>
> Derived solution (plates at z = -/+ L/2 grounded, charge at z = d from the centre, box A = 32^2, Gaussian of std sigma):
>   phi = phi_q(rho,z) + phibar [1 - cosh(kappa z)/cosh(kappa L/2)]
>   phi_q = (q/(pi L eps)) sum_m c_m xi_m(z) I_m(rho) + lateral images,
>   I_m(rho) = int_0^inf k J0(k rho) exp(-k^2 sigma^2/2)/(k^2 + g_m^2) dk,  alpha_m = m pi/L,  g_m^2 = alpha_m^2 + kappa^2
>   c_m = int_{-L/2}^{L/2} g(s) xi_m(s) ds   (projection of the z profile g on the plate eigenfunctions)
>   phibar  = (kappa L/2) coth(kappa L/2) <phi_q>
>
> z profile g (CHANGE 20261030): with TRUNC = True the Gaussian is truncated to [-L/2, L/2] and renormalised,
> … (9 líneas más en el archivo)

Uso (del archivo):
    Reduced units q/eps = 10.   Usage: python3 grafico_C_electrolito_gauss.py

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_error_capas_desplazada.py`
> Error of psi and E at the lattice nodes, displaced charge, three lattices (2026-11-09, user's request): dx = 2 (Hann-4),
> 1 (Hann-8), 1/2 (Hann-16), delta = 1.0, 4.7, 8.3, 12.0 (the positions with the three lattices).
>
> Pointwise error at every fluid node r_i (no interpolation): e_psi(r_i) = psi_L(r_i) - phi_DH(r_i),
> e_E(r_i) = |E_L(r_i) - E_DH(r_i)| (modulus of the vector difference). psi_L, E_L: last dump (efield/kT, as elsewhere).
> Reference: linear Debye-Hueckel solution with the continuous Hann source (grafico_C_hoja_hann.py) at the nodes:
>   phi_DH = cap(z) + U(x, y, z) + fac <U> h(z),  U = (2 pi/A) sum_{|n|,|m| <= NK} psit(kx, ky, z) e^{i(kx x + ky y)}.
>   On the nodes x = j dx the modes n and n + N (N = 32/dx) coincide, so the sum is folded mod N and done with a 2D FFT
>   (exact on the nodes). E_DH,x/y: modes times -i kx, -i ky; E_DH,z: closed-form dpsit (HH.dU), cap' and h'.
>   Check printed: U and dU/dz on the axis through the charge vs HH.U, HH.dU.
> 1) Shells of distance d = |r_i - r0| (minimum image in x, y), width W_SHELL (physical units, parameter w=):
>    RMS(e)/RMS(ref) in each shell (absolute error divided by the RMS of the reference in the same shell; user's choice),
> … (6 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_error_capas_desplazada.py [w=1]

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_error_capas_p4.py`
> Error of psi and E at the lattice nodes near the kernel edge: Hann cos^4 (physical width 10) vs Hann cos^2 (width 8),
> delta = 1, dx = 2 and 1 (2026-11-09, user's request: convergence at the edge of the kernel).
>
> Same metric as grafico_C_error_capas_desplazada.py (whose functions dh_on_nodes, distances, shells are used):
>   e_psi(r_i) = psi_L(r_i) - phi_DH(r_i),  e_E(r_i) = |E_L(r_i) - E_DH(r_i)| at every fluid node, last dump, no interpolation;
>   phi_DH: linear Debye-Hueckel solution of the charge spread with the SAME kernel as the run (cos^4: grafico_C_hoja_hann_pot);
>   shells of d = |r_i - r0| (minimum image in x, y) of width W_SHELL (physical, parameter w=):
>   RMS(e)/RMS(ref) per shell and the order p = log2 RMS(e, dx = 2)/RMS(e, dx = 1) in the same shell.
> Kernel edge: the support is the cube |x|, |y|, |z - z0| < n/2 (tensor product), so in d it starts at n/2 (faces) and ends
> at n sqrt(3)/2 (corners): 4 ... 6.93 for cos^2, 5 ... 8.66 for cos^4 (lines in the figure).
> Smoothness at the edge: cos^2 window C^1 (W'' jumps), cos^4 C^3 (W'''' jumps) (subgrid.c d_hann comment).
> Global norms (all fluid nodes, ||e||_2 = sqrt(sum dx^3 e^2), ||e||_inf) printed with their orders.
> … (2 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_error_capas_p4.py [w=1]   -> graficos/C_error_capas_p4.png

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_formulas.py`
> Ludwig potential vs the analytical formulas (plates at z = -/+ L/2, charge at z = d measured from the centre).
>
> (a) No electrolyte, plates at psi_1 = 0.1 (lower) and psi_2 = 0.3 (upper), Gaussian charge (sigma of Hann-8):
>       phi_TOT = phi_CAP(z) + phi_G(rho, z),   phi_CAP = (p1+p2)/2 + (p2-p1) z/L
>       phi_G   = l_B q sum_n [ erf(r_n^+/(sqrt2 sigma))/r_n^+ - erf(r_n^-/(sqrt2 sigma))/r_n^- ]
>       r_n^+ = sqrt(rho^2 + (z-(2nL+d))^2),  r_n^- = sqrt(rho^2 + (z-((2n+1)L-d))^2)
> (b) 1:1 electrolyte, grounded plates (as written on the sheet, NO mean-potential term):
>       phi_CAP = [p1 sinh k(L/2-z) + p2 sinh k(L/2+z)]/sinh kL ,  phi_TOT = phi_CAP + sum_n [phi^A_n + phi^B_n]
> Data: C-dospotenciales-hann08 and C-electrolito-hann8-LD6 (Hann-8, d = 8 from the lower plate, L = 32, box 32x32).
> Reduced units (psi = e phi / k_B T, 4 pi l_B = 10).   Usage: python3 grafico_C_formulas.py

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_formulas_closed.py`
> Electrolyte capacitor with a point charge: Ludwig potential vs the CLOSED-system (fixed ion number) formula
>
>   phi(rho,z) = phi0(rho,z) + [kL / (2 tanh(kL/2))] <phi0> [1 - cosh(k z)/cosh(k L/2)],     k = kappa
>   phi0 = phi_CAP(z) + sum_n [phi^A_n + phi^B_n]       (open-system solution, point charge)
>   <phi0> = <phi_q> + <phi_CAP>,  <phi_CAP> = (p1+p2) tanh(kL/2)/(kL),
>   <phi_q> = (4q/(pi eps V)) sum_{n>=0} (-1)^n/(2n+1) cos((2n+1) pi d/L)/(((2n+1)pi/L)^2 + k^2),  V = l^2 L
> Dashed: phi0 of the single cell (infinite plates); solid: phi0 summed over lateral images (10 per side).
> Reduced units, q/eps = 10.  Data: C-electrolito-hann8-LD6 (Hann-8, d = 8 from the lower plate).
> Usage: python3 grafico_C_formulas_closed.py -> ../graficos/C_formula_electrolyte_closed_vs_ludwig.png

Lee argumentos de línea de comandos: sys.argv[1, 2, 3, 4]

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_fuerza_desplazada.py`
> Spurious force on the displaced charge (C-electrolito-hann8-LD6-desplazada[-cont], crear_desplazada.py /
> crear_continuacion.py): Hann-8, dx = 1, q = +1 at x = y = 16 (on nodes), z0 = delta from the mid-plane,
> delta = 0, 0.10, 0.25, 0.40, 0.50. Units: reduced field e E/kT per physical length; with q = 1 it is also the force
> in kT per length.
>
> Two quantities, z component (x, y vanish by symmetry):
>   F_L      Ludwig: Esub_z of proceced_data/particle_Esub.csv at the last step. subgrid.c subgrid_update_Esub:
>            Esub = sum_i W(r_i - r0) e(r_i), e = psi_electric_field (27-point stencil, wall links doubled), and the force
>            is q kT/e Esub (subgrid_update_forces_electrokinetics).
>            Check printed: the same sum recomputed here from the efield dump (efield/kT) at that step.
>   F_kernel analytic force distributed by the kernel, the continuum limit of F_L (since 2026-11-06 computed directly,
>            Fbar = -int W(z - z0) dPhibar/dz dz with the closed-form derivative HH.dpsit; the form below, integrated
> … (14 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_fuerza_desplazada.py [cont | celda | placa | ref2]   (same run selection as grafico_C_desplazada.py;

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_fuerza_p4.py`
> Force error of the displaced charge with the Hann cos^4 kernel vs the original Hann cos^2, delta = 1, dx = 2 and 1
> (2026-11-09, user's request).
>
> Runs (crear_desplazada_p4.py; the cos^4 runs are copies of the cos^2 ones with only the kernel changed):
>   cos^2, physical width 8:  dx = 2 C-electrolito-hann4-LD3-ref0.5-desplazada/z_10.00 (Hann-4)
>                             dx = 1 C-electrolito-hann8-LD6-desplazada/z_18.50          (Hann-8)
>   cos^4, physical width 10: dx = 2 C-electrolito-hann5p4-LD3-ref0.5-desplazada/z_10.00 (n = 5)
>                             dx = 1 C-electrolito-hann10p4-LD6-desplazada/z_18.50        (n = 10)
> Each kernel against ITS OWN continuum target (same definitions as grafico_C_fuerza_desplazada.py):
>   F_L      = Esub_z / dx (last line of proceced_data/particle_Esub.csv, physical units)
>   F_kernel = int W(r - r0) E_DH,z(r) d^3r, E_DH the linear Debye-Hueckel field of the charge spread with the same W
>              (grafico_C_hoja_hann_pot.py: grafico_C_hoja_hann.py generalised to cos^p)
> … (11 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_fuerza_p4.py   (text output only)

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_gaussiana.py`
> Prueba C (sin fondo): potencial sobre la linea de la particula y carga
> inducida, medidos, contra una carga GAUSSIANA entre placas a tierra infinitas.
> Gaussiana isotropa con la misma varianza que el kernel Hann-n:
>   sigma^2 = n^2 (1/12 - 1/(2 pi^2)) -> sigma = 1.446 para n = 8.
> Potencial de una gaussiana en el vacio: psi = l_B q erf(r/(sqrt2 sigma))/r
> (valor central l_B q sqrt(2/pi)/sigma). Imagenes: la misma gaussiana, +q en
> d+2nL y -q en -d+2nL. Referencia adicional: caja periodica (Fourier, cada modo
> lateral multiplicado por exp(-k^2 sigma^2/2), convolucion con la gaussiana en z).
> Uso: python3 grafico_C_gaussiana.py -> ../graficos/C_gaussiana_vs_placas_infinitas.png

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_gaussiana_dospotenciales.py`
> Prueba C (sin fondo): potencial sobre la linea de la particula y carga
> inducida, medidos, contra una carga GAUSSIANA entre placas a tierra infinitas.
> Gaussiana isotropa con la misma varianza que el kernel Hann-n:
>   sigma^2 = n^2 (1/12 - 1/(2 pi^2)) -> sigma = 1.446 para n = 8.
> Potencial de una gaussiana en el vacio: psi = l_B q erf(r/(sqrt2 sigma))/r
> (valor central l_B q sqrt(2/pi)/sigma). Imagenes: la misma gaussiana, +q en
> d+2nL y -q en -d+2nL. Referencia adicional: caja periodica (Fourier, cada modo
> lateral multiplicado por exp(-k^2 sigma^2/2), convolucion con la gaussiana en z).
> Paredes a DOS potenciales (V1 abajo, V2 arriba): por linealidad, psi = V1 + (V2-V1) z/L +
> psi_carga(paredes a tierra). Corridas reales en C-dospotenciales-hann08/.
> Uso: python3 grafico_C_gaussiana_dospotenciales.py -> ../graficos/C_gaussiana_dospotenciales.png

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_hoja_gauss.py`
> Electrolyte capacitor (closed system) with a GAUSSIAN charge: the hand-derived solution of the sheet (2026-10-30),
> evaluated term by term, vs Ludwig.
>
>   phi   = phi_C(z) + U(x,y,z) + [kL/(2 tanh(kL/2))] <U> [1 - cosh(k z)/cosh(k L/2)],   k = kappa
>   phi_C = V_m + (DV/2) sinh(k z)/sinh(k L/2),  V_m = (p1+p2)/2, DV = p2 - p1   (capacitor of the CLOSED system, see below)
>   U     = (2 pi/A) sum_nm  psit(k_nm, z) exp(i k_nm . r),     k_nm = (2 pi/l)(n, m)      (exact periodic images)
>   <U>   = (2 pi/(A L)) int_{-L/2}^{L/2} psit(0, z) dz                                    (k = 0 mode, by quadrature)
>   psit  = q exp(-s^2 k^2/2) / (2 pi eps a sinh(aL)) [sinh(a(L/2 - z)) I_< + sinh(a(L/2 + z)) I_>],   a^2 = k^2 + kappa^2
>   I_<   = (1/C1) int_{-L/2}^{z}  exp(-(s'-z0)^2/(2 s^2)) sinh(a(s'+L/2)) ds'
>   I_>   = (1/C1) int_{z}^{L/2}   exp(-(s'-z0)^2/(2 s^2)) sinh(a(L/2-s')) ds'
>   C1    = int_{-L/2}^{L/2} exp(-(s'-z0)^2/(2 s^2)) ds'        (Gaussian truncated to the channel, renormalised)
>   phi_CAP = [p1 sinh k(L/2-z) + p2 sinh k(L/2+z)]/sinh kL   (open system: phi'' - k^2 phi = 0)
> … (15 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_hoja_gauss.py            (grounded + two-potential figures)
    LAYOUT=left python3 grafico_C_hoja_gauss.py  (slide version, formulas at the left -> *_ppt.png)

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_hoja_hann.py`
> Electrolyte capacitor (closed system): the solution of the sheet (grafico_C_hoja_gauss.py) with the charge spread by
> Ludwig's Hann-8 kernel instead of a Gaussian, vs Ludwig and vs the Gaussian solution.
>
> Charge density (continuous version of Ludwig's spread, subgrid.c d_hann):
>   rho_part = q W(x) W(y) W(z - z0),   W(u) = (1/n)[1 + cos(2 pi u/n)] for |u| <= n/2, 0 otherwise,  n = 8
> Same derivation as the sheet; only the source transform changes:
>   phi   = phi_C + U + [kL/(2 tanh(kL/2))] <U> [1 - cosh(k z)/cosh(k L/2)]    (phi_C: closed-system capacitor, as there)
>   U     = (2 pi/A) sum_nm psit(kx, ky, z) exp(i k.r)
>   psit  = q What(kx) What(ky) / (2 pi eps a sinh(aL)) [sinh(a(L/2-z)) I_< + sinh(a(L/2+z)) I_>],  a^2 = kx^2+ky^2+kappa^2
>   What(k) = [sin(x)/x] / [1 - (x/pi)^2],  x = k n/2       (2D transform is NOT radial: depends on kx, ky separately)
>   I_<   = int_{-L/2}^{z} W(s-z0) sinh(a(s+L/2)) ds,   I_> = int_{z}^{L/2} W(s-z0) sinh(a(L/2-s)) ds
>   int e^{b s} [1 + cos(w(s-z0))] ds = e^{b s} [1/b + (b cos(w(s-z0)) + w sin(w(s-z0)))/(b^2 + w^2)],   w = 2 pi/n
> … (8 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_hoja_hann.py ;  LAYOUT=left python3 grafico_C_hoja_hann.py  (slide version, *_ppt.png)

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_hoja_hann_pot.py`
> Analytic linear Debye-Hueckel reference of grafico_C_hoja_hann.py generalised to the Hann window cos^p (2026-11-09,
> user's request: reference for the cos^4 kernel, subgrid.c d_hann with subgrid_hann_power 4).
>
> Kernel (subgrid.c, CHANGE 20261104): w(u) = cos^p(pi u/n) / (n C_p), |u| < n/2, C_p = (p-1)!!/p!!.
> For even p, cos^p(theta) = 2^-p [binom(p, p/2) + 2 sum_{m=1}^{p/2} binom(p, p/2 - m) cos(2 m theta)], so with
> omega = 2 pi/n
>   W(u) = (1/n) sum_{m=0}^{p/2} c_m cos(m omega u),   c_0 = 1,  c_m = 2 binom(p, p/2 - m)/binom(p, p/2)
>   p = 2: c = (1, 1)   [(1/n)(1 + cos(omega u)), the original Hann]
>   p = 4: c = (1, 4/3, 1/3)
> Only the source changes in the derivation of grafico_C_hoja_hann.py; the three kernel-dependent pieces become sums over
> the harmonics m omega:
>   What(k) = sum_m (c_m / (1 + [m > 0])) [sinc(x - m pi) + sinc(x + m pi)] with x = k n/2, sinc(y) = sin(y)/y
> … (8 líneas más en el archivo)

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_hoja_hannp.py`
> Electrolyte capacitor (closed system): the solution of grafico_C_hoja_hann.py with Ludwig's generalised Hann window
> w ~ cos^p (subgrid_hann_power p; subgrid.c d_hann), support n = 8 physical, p even (here p = 4: C^3).
>
> Charge density (continuous version of Ludwig's spread):
>   rho_part = q W(x) W(y) W(z - z0),   W(u) = cos^p(pi u/n) / (n C_p) for |u| < n/2,   C_p = (p-1)!!/p!!
> For even p, cos^p(theta) = 2^-p [binom(p, p/2) + 2 sum_{m=1}^{p/2} binom(p, p/2 - m) cos(2 m theta)], so
>   W(u) = sum_{m=0}^{p/2} c_m cos(m w u),   w = 2 pi/n
> (p = 2: c = (1, 1)/n, the Hann sheet; p = 4: c = (3, 4, 1)/(3n)). Every term has elementary integrals, so the sheet's
> closed forms carry over term by term:
>   What(k) = sum_m c_m (n/2) [sinc((k - m w) n/2) + sinc((k + m w) n/2)],   sinc(x) = sin(x)/x
>   int e^{b s} cos(m w (s - z0)) ds = e^{b s} (b cos(m w (s - z0)) + m w sin(m w (s - z0)))/(b^2 + m^2 w^2)
> and psit, U, <U>, phi exactly as in grafico_C_hoja_hann.py (Green function in z: Jackson 3rd ed., Problem 3.17).
> … (5 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_hoja_hannp.py   (checks; used as a module, P = 4 by default: set_power(p) to change)

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_hoja_peskin6.py`
> Electrolyte capacitor (closed system): the solution of grafico_C_hoja_hann.py with the charge spread by Ludwig's
> STRETCHED Peskin 6-point kernel (subgrid_kernel peskin6, subgrid_peskin6_stretch s; subgrid.c d_peskin6).
>
> Charge density (continuous version of Ludwig's spread, in physical units = lattice units of the dx = 1 run):
>   rho_part = q W(x) W(y) W(z - z0),   W(u) = phi6(u/S)/S,   S = 2 (s = 2r lattice units = 2 physical for every r)
>   phi6 = Bao, Kaplan & Peskin (2016) 6-point kernel, support |v| < 3, C^3 (J. Comput. Phys. 347 (2017) 183;
>   arXiv:1505.07529, Eqs. 2.16-2.23 as coded in subgrid.c peskin6_phi_segment). Support |u| < 3S = 6, sigma = S sqrt(K)
>   = 1.69006 with K = 59/60 - sqrt(29)/20 (its second-moment condition).
> Same derivation as the Hann sheet; only the source transform changes:
>   phi   = phi_C + U + [kL/(2 tanh(kL/2))] <U> [1 - cosh(k z)/cosh(k L/2)]
>   U     = (2 pi/A) sum_nm psit(kx, ky, z) exp(i k.r)
>   psit  = q What(kx) What(ky)/(2 pi eps) int W(s - z0) g_a(s, z) ds,   a^2 = kx^2 + ky^2 + kappa^2
> … (10 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_hoja_peskin6.py   (checks only; used as a module by grafico_C_refinamiento_peskin6.py)

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_hoja_potencial.py`
> Potencial electrostatico de Ludwig (carga q = +1 subgrid, Hann-8, entre placas a tierra, electrolito 1:1)
> contra la solucion analitica de la hoja (separacion de variables, Fourier en z + K0 en rho), CORREGIDA:
>
>   phi_q(rho, z) = (q/(pi L eps)) sum_{m>=1} sen[m pi (d+L/2)/L] sen[m pi (z+L/2)/L] K0( sqrt((m pi/L)^2 + kappa^2) rho )
>   phi_TOT(abierta)  = phi_1 sinh k(L/2-z)/sinh kL + phi_2 sinh k(L/2+z)/sinh kL + phi_q     (placas en +-L/2)
>   phi_TOT(cerrada)  = phi_TOT(abierta) + psibar * [1 - cosh(kappa z)/cosh(kappa L/2)],
>         psibar = (kappa L/2) coth(kappa L/2) <phi_q>,  <phi_q> = (1/A L) int phi_q dV
>                = (q/(pi L eps)) (2 pi/(A L)) sum_{m impar} sen(m pi (d+L/2)/L) (2L/(m pi)) / g_m^2,  g_m^2 = (m pi/L)^2 + kappa^2
> La "abierta" supone sal en contacto con un reservorio (n0 a potencial cero); Ludwig conserva el numero de iones en la caja
> de area A = 32^2, y eso agrega el desplazamiento psibar (orden 1/A) y h(z). Mas replicas laterales (|i|,|j| <= 2).
> Unidades reducidas: q/eps = 4 pi l_B q = 10.   Uso: python3 grafico_C_hoja_potencial.py

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_ion_conc.py`
> Ion concentrations n+/n0 and n-/n0 of the closed-system electrolyte capacitor (plates at psi=0.1 and 0.3, lambda_D=6,
> point charge q/eps=10 at d=8 from the lower plate) -- Ludwig (qsi file, step 5000) vs the linearised closed-system theory
>
>    n_pm = n0 [ 1 -/+ (phi - <phi>) ],   phi = phi0 + sum of lateral images + closed-system correction  (same phi as the potential figure)
>
> Reduced units psi = e phi / kBT.  n0 = mean ion density (= initial rho_el, ion number conserved).
> Usage: python3 grafico_C_ion_conc.py -> ../graficos/C_ion_concentration_closed_2pot.png

Lee argumentos de línea de comandos: sys.argv[0]

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_placas_infinitas.py`
> Prueba C (sin fondo): campo medido sobre la particula vs el del sistema ideal
> de una carga puntual entre dos placas a tierra INFINITAS (sin periodicidad
> lateral): serie de imagenes a lo largo de z, +q en d+2nL (n!=0), -q en -d+2nL.
> E_z = l_B * sum q_i (d - z_i)/|d - z_i|^3, unidades de Ludwig (4 pi l_B = 10).
> Uso: python3 grafico_C_placas_infinitas.py   -> ../graficos/C_vs_placas_infinitas.png

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_potencial_carga.py`
> Prueba C (sin fondo): potencial sobre la linea de la particula y carga
> inducida en las paredes, medidos, contra el sistema ideal de una carga puntual
> entre dos placas a tierra INFINITAS (imagenes solo en z). Se agrega la teoria
> periodica en x,y de la caja (Fourier) como referencia.
> psi reducido = l_B * sum q_i/|r - r_i|, l_B = 10/(4 pi). Particula en x=y=16 (nodo).
> Uso: python3 grafico_C_potencial_carga.py -> ../graficos/C_potencial_carga_vs_placas_infinitas.png
>
> 2026-11-16 (pedido del usuario): en el panel (a) la referencia teórica es ahora la carga repartida con el kernel
> Hann-8 CONTINUO (psi_hann), no la carga puntual: misma serie de Fourier en x, y que psi_periodic (placas a tierra,
> caja periódica 32x32, sin fondo), con la fuente transformada  W^(kx) W^(ky),  W^(k) = [sin(kn/2)/(kn/2)] / [1-(kn/2pi)^2]
> (la misma de grafico_C_hoja_hann.py), y convolución en z con W(s - d):
>    psi(z) = (4 pi l_B / A) sum_k W^(kx) W^(ky) int W(s-d) G_k(z, s) ds,
> … (3 líneas más en el archivo)

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_refinamiento.py`
> Lattice refinement of the electrolyte capacitor with one Gaussian-like charge: Ludwig with Hann-8, -10, -12, -16 on
> lattices dx = 1, 1/1.25, 1/1.5, 1/2 (same physical system, built by crear_refinado.py) vs the analytic solution with the
> continuous Hann charge (grafico_C_hoja_hann.py), whose physical width n dx = 8 is the same for all runs.
>
> Convergence panel: residual Ludwig - analytic on the z axis through the charge (rho = 0) at FIXED physical distances
> z - z0 (cubic spline through the nodes: inside the support with the nodes |z - z0| < 4, outside with the nodes
> z - z0 > 4 on the channel side), so every run is compared at the same points.
>
> Physical units: L = 32, lambda_D = 6, z0 = -8 (d = 8 from the lower plate), q/(eps kT) = 10, plates grounded.
> Ludwig node k (1-based, z) is at physical z = (k - 1.5)/r - L/2; the charge sits at x = y = 16 r.
> Field line y = 0 (through the charge), x = rho: rho_lat = round(r rho), evaluated at the exact rho_lat/r.
> Usage: python3 grafico_C_refinamiento.py ;  LAYOUT=left python3 grafico_C_refinamiento.py  (slide version, *_ppt.png)

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_refinamiento_anidado.py`
> Nested lattice refinement with the Hann kernel: dx = 2, 1, 1/2 (Hann-4, -8, -16; same physical system, built by
> crear_refinado.py 4 / the Hann-8 base / crear_refinado.py 16), local convergence order from the RMS of the nodal
> residual on bands (rho = 0, z axis through the charge), vs the analytic continuous-Hann solution (grafico_C_hoja_hann.py).
>
> Why bands and not nodes: the wall sits at half a link, so nodes are cell centres z + 16 = (m + 1/2) dx and never
> coincide between lattices (halving dx splits each cell in two, with new centres at 1/4 and 3/4 of the old cell).
> What does nest are the CELL FACES z + 16 = m dx: every face of dx = 2 is a face of dx = 1 and of dx = 1/2.
> A band whose ends are faces of the coarse lattice (even z + 16) is tiled by whole cells in all three lattices, so
>     RMS_[a,b] = ( (1/(b-a)) sum_{z_k in [a,b]} d_k^2 dx )^(1/2)
> is the midpoint rule of (1/(b-a)) int d^2 dz in every lattice, with no clipped cells and no interpolation; the same
> band covers exactly 1 + ... coarse cells, twice as many medium and four times as many fine ones.
> Two families: width 2 every 2 (one coarse cell per band: contiguous) and width 4 every 2 (overlap 2).
> … (1 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_refinamiento_anidado.py

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_refinamiento_bandas.py`
> Local convergence order from the RMS of the nodal residual on NARROW overlapping bands (rho = 0, z axis through
> the charge), for the Hann-n and Gaussian refinement series (grafico_C_refinamiento.py, grafico_C_refinamiento_gauss.py).
>
> Bands [s, s + W] with W = 2 physical, s every STEP = 1 (consecutive bands overlap by 1).
> Each node stands for its cell [z_k - dx/2, z_k + dx/2] (cell-centred lattice, wall at half a link). The band integral
> is the sum over cells of d_k^2 times the length of the cell that lies inside the band (exact for whole cells, i.e. the
> midpoint rule, and clipped at the two ends):
>     RMS_[s, s+W] = ( (1/W) sum_k d_k^2 |cell_k ∩ [s, s+W]| )^(1/2)
> Unlike the 4-unit bands, the ends of a 2-unit band are not cell faces of every lattice (common faces of dx = 1, 0.8,
> 2/3, 1/2 are 4 apart), so the clipped end cells add a quadrature error O(dx^2) that depends on where the end falls in
> the cell. Its size is measured with a smooth test function (QUAD check below) and printed: up to 10 % for dx = 0.8 and
> 2/3, enough to fake orders of 3 or -4 where the residual changes by 10x inside a band (support edges).
> … (3 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_refinamiento_bandas.py

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_refinamiento_gauss.py`
> Smoothness control for the order ~1 near the edge of the Hann support (grafico_C_refinamiento.py): the same lattice
> refinement with Ludwig's truncated Gaussian kernel (subgrid_kernel gauss; runs built by "crear_refinado.py gauss"),
> sigma = 1.44612 physical (same variance as Hann-8), cut at 5.5 sigma, vs the analytic Gaussian solution
> (grafico_C_hoja_gauss.py, full Gaussian; the cut changes the kernel by exp(-5.5^2/2) = 2.7e-7 relative).
> Hann is only C^1 (W'' jumps at |u| = 4); the Gaussian is smooth. If the order-1 bands of Hann become order 2 with the
> Gaussian, the loss of order comes from the non-smooth edge of the Hann kernel.
>
> Metric: RMS of the nodal residual on the same bands as grafico_C_refinamiento.py (ends on cell faces common to every
> lattice, so no interpolation), rho = 0, z axis through the charge.
> Usage: python3 grafico_C_refinamiento_gauss.py ;  LAYOUT=left python3 grafico_C_refinamiento_gauss.py

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_refinamiento_peskin6.py`
> Nested refinement dx = 2, 1, 1/2 with Ludwig's stretched Peskin 6-point kernel (C^3; subgrid_peskin6_stretch
> s = 2r, fixed physical support |u| < 6, sigma = 1.690), runs built by "crear_refinado.py peskin6 [q=0.25] 4 8 16",
> against the analytic Peskin-6 solution (grafico_C_hoja_peskin6.py); side by side with the Hann-8 runs
> (C^1, support |u| <= 4) of grafico_C_refinamiento_q.py, and with the generalised Hann window cos^4 (subgrid_hann_power 4,
> C^3, the SAME support |u| <= 4 and lattices as Hann-8, sigma = 1.132; "crear_refinado.py p=4 [q=0.25] 4 8 16", analytic
> solution grafico_C_hoja_hannp.py).
>
> Question: is the loss of order near the support edge (order ~1 instead of 2) caused by the jump of W'' of Hann?
> Peskin-6 is C^3 (W'''' is the first derivative that jumps), so if the edge order is a kernel-smoothness effect it
> must go back to 2 here.
>
> Per lattice, band RMS of the nodal residual on bands of width 2 with ends on faces of dx = 2 (whole cells in all
> … (5 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_refinamiento_peskin6.py

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_refinamiento_q.py`
> Linearity test of the dx-independent floor C(z) in the refinement residual (Hann, nested lattices dx = 2, 1, 1/2):
> same runs with charge q = 1 (grafico_C_refinamiento_anidado.py) and q = 0.25 ("crear_refinado.py q=0.25 4 8 16").
>
> The analytic solution is linear Debye-Hueckel, so it scales with q: residual d_q = psi_Ludwig(q) - q phi_lin.
> Ludwig solves Nernst-Planck, whose equilibrium is nonlinear Poisson-Boltzmann (ion densities ~ exp(-+psi)):
> psi = q phi_lin + q^3 phi_3 + ..., so the nonlinear part of d_q / q scales as q^2. If the floor is nonlinearity,
>     C_{0.25}/0.25 = C_1/1 / 16.
> Discretisation errors are linear in the source: (d_q/q) at a fixed dx is the same for both q except for the floor.
>
> Band means (signed): M_[a,b] = (1/(b-a)) sum_{z_k in [a,b]} d_k dx, bands of width 2 with ends on faces of dx = 2
> (nested: whole cells in all three lattices, exact midpoint rule, no interpolation). M is linear in d, so the model
> M = C + a dx^2 holds band by band. C from Richardson with dx = 1 and 1/2: C = (4 M_{1/2} - M_1)/3; the 3-lattice
> … (2 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_refinamiento_q.py

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_refinamiento_zoom.py`
> Zoom of grafico_C_refinamiento.py on |z - z0| <= 2 (physical units) along the z axis through the charge (rho = 0):
> potential psi and field E_z = -dpsi/dz for Ludwig Hann-8, -10, -12, -16 vs the analytic continuous-Hann solution
> (grafico_C_hoja_hann.py), with the residuals below. With the argument "full" the whole channel -16 <= z <= 16.
>
> On the axis E_x = E_y = 0 by symmetry, so E_z is the whole field.
> Ludwig E_z: the field Ludwig itself writes (efield-*, psi_electric_field in psi_gradients.c): 27-point gradient stencil
> (electrokinetics_solver_stencil 27), difference form with the wall links doubled (Asta walls). The file is in lattice
> units of the potential (kT/e) per lattice spacing: E_reduced, physical = efield / kT * r, kT = 1e-5 (input "temperature").
> Third row: relative error |Ludwig - analytic| / |analytic| (%, log scale). It grows where psi -> 0 at the grounded
> plates; at z0 it stays finite because E_z and its error are both odd about the charge (both vanish linearly).
> Fourth row: the same relative error with its sign, (Ludwig - analytic) / analytic in %, on a linear scale.
> For reference the residual panel also shows the plain central difference -(psi_{k+1} - psi_{k-1})/(2 dx) (open, faint).
> … (2 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_C_refinamiento_zoom.py [full]   (-> C_refinamiento_hann_zoom.png / C_refinamiento_hann_full.png)

#### `capacitor/scripts/WALL-DIRICHLET/grafico_C_subkernel.py`
> Prueba C, particula con el soporte del kernel cruzando la pared (d < n/2 - 1/2 = 3.5, Hann-8):
> Ludwig (C-sin-fondo-hann08) contra DOS tratamientos teoricos del kernel que cruza el plano medio:
>   T0 "truncado": solo los nodos fluidos reciben carga y fuerza; la carga que cae en nodos de pared se pierde.
>   T1 "imagen espejada" (DISCOS, Ladiges et al., PRE 106, 035104 (2022), Sec. II.B.1): se agrega una particula
>      imagen de carga opuesta espejada en el plano medio, cuyos pesos en el dominio se restan a los de la real;
>      el campo se interpola con celdas fantasma pares, es decir F = q sum_n [w(z_n-d) + w(z_n+d)] E(z_n).
> Los dos usan el kernel Hann-8 MUESTREADO en los nodos (x, y, z), Fourier lateral con ell = 32 y G_k con Dirichlet en
> z = 0 y L = 32. Carga inducida: reciprocidad de Green (Jackson 1.13) con las cargas nodales de cada tratamiento.
> Uso: python3 grafico_C_subkernel.py -> ../graficos/C_subkernel_pared.png (y tabla por stdout)

#### `capacitor/scripts/WALL-DIRICHLET/grafico_hasta_la_pared.py`
> Free particle (no electrolyte, eta = 0.0005) falling onto the lower grounded wall: full trajectory (steps 0-30000 from
> C-particula-libre-hann08-eta0.0005 and the continuation C-particula-libre-hann08-eta0.0005-hastapared until d < 0.1).
> Panels: (a) d(t); (b) reduced field on the particle vs d, with the fixed-particle values of C-sin-fondo-hann08 (static);
> (c) induced charge on the walls vs d (sum = -q_eff: charge lost when the kernel reaches wall nodes).
> Usage: python3 grafico_hasta_la_pared.py -> ../graficos/particula_libre_hasta_la_pared.png

#### `capacitor/scripts/WALL-DIRICHLET/grafico_paper_asta.py`
> Comparacion DIRECTA con la Fig. 3 de Asta et al. 2019 (capacitor de placas
> paralelas, seccion III.A.1): mismos parametros reducidos del paper
> (L = 76 Dx, lambda_D = 6 Dx, beta*e*Dpsi = 0.1, D = 0.05), corridos con nuestra
> implementacion de paredes a potencial constante. Mismo formato de la figura
> del paper: panel (a) perfil de potencial, panel (b) concentracion de iones.
>
> A diferencia de graficos_paredes.py (pruebas A y B), aca NO se calcula la
> solucion discreta de referencia ni un panel de desviacion: se compara solo
> contra la formula analitica continua (ec. 17 del paper) y la distribucion de
> Boltzmann que de ella se deriva, que es lo unico que el paper mismo grafica
> en su Fig. 3.
>
> … (2 líneas más en el archivo)

Uso (del archivo):
    Uso: python3 grafico_paper_asta.py

#### `capacitor/scripts/WALL-DIRICHLET/grafico_paper_asta_carga.py`
> Comparacion DIRECTA con la Fig. 4a de Asta et al. 2019 (dinamica de carga del
> capacitor de placas paralelas, seccion III.A.3): Q(t)/Q_inf vs t/tau_norm,
> con tau_norm = L*lambda_D/(2D), y el ajuste exponencial
> Q(t) = Q_inf + (Q0 - Q_inf) * exp(-t/tau).
>
> Q(t) se saca de las lineas "Wall charge lower/upper" que Ludwig imprime cada
> freq_psi_resid pasos (ver E-paper-asta-carga/input). Se usa la carga de la
> pared superior (mismo signo que Dpsi>0); Q0 es el primer valor (capacitor
> neutro, antes de que los iones se muevan) y Q_inf el ultimo (estacionario).
>
> Lee <carpeta>/run.log. Escribe graficos/<nombre>.png.
> Uso: python3 grafico_paper_asta_carga.py [carpeta] [nombre_salida] [titulo_extra]
> … (1 líneas más en el archivo)

Lee argumentos de línea de comandos: sys.argv[1, 2, 3]

#### `capacitor/scripts/WALL-DIRICHLET/grafico_paper_asta_carga_sinfit.py`
> Version de grafico_paper_asta_carga.py con SOLO el panel (a) (Q(t)/Q_inf vs
> t/tau_norm) y SIN el ajuste exponencial superpuesto -- solo los datos crudos
> de Ludwig. Reutiliza la lectura de datos de grafico_paper_asta_carga.py.
>
> Uso: python3 grafico_paper_asta_carga_sinfit.py [carpeta] [nombre_salida] [titulo_extra]
>   (por defecto: E-paper-asta-carga-A0 / E_paper_asta_carga_A0_sinfit / "")

Lee argumentos de línea de comandos: sys.argv[1, 2, 3]

#### `capacitor/scripts/WALL-DIRICHLET/grafico_pared_blanda.py`
> Free particle (no electrolyte, eta = 0.0005) with and without the wall soft-sphere repulsion (wall_ss_cut, R = 4, nu = 12,
> eps = 1.044e-12, hc = 3), both restarted from step 30000 (d = 5.04).  (a) d(t); (b) net force on the particle vs d;
> (c) induced charge on the walls (sum) vs d.   Usage: python3 grafico_pared_blanda.py   (STEELE=1: with the Steele wall)
> 2026-11-16 (user's request): panel (c) removed (the wall charge range is still printed); legend of (a) moved off the
> curves (lower left); the formula/parameter text is split into lines so that it is no wider than the two panels.

#### `capacitor/scripts/WALL-DIRICHLET/grafico_particula_libre.py`
> Free charged particle between grounded plates (L = 32, Hann-8, q = +1, eta = 0.0005, initial d = 8 from the lower plate).
> Three runs: no electrolyte; electrolyte (lambda_D = 6) started from scratch (uniform ions, particle free from step 0);
> electrolyte (lambda_D = 6) released from the equilibrated state (restart at step 5000 of the fixed-particle run).
> Reads config.cds* (position/velocity every 100 steps) and proceced_data/particle_force.csv.
> Usage: python3 grafico_particula_libre.py -> ../graficos/particula_libre_eta0.0005.png

#### `capacitor/scripts/WALL-DIRICHLET/graficos_paredes.py`
> Figuras de validacion de las paredes a potencial constante (PETSc, Asta et al.
> 2019). Lee las salidas de Ludwig de A-sin-iones/ y B-debye-huckel/ (con sus
> subcarpetas previo/, que guardan el campo ANTES de corregir psi_electric_field)
> y recalcula las pruebas numericas con verificar_carga_pared.py y
> convergencia_campo.py. Escribe PNG en ../graficos/.
>
> Textos de las figuras en ingles (para presentaciones).
> Uso:  python3 graficos_paredes.py

#### `capacitor/scripts/WALL-DIRICHLET/potencial_kernel_imagenes.py`
> Potencial sobre la linea de la particula (rho = 0) de la carga q con kernel Hann-n entre placas a tierra
> (placas infinitas, sin imagenes laterales): potencial LIBRE del kernel (representacion de una integral 1D) +
> imagenes con la correccion de hexadecapolo.  Placas en z = 0 y z = L, carga en z = d.
>
>  phi_libre(r) = l_B q (2/sqrt(pi)) int_0^inf g(x;t) g(y;t) g(z;t) dt,   [1/r = (2/sqrt(pi)) int_0^inf e^{-r^2 t^2} dt]
>  g(x;t) = int w(x') e^{-t^2 (x-x')^2} dx'  (cerrada con erf compleja, kappa = 2 pi/n):
>    g = (1/n)(sqrt(pi)/(2t)) { erf(t(n/2-x)) + erf(t(n/2+x))
>          + Re[ e^{i kappa x - kappa^2/(4t^2)} ( erf(t(n/2-x) - i kappa/(2t)) + erf(t(n/2+x) + i kappa/(2t)) ) ] }
>  phi_img(z) = l_B q sum_j s_j [ 1/|z-z_j| + (7/4) kappa4 / |z-z_j|^5 ]     (imagenes, sobre el eje)
> Uso: python3 potencial_kernel_imagenes.py -> tabla contra la teoria de Fourier con el kernel real (ell = 200)

#### `capacitor/scripts/WALL-DIRICHLET/verificar_carga_pared.py`
> Capacitor de placas: electrodos solidos (ns capas) en ambos extremos de z,
> fluido en el medio (nf capas), periodico en x,y. Stencil D3Q27, laplaciano
> nabla^2 phi = (2/cs^2) sum_i w_i [phi(r+c_i) - phi(r)] * factor_i.
> Se resuelve SOLO en nodos fluidos (Dirichlet en solidos, columnas eliminadas),
> y se comparan formulas de carga inducida por electrodo.
> Unidades: nabla^2 phi = -4 pi lB rho_el ; tomamos 4 pi lB = 1.


## Proyecto microgel (microgel/scripts/)


### `microgel/scripts/`

#### `microgel/scripts/analizar_campo_electrico.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

#### `microgel/scripts/calc_mui.py`
> Calcula el potencial químico por nodo y especie:
>     μ_i(x) = kT * ln(ρ_i(x)) + z_i * e * ϕ(x)
>
> donde:
>   - ρ_i: densidad de la especie i (del archivo qsi, columnas 0 y 1)
>   - ϕ:   potencial electrostático (del archivo psi, escalar por nodo)
>   - z_i: carga de la especie i (+1 para cationes, -1 para aniones)
>   - e=1, kT: parámetros de entrada
>
> Salidas por paso:
>   - Por pantalla: μ_i mínimo, máximo y delta por especie
>   - mui-NNNNNNNNN.001-001 : mismo formato que qsi (N filas, 2 columnas)
> … (23 líneas más en el archivo)

Uso (del archivo):
    Uso (un paso):
    python calc_mui.py -n 100 -L 32 -kt 0.00001 [-d directorio]
    Uso (serie temporal):
    python calc_mui.py --n-start 100 --n-end 5600 --n-step 100 -L 32 -kt 0.00001

Argumentos (argparse): `-n / --nstep`: Paso de tiempo único; `--n-start`: Paso inicial de la serie; `--n-end`: Paso final de la serie; `--n-step`: Incremento entre pasos de la serie; `-L / --grid-size`: Tamaño de la malla cúbica L×L×L; `-kt / --kt`: Energía térmica kT; `-d / --directory`: Directorio con los archivos de Ludwig (default: .); `--z0`: Carga de la especie 0 (cationes, default: +1); `--z1`: Carga de la especie 1 (aniones, default: -1); `--rho-min`: Densidad mínima para evitar log(0) (default: 1e-20); `--theta`: Calcular θ_k = ρ_k·exp(z_k·ϕ) y guardar theta-NNNNNNNNN.001-001; `--theta-only`: Gráfico de θ_k por especie en ejes separados (implica --theta); `--boltzmann-plot`: Scatter de exp(βz_k eψ) vs ρ_k por nodo y especie; `--ln-rho`: Scatter ln(ρ_k) vs ϕ con recta de Boltzmann (figura independiente); `--csv`: Nombre del archivo CSV de salida (default: mui_series.csv); `--plot`: Nombre del gráfico de salida (default: mui_delta_evolution.png)

#### `microgel/scripts/calcoffset.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-l`: side length

#### `microgel/scripts/calculos.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-nciclos`: number of cycles; `-npaso`: cycles per step; `-o`: output file name; `-i`: input directory (default: current directory); `-outdir`: output directory (default: current directory)

#### `microgel/scripts/calculosvel.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-nciclos`: number of cycles; `-npaso`: cycles per step; `-o`: output file name

#### `microgel/scripts/calculosvelfluid.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-nciclos`: number of cycles; `-ninicio`: Number of cycle to start from; `-npaso`: cycles per step; `-o`: output file name; `--idir`: input directory

#### `microgel/scripts/calculosvelfluidonly.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-nciclos`: number of cycles; `-ninicio`: Number of cycle to start from; `-npaso`: cycles per step; `-o`: output file name

#### `microgel/scripts/coloideacsv.sh`
> Convert config colloid files by batch to .csv deleting vtk files
> Parameters

#### `microgel/scripts/compare_field_theory.py`
> Script para comparar el campo eléctrico simulado con el campo teórico
> de una carga puntual.
>
> Para una carga puntual q en el vacío (o medio uniforme), el campo eléctrico es:
> E = (q / (4πε)) * (r / |r|³)
>
> donde r es el vector desde la carga al punto de observación.
>
> En Ludwig, el sistema está discretizado en una malla y se usa el sistema de
> unidades de lattice Boltzmann.

Argumentos (argparse): `-f / --file`: Archivo psi a leer; `-s / --size`: Tamaño de la malla; `--charge-pos`: Posición de la carga puntual; `--charge`: Magnitud de la carga (default: 1.0); `--epsilon`: Permitividad del medio (default: 1.0); `--kt`: Energía térmica k_B*T (default: 1.0); `--kappa`: Parámetro de Debye κ (inverso de longitud de Debye) para screening iónico. kappa=0 usa Coulomb simple (vacío). kappa>0 usa Debye-Hückel (electrolito). En unidades de lattice. (default: 0.0); `--ionic-strength`: Fuerza iónica en mol/L. Si se proporciona, calcula kappa automáticamente. Requiere conocer la escala de longitud del sistema.; `-m / --mode`: Modo de comparación: line o plane; `--start`: Punto inicial de la línea (obligatorio para mode=line); `--end`: Punto final de la línea (obligatorio para mode=line); `--num-points`: Número de puntos a lo largo de la línea (default: 100); `--show-nodes`: Mostrar solo los nodos de Ludwig con marcadores (sin interpolación); `-p / --plane`: Plano a comparar (obligatorio para mode=plane); `--position`: Posición del plano (índice). Default: centro; `-o / --output`: Archivo de salida (default: mostrar en pantalla); `--external-field`: Componentes del campo externo constante a restar (Ex, Ey, Ez). Ejemplo: --external-field 0.1 0.0 0.0

#### `microgel/scripts/compare_force_coulomb.py`
> Compara la fuerza simulada entre dos colloides contra tres curvas teoricas:
>
>   - Coulomb puro (sin apantallamiento):
>         F(r) = q1 * q2 * e_unit^2 / (4 * pi * epsilon * r^2)
>
>   - Debye-Huckel apantallado, carga puntual (nube ionica linealizada, con
>     el kappa que ya calcula el propio codigo para esta corrida):
>         f(d)    = q1 * q2 * e_unit^2 / (4 * pi * epsilon)
>                   * exp(-kappa * d) * (kappa / d + 1 / d^2)
>
>   - Debye-Huckel apantallado, carga gaussiana (--gaussian-sigma sigma):
>     cada colloide modelado como una nube de carga gaussiana de ancho sigma
> … (41 líneas más en el archivo)

Argumentos (argparse): `base_dir`: Directorio con los subdirectorios pos_X_Y_Z (default: directorio actual); `-n / --step`: Paso de tiempo a comparar (usa el ultimo disponible <= step si no hay coincidencia exacta); `-q / --charge`: Carga q0 de cada colloide (override; si no se da, se lee del log de cada simulacion); `--epsilon`: Permitividad del medio (override; si no se da, se lee del log de cada simulacion); `--kappa`: Parametro de Debye kappa=1/lambda_D (override; si no se da, se toma la mediana de "Kappa calculado" entre los logs de cada simulacion); `--no-dh`: No graficar la curva de Debye-Huckel apantallada; `-L / --box-size`: Tamano de la caja periodica en x (override; si no se da, se lee "System size" del log); `--n-images`: Cantidad de imagenes periodicas a cada lado sumadas en la curva de Debye-Huckel (default: 2; la serie converge rapido por el apantallamiento exponencial); `--no-images`: Graficar Debye-Huckel de una carga aislada, sin sumar las imagenes periodicas de la caja; `--gaussian-sigma`: Si se da, agrega la curva de Debye-Huckel para dos cargas gaussianas de este ancho (sigma, por particula, en unidades de red) en vez de puntuales. Ancho combinado del par = sigma*sqrt(2); `-o / --output`: Archivo de salida del grafico (default: force_vs_distance_step<N>.png dentro del directorio base, junto a las simulaciones); `--kernel`: kernel label shown in the plot title, e.g. "Hann-8"; `--relerr`: bottom panel: relative error (%%) against the Gaussian-Gaussian force, instead of the sim/theory ratio; `--linear`: Usar escala lineal en vez de log-log (default: log-log)

#### `microgel/scripts/compare_force_two_runs.py`
> Two-particle force: corrected and uncorrected runs on the same axes.
>
> Both runs must be a set of pos_* directories at matching separations. The
> theory curves and the reference are shared, so the only difference between
> the two point sets is pm_sr_correction.
>
> Reuses the readers and the screened-Gaussian potential of
> compare_force_coulomb.py, so the numbers match those single-run figures.

Argumentos (argparse): `--no-corr`: base dir of the uncorrected run; `--corr`: base dir of the corrected run; `-n / --step`; `-q / --charge`; `--gaussian-sigma`; `--rcut`; `--images`; `--kernel`; `-o / --output`

#### `microgel/scripts/compare_kernels_force.py`
> Two-particle force error for two kernels, with and without the pm_sr
> correction, on the same axes. Same reference as compare_force_two_runs.py
> (Gaussian-Gaussian Debye-Huckel, 3D periodic image sum); points where the
> target is ~0 by symmetry are dropped.

Argumentos (argparse): `--a-no`; `--a-yes`; `--b-no`; `--b-yes`; `--a-label`; `--b-label`; `--sigma`; `--rcut`; `-o / --output`

#### `microgel/scripts/copy_plots.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-origen`: Carpeta origen; `-destino`: Carpeta destino

#### `microgel/scripts/correr_d8a10.sh`
> Corre las 18 corridas d = 8.2 ... 9.8 (dx = 2 primero, son cortas; despues dx = 1), crear_2part_d8a10.py, 2026-11-11.
> Lanza en paralelo las que entren dejando siempre >= RESERVA MiB libres en la GPU (criterio del usuario, 2026-11-10);
> cada corrida usa <= 1014 MiB (48^3, medido). Espera 60 s entre lanzamientos. No borra nada; log en run.log de cada una.

#### `microgel/scripts/crear_2part_d8a10.py`
> Two fixed charges, intermediate separations d = 8.2, 8.4, ..., 9.8 for dx = 2 and dx = 1 (2026-11-11, user's request:
> look for a fluctuation of the force with the position of the charge inside the cell; the dx = 2 error jumps between
> integer d and dx = 1 does not).
>   dx = 1: copies of SLIDES-LD8-hann08-corr_no-N10000/pos_16.00_24.00_24.00 (d = 8; box 48, Hann-8, 10000 steps);
>           only the x of the moving charge changes: 16 -> 8 + d (the other stays at x = 8; y = z = 24).
>           -> SLIDES-LD8-hann08-corr_no-N10000-d8a10/pos_XX.XX_24.00_24.00
>   dx = 2: copies of SLIDES-LD8-hann04-corr_no-dx2/r_08 (d = 8; box 24, Hann-4, 3000 steps); the moving charge
>           8.25 -> 4.25 + d/2 (same physical position as dx = 1; the other stays at 4.25; y = z = 12.25).
>           -> SLIDES-LD8-hann04-corr_no-dx2-d8a10/r_XX.X
>   Inside one coarse cell: dx = 2 moves the charge by 0.1 lattice units per point, dx = 1 by 0.2.
> Executable: the CURRENT src/Ludwig.exe (md5 45dc94bc; with hydrodynamics on it gives results byte-identical to 4f0c20e7,
> checked 2026-11-11 with TEST-HIDRO-dx2/hidro_si), as the user chose for d = 4.1 ... 4.9. The d = 8 and d = 10 runs
> … (2 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 crear_2part_d8a10.py   and then   ./correr_d8a10.sh > bateria_d8a10.log

#### `microgel/scripts/crear_2part_dx1_intermedios.py`
> Two fixed charges, dx = 1, intermediate separations d = 4.1, 4.2, ..., 4.9 (2026-11-11, user's request: points every 0.1
> between d = 4 and 5). Copies of SLIDES-LD8-hann08-corr_no-N10000/pos_12.00_24.00_24.00 (box 48^3, lambda_D = 8, Hann-8,
> 10000 steps, freq_config 1000) where ONLY the x of the moving charge changes: 12 -> 8 + d (the other stays at x = 8;
> y = z = 24). With d non-integer the moving charge is between nodes in x.
> Executable (user's choice): the CURRENT src/Ludwig.exe (md5 4f0c20e7, 2026-11-04, with the halo fix of 2026-09-26), NOT
> the 1c30df4b of the d = 4 and d = 5 runs of SLIDES-LD8-hann08-corr_no-N10000 (those are not redone).
> Creates SLIDES-LD8-hann08-corr_no-N10000-d4a5/pos_XX.XX_24.00_24.00; an existing directory is not touched.
> Usage: python3 crear_2part_dx1_intermedios.py   and then   cd SLIDES-LD8-hann08-corr_no-N10000-d4a5 && ./correr.sh 3 > bateria.log

#### `microgel/scripts/crear_2part_dx2.py`
> Two fixed charges in the periodic box with dx = 2 (2026-11-09, user's request): same physical system as
> SLIDES-LD8-hann08-corr_no-N10000 (dx = 1, box 48, Hann-8, lambda_D = 8, q1 = q2 = +1), refined by r = 1/2 with the
> recipe of capacitor/local/WALL-DIRICHLET/crear_refinado.py:
>   - lengths x r: box 48 -> 24, Hann-8 -> Hann-4 (physical width 8)
>   - electrokinetics_epsilon / r (= 2e4): q/(eps kT) in lattice units scales with r, psi(physical x) unchanged
>   - electrokinetics_init_rho_el / r^3 (= 6.17504e-3): with eps/r gives lambda_D = 8 r = 4 nodes (8 physical)
>   - D per step unchanged -> diffusion time in steps scales as r^2: N = 10000 r^2 = 2500, rounded up to 3000
>   - freq_config 1000 r^2 = 250 (steady state: psi(3000) - psi(2750))
> Positions: the same PHYSICAL positions as dx = 1 (user's choice). Ludwig node i is at coordinate i, the box spans
> [0.5, N + 0.5], so the physical coordinate is X - 0.5 and X_r = 0.5 + r (X_1 - 0.5): the fixed charge at X_1 = 8 ->
> 4.25, the other at X_1 = 8 + d -> 4.25 + d/2, y = z = 24 -> 12.25. No charge sits on a node of the coarse lattice
> (with dx = 1 both are on nodes); psi_electroneutral therefore removes no "near" nodes here (it removes the 2 occupied
> … (4 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 crear_2part_dx2.py   and then   cd SLIDES-LD8-hann04-corr_no-dx2 && ./correr.sh 4 > bateria.log

#### `microgel/scripts/crear_2part_dx2_r4_shift.py`
> Two fixed charges, dx = 2, FIXED separation r = 4 (physical), both charges shifted together inside the cell
> (2026-11-11, user's request: force vs intra-cell position of the pair).
>   Shift delta (physical) = 0, 0.2, ..., 2.0 (one full cell of the dx = 2 lattice = 1 node; delta = 0 and 2.0 must be
>   equivalent by translation invariance: control), along x only (user's choice). In lattice units both charges move by
>   delta/2: charge 1 6.25 -> 6.25 + delta/2, charge 2 4.25 -> 4.25 + delta/2; y = z = 12.25 unchanged.
> Copies of SLIDES-LD8-hann04-corr_no-dx2/r_04 (box 24^3, Hann-4, lambda_D = 8 physical, 3000 steps); only the two charge
> positions change. Executable: the CURRENT src/Ludwig.exe (45dc94bc); delta = 0 is redone with it (r_04 used 1c30df4b).
> Creates SLIDES-LD8-hann04-corr_no-dx2-r4-shift/d_X.X; existing directories are not touched.
> Usage: python3 crear_2part_dx2_r4_shift.py   and then
>        cd SLIDES-LD8-hann04-corr_no-dx2-r4-shift && ./correr.sh 4 > bateria.log

#### `microgel/scripts/crear_2part_dxhalf.py`
> Two fixed charges in the periodic box with dx = 1/2 (2026-11-09, user's request): ONLY separation d = 4, 20000 steps,
> psi/efield every 2000 steps to see when the steady state is reached (user's choice). Same physical system as
> SLIDES-LD8-hann08-corr_no-N10000 (dx = 1) and SLIDES-LD8-hann04-corr_no-dx2 (dx = 2), refined by r = 2 with the recipe of
> crear_2part_dx2.py:
>   - lengths x r: box 48 -> 96, Hann-8 -> Hann-16 (physical width 8)
>   - electrokinetics_epsilon / r (= 5e3), electrokinetics_init_rho_el / r^3 (= 9.6485e-5): lambda_D = 16 nodes = 8 physical
>   - D per step unchanged: the dx = 1 runs needed ~9000-10000 steps -> ~36000-40000 here; 20000 = t 5000 (user's choice)
> Positions: same PHYSICAL positions, X_r = 0.5 + r (X_1 - 0.5): x = 8 -> 15.5, x = 8 + d -> 15.5 + 2d, y = z = 24 -> 47.5
> (between nodes in the three coordinates; with dx = 1 both charges are on nodes).
> Executable (user's choice): the CURRENT src/Ludwig.exe (md5 4f0c20e7..., 2026-11-04, with the halo double-count fix of
> 2026-09-26), NOT the one of the dx = 1 and dx = 2 runs (1c30df4b..., 2026-09-15).
> Creates SLIDES-LD8-hann16-corr_no-dxhalf/r_04; an existing directory is not touched.
> … (3 líneas más en el archivo)

Uso (del archivo):
    2026-11-10: optional arguments d and N (python3 crear_2part_dxhalf.py 5 12000 -> r_05 with 12000 steps); without
    Usage: python3 crear_2part_dxhalf.py [d N]   and then   cd SLIDES-LD8-hann16-corr_no-dxhalf && ./correr.sh 1 > bateria.log

Lee argumentos de línea de comandos: sys.argv[1, 2]

#### `microgel/scripts/crear_continuacion_dxhalf.py`
> Continuation of the dx = 1/2 two-charge runs SLIDES-LD8-hann16-corr_no-dxhalf/r_05, r_06, r_07 up to step 20000
> (2026-11-10, user's request: continue the runs already done, no new runs). Recipe of
> capacitor/local/WALL-DIRICHLET/crear_continuacion.py (Ludwig restarts from dist/psi/qsi/rho/vel of step N0 and
> config.cds000N0, "Re-starting simulation at step N0").
>   r_05: from step 12000 (end of its run), r_07: from 12000, r_06: from 10000 (its run died at 10775 by GPU out of memory;
>   10000 is its last dump).
> The originals are NOT touched: SLIDES-LD8-hann16-corr_no-dxhalf-cont/r_XX gets copies of input, config.cds.init,
> Ludwig.exe (same exe, 4f0c20e7), the *-metadata files and the files of step N0 only (not the earlier dumps), and the
> proceced_data/*.csv truncated to step <= N0 (Ludwig appends to existing csv files; r_06 had rows up to 10775).
> Input: N_start N0, N_cycles 20000 - N0 (freq_config 2000 unchanged).
> correr_cont.sh does not delete anything; log in run_cont.log; [OK] requires "finished normally" and the restart line.
> Usage: python3 crear_continuacion_dxhalf.py   and then   cd SLIDES-LD8-hann16-corr_no-dxhalf-cont && ./correr_cont.sh 3 > bateria.log

#### `microgel/scripts/csv_watch.sh`
> Convierte cada config.cds recien volcado a un colloids-*.csv para Paraview,
> a medida que la simulacion corre -- el mismo paso que runbg.sh hace via
> runplot.sh/coloideacsv.sh, pero sin el resto de ese flujo (que asume un
> input armado por runbg.sh). Pensado para correr en paralelo a un run.sh
> que ya lanzo Ludwig.exe por su cuenta.
> Uso:
> csv_watch.sh RUNDIR [--pid PID] [--poll SEGUNDOS]
> Sin --pid corre hasta que no aparezcan volcados nuevos durante --poll*3
> segundos seguidos (util para "engancharse" a una corrida ya lanzada, o
> para pasar sobre una que ya termino y dejo volcados sin convertir).

#### `microgel/scripts/delete_files.sh`
> Delete unnecesary files.
> Leave last file so that restart is possible.

#### `microgel/scripts/delete_files_all_leave_image copy.sh`
> Delete unnecesary files.
> Leave last file so that restart is possible.

#### `microgel/scripts/delete_files_all_leave_image.sh`
> Delete unnecesary files.
> Leave last file so that restart is possible.

#### `microgel/scripts/delete_files_dist_rho.sh`
> Delete unnecesary files.
> Leave last file so that restart is possible.

#### `microgel/scripts/energia_2part.py`
> Free energy of the two-charge system vs time and chemical-potential residual (2026-11-10, user's request: steps 0-2).
>
> Runs: separation r = 4 (physical) in the three lattices of the same physical system (box 48, lambda_D = 8, Hann of
> physical width 8, cos^2, q1 = q2 = +1):
>   dx = 2  : SLIDES-LD8-hann04-corr_no-dx2/r_04                 (24^3, Hann-4,  charges at x = 4.25, 6.25;  y = z = 12.25)
>   dx = 1  : SLIDES-LD8-hann08-corr_no-N10000/pos_12.00_24.00_24.00 (48^3, Hann-8, x = 8, 12;  y = z = 24)
>   dx = 1/2: SLIDES-LD8-hann16-corr_no-dxhalf/r_04              (96^3, Hann-16, x = 15.5, 23.5;  y = z = 47.5)
> Every psi/qsi dump (psi in e psi/kT; qsi = rho_+, rho_- of the fluid only, per node; valencies Z = +1, -1).
>
> Step 0 (checks): N_a = sum_n rho_a,n per dump (must be constant: NP is in flux form); ion charge Q_ion = N_+ - N_-;
>   fixed charge Q_fix = sum_n rho_fix,n (2 if the lattice Hann weights sum to 1); total charge Q_ion + Q_fix.
> Step 1: free energy in kT, the functional of fe_electro.c (Rotenberg et al. 2010), psi already in units of kT/e:
> … (11 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 energia_2part.py   (from microgel/local/scripts) -> table printed and ../FUERZA-2PART-MEAS/ENERGIA_r04.csv

#### `microgel/scripts/error_capas_2part.py`
> Error of psi and E at the lattice nodes, two fixed charges in the periodic box (2026-11-09, user's request; same metrics as
> capacitor/local/WALL-DIRICHLET/scripts/grafico_C_error_capas_desplazada.py), two lattices and the observed order.
>
> Runs (same physical system, box 48, lambda_D = 8, Hann of physical width 8, q1 = q2 = +1, separation r = 1 ... 16):
>   dx = 1: SLIDES-LD8-hann08-corr_no-N10000/pos_X.00_24.00_24.00 (box 48^3, Hann-8, charges at x = 8 and X = 8 + r,
>           y = z = 24), 10000 steps; stationarity psi(10000) - psi(9000).
>   dx = 2: SLIDES-LD8-hann04-corr_no-dx2/r_DD (box 24^3, Hann-4, charges at x = 4.25 and 4.25 + r/2, y = z = 12.25: the
>           same physical positions, crear_2part_dx2.py), 3000 steps; stationarity psi(3000) - psi(2750).
>   (The first version used SLIDES-LD8-hann08-corr_no, 3000 steps: not stationary; and kappa of its run.log, 0 at r = 16.)
> Units: psi dump in e psi/kT; efield dump / kT in lattice units, divided by dx -> per physical length. Distances, shells and
> norms in physical units (dx = 1 lattice).
>
> … (17 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 error_capas_2part.py [w=1]   (from microgel/local/scripts)

#### `microgel/scripts/error_fuerza_2part.py`
> Force between two fixed charges in the periodic box vs separation, dx = 2 and 1 (2026-11-09, user's request).
>
> Runs (as error_capas_2part.py): dx = 1 SLIDES-LD8-hann08-corr_no-N10000 (box 48, Hann-8), dx = 2 SLIDES-LD8-hann04-corr_no-dx2
> (box 24, Hann-4, same physical positions); q1 = q2 = +1, lambda_D = 8; charge 2 fixed, charge 1 at separation r along +x.
> Ludwig force on charge 1 (physical units, kT per physical length times kT): F_L = kT Esub_x / dx, Esub_x of particle 1 at
> the last step (proceced_data/particle_Esub.csv). Esub, not particle_force.csv: in the dx = 2 runs (old exe 1c30df4b) the
> force of particle 2 (x = 4.25, edge colloid cell) is doubled by the halo bug fixed 2026-09-26; particle 1 and Esub are right.
>
> Reference (user's choice: option 2, the same as for psi and E): linear Debye-Hueckel, periodic box L = 48, both charges spread
> with the continuous Hann window of width 8; force = kernel-weighted field at charge 1, F_ref = kT q int W(r - r1) E_DH,x(r) d^3r.
> With phi(r) = sum_k A(k) sum_p e^{i k.(r - r_p)}, A = (4 pi l_B/V) What3(k)/(k^2 + kappa^2), What3 = What(kx) What(ky) What(kz),
> and int W(r - r1) e^{i k.r} d^3r = What3(k) e^{i k.r1}:
> … (8 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 error_fuerza_2part.py   -> ../FUERZA-2PART-MEAS/FUERZA_LD08.0_Hann-8_corr-no_dx2-1.png (+ table printed)

#### `microgel/scripts/extraer_posicion.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-nciclos`: number of cycles; `-ninicio`: Number of cycle to start from; `-npaso`: cycles per step; `-o`: output file name

#### `microgel/scripts/grafico_fuerza_shift_r4_dx2.py`
> Force between two fixed charges at FIXED separation r = 4 vs the shift delta of the pair inside the cell, dx = 2
> (2026-11-11, user's request; runs of crear_2part_dx2_r4_shift.py: SLIDES-LD8-hann04-corr_no-dx2-r4-shift/d_X.X,
> both charges moved along x by delta (physical) = 0 ... 2.0, one full cell; exe 45dc94bc; 3000 steps).
> Force on each charge, physical units: F = kT Esub_x / dx (dx = 2), last step (proceced_data/particle_Esub.csv;
> index 1 = charge at larger x, index 2 = the other). Plotted: F_x on charge 1 and -F_x on charge 2 (equal for an
> isolated pair), with the reference F_ref(4) of error_fuerza_2part.py (periodic DH, continuous Hann of width 8; it does not
> depend on delta). Panels: force, relative error (F - F_ref)/F_ref in %, and (2026-11-11) variation about the mean of one period
> (F - <F>)/<F> in %, <F> over delta = 0 ... 1.8.
> Printed: table, and the control delta = 0 vs delta = 2.0 (same configuration translated by one node).
> Usage: python3 grafico_fuerza_shift_r4_dx2.py   -> ../FUERZA-2PART-MEAS/FUERZA_shift_r4_dx2.png

#### `microgel/scripts/grafico_psi_E_2part.py`
> Potential and electric field of two fixed charges (r = 4 physical), same figure as
> TEST-HIDRO-dx1/scripts/grafico_psi_E_hidro_no.py but for any lattice (2026-11-11, user's request: Hann-16 and Hann-8 WITH
> hydrodynamics, to compare with the Hann-8 run without hydrodynamics).
>
> Cases (python3 grafico_psi_E_2part.py <case>):
>   hann8_hidro : SLIDES-LD8-hann08-corr_no-N10000/pos_12.00_24.00_24.00, dx = 1, box 48^3, Hann-8, charges x = 8, 12,
>                 y = z = 24 (on nodes), step 10000, exe 1c30df4b (old), hydrodynamics on.
>   hann16_hidro: SLIDES-LD8-hann16-corr_no-dxhalf/r_04, dx = 1/2, box 96^3, Hann-16, charges x = 15.5, 23.5,
>                 y = z = 47.5 (BETWEEN nodes), step 20000, exe 4f0c20e7, hydrodynamics on.
> Everything is plotted in PHYSICAL units (dx = 1 lattice): x_phys = 0.5 + (X - 0.5) dx (X = Ludwig coordinate),
> E per physical length = (efield / kT) / dx; psi in e psi/kT (lattice independent).
> Line: the node line y = z = YL closest to the charge axis. For dx = 1 it is the axis itself; for dx = 1/2 the axis
> … (11 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 grafico_psi_E_2part.py hann16_hidro   -> ../FUERZA-2PART-MEAS/PSI_E_<case>.png

Lee argumentos de línea de comandos: sys.argv[1]

#### `microgel/scripts/graficos.sh`
> Script to collect all velocidades.png files into a single directory

#### `microgel/scripts/graficos_neutro_vs_cargado.py`
> Compara el microgel neutro en equilibrio con el mismo microgel cargado en la
> superficie, sin sal, tambien en equilibrio (misma estructura de red de
> partida: MGEL-RED-neutro/config.cds00060000 == MGEL-RED-SURF-NOSALT
> config.cds.init, antes de aplicar la carga).
>
> (a) Trayectoria de Rg: relajacion neutra (0->60000) y relajacion bajo carga
>     (0->500000), cada una con su propia asintota extrapolada.
> (b) Radio de cada monomero al centro de masa, antes (=neutro) y despues
>     (=cargado en equilibrio), separado en capa cargada (65 monomeros
>     externos) y nucleo sin carga (93 monomeros), para ver DONDE hincha.
>
> Escribe graficos/neutro_vs_cargado.png

#### `microgel/scripts/lattice_error_gauss_vs_hann.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

#### `microgel/scripts/lattice_error_gauss_width.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

#### `microgel/scripts/ludwig_redo.sh`
> Compile ludwig and copy program to current dir

#### `microgel/scripts/make_dimer.py`
> Genera ficheros de configuracion de colonias (config.cds.init.001-001) para un
> DIMERO: dos particulas subgrid cargadas, enlazadas por un bond harmonico.
>
> Escribe el estado en el orden EXACTO de colloid_state_write_ascii (src/colloid.c):
> un valor por linea, primero los enteros (32) y luego los doubles.
>
> Uso:
>     python3 make_dimer.py
> Produce:
>     config.cds.init.001-001.dimer_pp   (cargas +1 / +1, repulsivo)
>     config.cds.init.001-001.dimer_pm   (cargas +1 / -1, atractivo)

#### `microgel/scripts/make_plots.sh`
> Run multiple plots to compare teoretical Electric field to simulations
> ./make_plots.sh -d 3 -l 32 -e 10000 -k 0.200761489786388 -o psi_field -p 2000 -t 0.00001

#### `microgel/scripts/mgel_equilibrio.py`
> Estado de equilibrio de una corrida de microgel a partir de su serie de Rg.
>
> Ajusta la pendiente dRg/dpaso de los ultimos volcados a una exponencial
>   s(t) = s0 exp(-(t - t0)/tau)
> y reporta:
>   - Rg actual y pendiente (por 1000 pasos, promediada en 10000 pasos)
>   - tau y el Rg asintotico proyectado  Rg_inf = Rg + s * tau
>   - si se cumple el criterio de equilibrio y, si no, en que paso se cumpliria
>
> Criterio (por defecto): pendiente < 2e-4 por 1000 pasos, lo que con tau ~1.3e5
> deja menos de ~0.4% de hinchamiento por delante.
>
> … (2 líneas más en el archivo)

Uso (del archivo):
    Uso:

Argumentos (argparse): `dir`; `--umbral`: pendiente dRg/1000 pasos que se considera equilibrio; `--ventana`: pasos finales usados en el ajuste; `--cada`: usar solo volcados multiplos de este paso; `-o / --out`

#### `microgel/scripts/mgel_serie.py`
> Sigue la evolucion estructural de un microgel a lo largo de una corrida.
>
> Lee los config.cds volcados por Ludwig y, para cada paso, mide:
>   Rg     radio de giro (RMS de la distancia al centro de masa)
>   diam   2 * distancia maxima al centro de masa
>   bond   longitud media de enlace (y su minimo/maximo)
>   hmin   separacion minima entre monomeros no ligados
>   vmax   modulo de velocidad maxima
>   z      enlaces por monomero (constante: la topologia no cambia)
>
> El equilibrio se juzga por dRg/dpaso: mientras crezca de forma sostenida el
> microgel se sigue hinchando. Con la estructura de cadena ramificada (z=1.92)
> … (5 líneas más en el archivo)

Uso (del archivo):
    Uso:

Argumentos (argparse): `dirs`: directorio(s) de corrida; `-o / --out`: nombre del json dentro de cada directorio; `--every`: procesar 1 de cada N volcados

#### `microgel/scripts/microgel_redo.sh`
> Compile microgel_poly_cross_density from ../util and copy program to current dir

#### `microgel/scripts/microgel_run.sh`
> Run microgel_poly_cross_density

#### `microgel/scripts/mu_capas_2part.py`
> Chemical-potential residual around the charges, two fixed charges r = 4, three lattices (2026-11-10, user's request).
>
> Same runs and definitions as energia_2part.py: mu_a/kT = ln rho_a + Z_a psi (qsi = fluid only, psi in kT/e).
> delta mu_a,n = mu_a,n - <mu_a>, <.> = plain mean over all nodes (choice of this script; the residual max - min of
> energia_2part.py does not depend on it).
> Per shell of physical distance d to the NEAREST charge (minimum image, width 1 physical, as error_capas_2part.py):
> signed mean of delta mu, RMS and max |delta mu| (kT), for each species.
> Dumps: physical time t = 5000 (steps 1250, 5000, 20000 for dx = 2, 1, 1/2: same time in the three lattices) and the last
> dump of each run (dx = 2: t = 12000, dx = 1: t = 10000; dx = 1/2: t = 5000, the same).
> Usage: python3 mu_capas_2part.py   (from microgel/local/scripts) -> tables printed

#### `microgel/scripts/organize_by_parameter.sh`
> Script para organizar carpetas de simulaciones por parámetro
> Uso: ./organize_by_parameter.sh <directorio_base> <parametro> [prefijo]
> Ejemplo: ./organize_by_parameter.sh /path/to/simulations rhoel PB_Sn

#### `microgel/scripts/organize_by_parameter_dryrun.sh`
> Script para VISUALIZAR cómo se organizarían las carpetas de simulaciones por parámetro
> (DRY RUN - no mueve nada, solo muestra qué haría)
> Uso: ./organize_by_parameter_dryrun.sh <directorio_base> <parametro> [prefijo]

#### `microgel/scripts/plot.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-i`: input file; `-outdir`: output directory (default: current directory); `-o`: output file name (if not set, show plot interactively)

#### `microgel/scripts/plot_Esub_vs_pos.py`
> Script para graficar las componentes del campo E_sub sobre la partícula en
> función de su posición, recorriendo las subcarpetas de un directorio padre.
>
> Cada subcarpeta corresponde a una simulación con la partícula en una posición
> distinta (codificada en el nombre como pos_X_Y_Z). De cada una se lee el
> archivo proceced_data/particle_Esub.csv (separador ';', decimal ',') y se
> extraen Esub_X, Esub_Y, Esub_Z (y Emod) para el step indicado.
>
> Uso:
>     python plot_Esub_vs_pos.py --step 1000 --output Esub_vs_pos.png
>     python plot_Esub_vs_pos.py --step 500 -o salida.png --parent-dir /ruta/al/directorio

Argumentos (argparse): `--step / -s`: Número de step a extraer de particle_Esub.csv; `--output / -o`: Nombre del archivo de salida (png); `--parent-dir / -d`; `--axis`: Eje de posición para el plot (default: el que más varía); `--show-emod`: Incluir también |E_sub| (Emod) en el gráfico

#### `microgel/scripts/plot_campo_electrico_rhoel.py`
> Script para graficar campo_electrico_medio.csv de múltiples simulaciones
> variando el parámetro rhoel

Lee argumentos de línea de comandos: sys.argv[1, 2]

#### `microgel/scripts/plot_charge_distribution.py`
> Script para graficar la distribución de carga desde archivos qsi de Ludwig.
>
> Los archivos qsi contienen la densidad de carga en cada punto de la malla.
> Formato típico (2 columnas):
>   - Columna 0: Densidad de carga positiva ρ₊ (cationes, carga +1)
>   - Columna 1: Densidad de carga negativa ρ₋ (aniones, carga -1)
>
> IMPORTANTE: Las densidades son siempre ≥ 0 (representan concentraciones).
> La carga neta en cada nodo es: ρ_net = ρ₊ - ρ₋
>
> Modos de visualización:
> 1. Plano 2D: Muestra la distribución de carga en un plano específico (XY, XZ, o YZ)
> … (11 líneas más en el archivo)

Argumentos (argparse): `-f / --file`: Archivo qsi a leer; `-s / --size`: Tamaño de la malla (nx ny nz); `-m / --mode`: Modo: plane (2D), line (1D), plane3d (superficie 3D), stats (solo estadísticas); `-p / --plane`: Plano a graficar (obligatorio para mode=plane/plane3d); `-pos / --position`: Posición del plano (índice). Por defecto: centro; `--start`: Punto inicial de la línea (obligatorio para mode=line); `--end`: Punto final de la línea (obligatorio para mode=line); `--num-points`: Número de puntos a lo largo de la línea (por defecto: 100); `-c / --component`: Componente a graficar: net (carga neta), species0/1 (especie individual), total (densidad total), all (todas en modo line). Default: net; `--colormap`: Mapa de colores para modo plane3d (default: RdBu_r); `--elevation`: Ángulo de elevación de la vista 3D (default: 30); `--azimuth`: Ángulo azimutal de la vista 3D (default: -60); `--no-interp`: Usar pcolormesh (sin interpolación, valores exactos por nodo) en lugar de contourf; `--num-species`: Número de especies iónicas (default: 2); `-o / --output`: Archivo de salida (por defecto: mostrar en pantalla)

#### `microgel/scripts/plot_charge_vs_r.py`
> Radial profile of the ionic (fluid) charge around one subgrid charge, with
> and without the pm_sr short-range correction, against theory.
>
> The qsi dump holds the two ion densities (rho0, rho1); the net node charge is
> rho0 - rho1, in the same units as the particle charge (verified: the sum over
> the box is exactly -q). In linearised Poisson-Boltzmann the induced charge is
> rho_ion = -eps * kappa^2 * phi, so the same screened potential used for the
> field gives the charge with no extra calibration.

Argumentos (argparse): `--no-corr`: qsi dump of the uncorrected run; `--corr`: qsi dump of the corrected run; `-n / --ntotal`; `--pos`; `-q / --charge`; `--epsilon`; `--ld`; `--sigma`; `--rcut`; `--rmax`; `--images`; `--kernel`; `-o / --output`

#### `microgel/scripts/plot_componentes_ewald.py`
> Script para plotear por separado las componentes Real y Recíproca del campo eléctrico
> calculadas por Ewald en Ludwig.
>
> Este script ayuda a diagnosticar el problema de desconexión entre las dos partes.
>
> Para usar este script, necesitamos modificar temporalmente ewald_charge.c para
> que genere dos archivos separados:
>   - efield_real: Componente de espacio real
>   - efield_fourier: Componente de espacio recíproco
>
> Autor: Script para análisis de Ewald

Argumentos (argparse): `-n / --nstep`: Número de paso temporal (ej: 1, 1000, etc); `-d / --directory`: Directorio donde están los archivos (default: directorio actual); `-L / --grid-size`: Tamaño de la malla (asume cubo LxLxL); `-q / --charge`: Carga de la partícula; `-eps / --epsilon`: Permitividad relativa ε_r; `-kt / --kt`: Energía térmica kT (default: 1.0); `--kappa`: Parámetro de Debye κ = 1/λ_D (si no se especifica, se calcula de rho_el); `--rho-el`: Densidad electrolítica ρ_el para calcular κ (si no se da kappa); `--alpha`: Parámetro de splitting de Ewald α (default: 0.5); `--rc`: Radio de corte rc (default: 5.0); `-o / --output`: Archivo de salida (default: ewald_components.png); `--rmin`: Radio mínimo (default: 0.1); `--rmax`: Radio máximo (default: L/2); `--max-points`: Máximo número de puntos (default: 10000); `--linear`: Usar escala lineal (default: logarítmica); `--log-log`: Usar doble eje logarítmico log-log (default: semi-log); `--log-log-error`: Usar doble eje logarítmico en el panel de error relativo; `--xmin`: Valor mínimo del eje X en ambos paneles (default: 0 o rmin en log-log); `--ymin`: Valor mínimo del eje Y (aplica a ax1 en log/semilog, y a ax2 en log-log-error); `--ymax`: Valor máximo del eje Y en el panel 1 (default: automático); `--show-ewald-components`: Mostrar componentes real y Fourier de Ewald por separado; `--show-ewald-total`: Mostrar la suma total de Ewald numérica (real + Fourier); `--show-ewald-theory-total`: Mostrar la suma teórica total de Ewald (real + recíproco teóricos); `--show-dh`: Mostrar la curva teórica de Debye-Hückel; `--show-dh-periodic`: Mostrar Debye-Hückel con imágenes periódicas; `--dh-periodic-shells`: Número de capas de imágenes periódicas (default: 1 = 26 vecinos a distancia ~L); `--direction`: Dirección para filtrar datos y evaluar DH periódico. Ejemplos: 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz', '1,0,0', '1,1,0'. Sin este argumento se usan todos los puntos (isotrópico).; `--angle-tol`: Tolerancia angular en grados para filtrar puntos por dirección (default: 15); `--pointwise-error`: Calcula el error de cada nodo contra el DH periódico evaluado en el mismo punto 3D exacto (no en r·d̂). Requiere --show-dh-periodic y --direction.; `--no-error-inf`: Ocultar la serie de error vs DH infinito (púrpura) en el panel 2; `--no-error-dir`: Ocultar la serie de error vs DH periódico direccional (naranja) en el panel 2; `--particle-radius`: Radio de la partícula a (lu). Activa curva DH radio finito en panel 1 y error vs DH(a) en panel 2 (default: 0 = desactivado); `--only-error-dh-per-rad`: Mostrar en el panel de error SOLO la curva de error relativo vs DH periódico con partícula de radio a. Requiere --particle-radius > 0 y --show-dh-periodic-rad. Suprime las demás curvas de error.; `--show-dh-periodic-rad`: Mostrar DH periódico con radio finito en panel 1, independientemente de --show-dh-periodic. Requiere --particle-radius > 0.; `--gaussian-source`: Usar distribución gaussiana ρ(r)=Q/(π^(3/2)σ³)exp(-r²/σ²) como fuente en la teoría DH. Activa curva |E| gaussiana-DH.; `--sigma`: Anchura σ de la distribución gaussiana de carga (lu). Solo se usa con --gaussian-source (default: 1.0). Si se especifican --sigma-particle y --sigma-fluid, este valor se ignora y se usa σ_eff = √(σ_p² + σ_f²).; `--sigma-particle`: Anchura σ_p de las partículas gaussianas (lu). Combinado con --sigma-fluid define la interacción par σ_eff = √(σ_p² + σ_f²). Override de --sigma.; `--sigma-fluid`: Anchura σ_f de los nodos de fluido gaussianos (lu). Combinado con --sigma-particle define la interacción par σ_eff = √(σ_p² + σ_f²). Override de --sigma.; `--gaussian-error`: Mostrar el error relativo vs la teoría gaussiana-DH (infinita) en el panel 2. Requiere --gaussian-source.; `--gaussian-error-periodic`: Mostrar el error relativo vs la teoría gaussiana-DH periódica en el panel 2. Requiere --gaussian-source y --show-dh-periodic.; `--no-gaussian-dh-inf`: No graficar la curva gaussiana-DH sistema infinito (no periódica).; `--no-ewald-params / --pm`: No mostrar los parámetros de Ewald (α, rc) en el título ni las líneas verticales de rc. Usar cuando la corrida es PM (sin Ewald).; `--csv`: Exportar datos a CSV (nombre base sin extensión; se añade -n_NNNNNN.csv); `--csv-points`: Exportar a CSV los mismos puntos que aparecen en el gráfico (con filtro de dirección), incluyendo error vs DH∞ y DH periódico. Nombre base sin extensión; se añade -n_NNNNNN.csv

#### `microgel/scripts/plot_dimero_perfil.py`
> Perfiles de campo electrico (E), potencial (psi) y carga (rho/qsi) para un DIMERO.
>
> - Perfil 1D: a lo largo de la recta que une las 2 particulas.
> - Planos 2D: los dos planos que contienen esa recta (el que tiene la normal en
>   cada una de las dos direcciones transversales).
>
> Formato de datos Ludwig (ver plot_electric_field.py como base):
>   - psi / rho / qsi : 1 columna, nx*ny*nz filas, orden z-major -> reshape (nx,ny,nz)
>   - efield          : 3 columnas (Ex,Ey,Ez)
>   - Nodo entero i (0-based) esta en la posicion fisica i + 0.5 (Lmin=0.5).
>   - Posiciones de particulas: colloid_data/colloids-XXXXXXXX.csv (id,x,y,z,...)
>
> … (3 líneas más en el archivo)

Uso (del archivo):
    Uso tipico:
    python3 plot_dimero_perfil.py -d <run_dir> -step 2000
    python3 plot_dimero_perfil.py -d <run_dir> -step 2000 -L 32 -o salida_

Argumentos (argparse): `-d / --dir`: directorio de la corrida; `-step`: time-step a graficar; `-L`: lado del grid; `-o / --out`: prefijo de salida; `--no-planes`; `--no-theory`: no superponer la curva teorica Debye-Hueckel; `-q / --charge`: carga de cada particula (default: del nombre de carpeta); `--full-box`: la recta cubre todo el ancho L de la caja (periodica), no solo el tramo entre particulas +/- margen; `--extend`: margen a cada lado de las particulas (si no --full-box); `--gaussian-sigma`: usar teoria D-H con cargas GAUSSIANAS de ancho SIGMA (finita en los centros) en vez de cargas puntuales

#### `microgel/scripts/plot_efield.py`
> Script para graficar el campo eléctrico sobre una línea y sobre un plano
> leyendo directamente los archivos efield de Ludwig.
>
> Soporta archivos: efield, efield_real, efield_fourier
>
> Ejemplos de uso:
>
>   # Campo total a lo largo del eje X en y=16.5, z=16.5
>   python plot_efield.py -d DIR -n 10000 -L 32 --mode line       --start 0.5 16.5 16.5 --end 31.5 16.5 16.5
>
>   # Campo en plano XY (z fijo) coloreado por magnitud
>   python plot_efield.py -d DIR -n 10000 -L 32 --mode plane       --plane xy --position 16
> … (9 líneas más en el archivo)

Uso (del archivo):
    python plot_efield.py -d DIR -n 10000 -L 32 --mode line       --start 0.5 16.5 16.5 --end 31.5 16.5 16.5       --charge-pos 16.1 16.5 16.5 --charge 1.0 --kappa 0.1 --epsilon 1e4 --kt 1e-5
    python plot_efield.py -d DIR -n 10000 -L 32 --mode plane       --plane xz --position 16 --quiver
    python plot_efield.py -d DIR -n 10000 -L 32 --mode line       --start 0.5 16.5 16.5 --end 31.5 16.5 16.5 --field-type real

Argumentos (argparse): `-d / --directory`: Directorio con los archivos efield (default: .); `-n / --nstep`: Número de paso temporal; `-L / --grid-size`: Tamaño de malla: un valor (LxLxL) o tres (NX NY NZ); `--field-type`: Tipo de archivo efield a leer (default: total); `--colloid-file`: Archivo colloids-*.csv para leer posición del coloide; `--mode`: Modo de visualización; `--start`: Punto inicial de la línea; `--end`: Punto final de la línea; `--num-points`: Puntos de interpolación en la línea (default: 200); `--show-nodes`: Mostrar solo nodos de la malla (sin interpolación); `--plane`: Plano de corte; `--position`: Coordenada física del plano de corte (ej: 16.0). Se selecciona el nodo más cercano. (default: centro); `--quiver`: Superponer flechas de campo en el plano; `--quiver-step`: Paso de submuestreo para quiver (default: 2); `--charge-pos`: Posición de la carga para comparación teórica; `--charge`: Carga q (default: 1.0); `--kappa`: Parámetro de Debye κ (default: 0 = Coulomb); `--epsilon`: Permitividad ε (informativo, default: 1.0); `--kt`: Energía térmica kT (informativo, default: 1.0); `-o / --output`: Archivo de salida (default: mostrar en pantalla)

#### `microgel/scripts/plot_efield_compare.py`
> Radial field around one subgrid charge, with and without the pm_sr
> short-range correction, on the same axes.
>
> Note the correction enters the stored efield only through
> pm_sr_apply_force_correction (the particle<->node term). The particle-particle
> pair correction never touches the field, and with a single particle it would
> do nothing anyway.

Argumentos (argparse): `--no-corr`; `--corr`; `-n / --ntotal`; `--pos`; `-q / --charge`; `--epsilon`; `--ld`; `--sigma`; `--rcut`; `--rmax`; `--images`; `--kernel`; `-o / --output`

#### `microgel/scripts/plot_efield_single.py`
> Electric field around a single subgrid charge.
>
> Reads Ludwig's ascii efield dump (3 components per site, z fastest, x slowest)
> and draws the field in the plane through the particle, plus the profile along
> the axis through it.
>
> The dump stores kT*E in lattice units (see psi_force.c), so the colour scale is
> that stored quantity, not E itself.

Argumentos (argparse): `efield`; `-n / --ntotal`; `--pos`: particle position x y z; `--ld`: Debye length, for the label; `--kernel`; `-o / --output`

#### `microgel/scripts/plot_efield_vs_r.py`
> Radial profile of the lattice electric field around a single subgrid charge,
> against theory.
>
> The dump stores kT*E (psi_force.c), which is also the force per unit charge:
> F = kt * reunit * q * Esub. So the stored field is compared directly against
> the force-per-unit-charge expressions, in the same units used for the
> two-particle plots.
>
> Three references:
>   - point Coulomb      q/(4 pi eps r^2)
>   - point Debye-Huckel q/(4 pi eps) exp(-kr) (k/r + 1/r^2)
>   - Gaussian cloud in a Debye-Huckel medium, width sigma (NOT sigma*sqrt(2):
> … (3 líneas más en el archivo)

Argumentos (argparse): `efield`; `-n / --ntotal`; `--pos`; `-q / --charge`; `--epsilon`; `--ld`: Debye length; `--sigma`: kernel sigma; `--rmax`; `--images`; `--kernel`; `-o / --output`

#### `microgel/scripts/plot_electric_field.py`
> Script para graficar el campo eléctrico desde archivos psi de Ludwig.
>
> El campo eléctrico se calcula como E = -∇ψ, donde ψ es el potencial eléctrico.
>
> Modos de visualización:
> 1. Plano 2D: Muestra el campo eléctrico en un plano específico (XY, XZ, o YZ)
> 2. Línea 1D: Muestra el campo eléctrico a lo largo de una línea en el espacio
> 3. Plano 3D: Muestra superficie 3D del potencial o campo eléctrico en un plano
>
> SISTEMA DE COORDENADAS:
> - Los nodos de la malla están centrados en: (0.5, 1.5, 2.5, ..., nx-0.5)
> - Los coloides en colloids-*.csv tienen coordenadas continuas (ej: 16.0, 16.0, 16.0)
> … (4 líneas más en el archivo)

Argumentos (argparse): `-f / --file`: Archivo psi a leer; `-s / --size`: Tamaño de la malla (nx ny nz); `-m / --mode`: Modo de graficación: plane (plano 2D), line (línea 1D), o plane3d (superficie 3D); `-p / --plane`: Plano a graficar (obligatorio para mode=plane); `-pos / --position`: Posición del plano (índice). Por defecto: centro; `-v / --vectors`: Mostrar vectores del campo (solo para mode=plane); `--vector-stride`: Espaciado entre vectores (por defecto: 2); `--start`: Punto inicial de la línea (obligatorio para mode=line); `--end`: Punto final de la línea (obligatorio para mode=line); `--num-points`: Número de puntos a lo largo de la línea (por defecto: 100); `-c / --component`: Componente a graficar. Plane/Plane3D: magnitude/x/y/z/psi (default: magnitude). Line: all/magnitude/x/y/z/psi (default: magnitude). Usar "all" en modo line muestra todas las componentes juntas; `--colormap`: Mapa de colores para modo plane3d (default: viridis). Opciones: viridis, plasma, inferno, magma, RdBu_r, coolwarm, seismic; `--elevation`: Ángulo de elevación de la vista 3D en grados (default: 30); `--azimuth`: Ángulo azimutal de la vista 3D en grados (default: -60); `--external-field`: Campo eléctrico externo constante a restar (Ex Ey Ez). Útil para visualizar solo el campo generado por cargas, sin el campo aplicado externamente.; `-o / --output`: Archivo de salida (por defecto: mostrar en pantalla)

#### `microgel/scripts/plot_field_comparizon.sh`
> Run multiple plots to compare teoretical Electric field to simulations

#### `microgel/scripts/plot_mui_plane.py`
> Grafica planos de potencial químico μ_i centrados en la partícula.
>
> El plano se define por su vector normal. El corte pasa por el nodo de la malla
> más cercano a la posición de la partícula proyectado sobre ese eje.
>
> Los archivos mui-NNNNNNNNN.001-001 tienen el mismo formato que qsi:
>   N filas, 2 columnas (μ₊, μ₋), orden x-major de Ludwig.
>
> La posición de la partícula se lee de config.cdsNNNNNNNN.001-001, línea 36.
>
> Uso:
>     python plot_mui_plane.py -n 530 -L 32 [-d directorio]
> … (6 líneas más en el archivo)

Uso (del archivo):
    python plot_mui_plane.py -n 530 -L 32 --normal 0,1,0   # plano XZ
    python plot_mui_plane.py -n 530 -L 32 --normal 1,0,0 --normal 0,1,0 --normal 0,0,1
    python plot_mui_plane.py -n 530 -L 32 --all-planes
    python plot_mui_plane.py --n-start 100 --n-end 500 --n-step 100 -L 32 --all-planes

Argumentos (argparse): `-n / --nstep`: Paso de tiempo único; `--n-start`; `--n-end`; `--n-step`; `-L / --grid-size`: Tamaño de la malla cúbica L×L×L; `-d / --directory`: Directorio con los archivos (default: .); `--normal`: Vector normal al plano. Se puede repetir. Use 'x','y','z' o 'nx,ny,nz'. (default: z); `--all-planes`: Generar los tres planos ortogonales (x, y, z); `--output-prefix`: Prefijo del nombre de los archivos de salida (default: mui_plane)

#### `microgel/scripts/plot_pin_to_grid.py`
> Lattice pinning: spurious force on an isolated particle as it is moved
> across one lattice cell. By symmetry the force on a single charge in an
> otherwise uniform periodic box must vanish; whatever is left measures how
> strongly the spread/gather kernel couples to the mesh.
>
> Reads PIN2GRID-hann{04,06,08}/x*/proceced_data/particle_force.csv.

#### `microgel/scripts/plot_pin_to_grid_all.py`
> Lattice pinning (pin-to-grid) for all spread/gather kernels, lambda_D = 8, L = 48.
>
> Left: F_x on the isolated particle, linear axis. Right: |F_x|, log axis (zeros
> by symmetry at the nodes and at the cell midpoint are round-off and are left out).

#### `microgel/scripts/plot_pin_to_grid_gauss.py`
> Lattice pinning for every spread/gather kernel available.
>
> Left: linear axis, Gaussian vs Hann-8 (same support, same width, same cost
> class). Right: |F_x| on a log axis for all four kernels. The points on nodes
> and at the cell midpoint are zero by symmetry (~1e-22, round-off) and are left
> out of the log panel.

#### `microgel/scripts/plot_psi_ewald.py`
> Script para graficar el potencial electrostático ψ calculado por Ludwig (Ewald),
> comparándolo con la solución teórica de Debye-Hückel.
> Soporta las componentes psi_real y psi_fourier de la misma manera que
> plot_componentes_ewald.py maneja efield_real y efield_fourier.
>
> Datos de entrada: archivos psi-NNNNNNNN.001-001 generados por Ludwig.
> Opcionalmente: psi_real-NNNNNNNN.001-001 y psi_fourier-NNNNNNNN.001-001.
>
> Uso básico:
>     python plot_psi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1.0e4 \
>         --rho-el 8.0e-3 --alpha 0.7 --rc 3.0
>
> … (20 líneas más en el archivo)

Argumentos (argparse): `-n / --nstep`: Numero de paso de tiempo; `-L / --grid-size`: Tamano de la malla cubica L x L x L; `-q / --charge`: Carga de la particula; `-kt / --kt`: Energia termica kT; `-eps / --epsilon`: Permitividad relativa; `--alpha`: Parametro de splitting de Ewald alpha; `--rc`: Radio de corte del espacio real; `--kappa`: Parametro de Debye kappa (si no se da, se calcula desde --rho-el); `--rho-el`: Densidad de electrolito para calcular kappa; `-d / --directory`: Directorio con los archivos de Ludwig (default: .); `--rmin`: Radio minimo a graficar (default: 0.5); `--rmax`: Radio maximo (default: L/2); `--show-components`: Mostrar psi_real y psi_fourier como series separadas; `--no-total`: Ocultar psi total de Ludwig; `--show-ewald-theory`: Mostrar curva teorica Ewald total (erfc+erf)/r; `--show-dh`: Mostrar curva DH infinita; `--show-dh-periodic`: Mostrar curva DH periodica; `--dh-periodic-shells`: Capas de imagenes periodicas (default: 1); `--pointwise-error`: Calcular error en la posicion exacta de cada nodo (verde); `--no-error-inf`: Ocultar serie de error vs DH infinito (purpura); `--no-error-dir`: Ocultar serie de error vs DH periodico direccional (naranja); `--particle-radius`: Radio de la partícula a (lu). Activa curva DH radio finito y error vs DH(a) en el panel de error (default: 0 = desactivado); `--only-error-dh-per-rad`: Mostrar en el panel de error SOLO la curva de error relativo vs DH periódico con partícula de radio a. Requiere --particle-radius > 0 y --show-dh-periodic. Suprime las demás curvas de error.; `--gaussian-source`: Usar distribución de carga gaussiana ρ(r)=Q/(π^(3/2)σ³)exp(-r²/σ²) como fuente en la teoría DH (en lugar de carga puntual). Activa la curva teórica ψ_gauss-DH en el panel 1.; `--sigma`: Anchura σ de la distribución gaussiana de carga (lu). Solo se usa con --gaussian-source (default: 1.0); `--gaussian-error`: Mostrar el error relativo vs la teoría gaussiana-DH (sistema infinito) en el panel 2. Requiere --gaussian-source.; `--gaussian-error-periodic`: Mostrar el error relativo vs la teoría gaussiana-DH periódica en el panel 2. Requiere --gaussian-source y --show-dh-periodic. Usa --dh-periodic-shells para el número de capas.; `--no-gaussian-dh-inf`: No graficar la curva gaussiana-DH sistema infinito (no periódica).; `--linear`: Escala lineal en Y (default: semilog); `--log-log`: Escala log-log en ambos ejes; `--log-log-error`: Escala log en eje Y del panel de error; `--xmin`; `--ymin`; `--subtract-background`: Restar el potencial de fondo de contraiones (promedio en nodos lejanos); `--bg-rmin`: Radio minimo para calcular el fondo (default: 12.0); `--rise-to-DHp`: Subir todos los valores de psi para que el minimo de los puntos seleccionados coincida con la curva DH-periodica en ese punto. Requiere --show-dh-periodic y --dh-periodic-shells.; `--direction`: Filtrar por direccion: 'x','y','z','xy','xyz' o 'dx,dy,dz'; `--angle-tol`: Tolerancia angular en grados (default: 15); `--output`: Nombre base del archivo de salida (default: psi_ewald); `--max-points`: Maximo numero de puntos a graficar (default: 10000); `--csv`: Exportar datos a CSV (nombre base; se anade -n_NNNNNN.csv)

#### `microgel/scripts/plot_qsi_ewald.py`
> Script para graficar la densidad de carga (qsi) calculada por Ludwig (Ewald),
> comparándola con la solución teórica de Debye-Hückel.
>
> Los archivos qsi contienen la densidad de carga por especie en cada nodo:
>   - Columna 0: densidad de cationes ρ₊ (≥ 0)
>   - Columna 1: densidad de aniones  ρ₋ (≥ 0)
>   Carga neta: ρ_net = ρ₊ - ρ₋
>
> La densidad de carga neta teórica de Debye-Hückel es:
>   ρ_net(r) = -κ² ε kT / e · ψ_DH(r) / (4π)
> equivalente a:
>   ρ_net(r) = -κ² q exp(-κ r) / (4π r)    [en unidades de Ludwig]
> … (13 líneas más en el archivo)

Uso (del archivo):
    Uso básico:
    python plot_qsi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1.0e4 \
    python plot_qsi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1.0e4 \
    python plot_qsi_ewald.py -n 2000 -L 32 -q 1.0 -kt 0.00001 -eps 1.0e4 \

Argumentos (argparse): `-n / --nstep`: Número de paso de tiempo; `-L / --grid-size`: Tamaño de malla cúbica L×L×L; `-q / --charge`: Carga de la partícula; `-kt / --kt`: Energía térmica kT; `-eps / --epsilon`: Permitividad relativa; `--kappa`: Parámetro de Debye kappa (si no se da, se calcula desde --rho-el); `--rho-el`: Densidad de electrolito para calcular kappa; `-d / --directory`: Directorio con los archivos de Ludwig (default: .); `--rmin`: Radio mínimo a graficar (default: 0.5); `--rmax`: Radio máximo (default: L/2); `--particle-radius`: Radio de la partícula (lu). Si > 0, agrega curva DH con radio finito: ρ(r)=-κ²q exp(-κ(r-a))/(4πr(1+κa)) (default: 0 = DH puntual); `--num-species`: Número de especies iónicas en el archivo qsi (default: 2); `--show-species0`: Mostrar densidad de especie 0 (cationes); `--show-species1`: Mostrar densidad de especie 1 (aniones); `--show-net`: Mostrar carga neta ρ_net = ρ₊ - ρ₋; `--no-net`: Ocultar carga neta; `--show-dh`: Mostrar curva ρ_net DH infinita; `--show-dh-periodic`: Mostrar curva ρ_net DH periódica (evaluada a lo largo de la dirección seleccionada); `--show-dh-periodic-pts`: Mostrar ρ_net DH periódico evaluado en la posición 3D exacta de cada nodo (más preciso; usa --dh-periodic-shells); `--dh-periodic-shells`: Capas de imágenes periódicas (default: 1); `--pointwise-error`: Calcular error en la posición exacta de cada nodo; `--no-error-inf`: Ocultar error vs DH infinito; `--no-error-dir`: Ocultar error vs DH periódico direccional; `--only-error-dh-per-rad`: Mostrar en el panel de error SOLO la curva de error relativo vs DH periódico con partícula de radio a. Requiere --particle-radius > 0 y --show-dh-periodic. Suprime las demás curvas de error.; `--gaussian-source`: Usar distribución gaussiana ρ(r)=Q/(π^(3/2)σ³)exp(-r²/σ²) como fuente en la teoría DH. Activa curva ρ_ions gaussiana-DH.; `--sigma`: Anchura σ de la distribución gaussiana de carga (lu). Solo se usa con --gaussian-source (default: 1.0); `--gaussian-error`: Mostrar error relativo vs teoría gaussiana-DH infinita en el panel 2. Requiere --gaussian-source.; `--gaussian-error-periodic`: Mostrar error relativo vs teoría gaussiana-DH periódica en el panel 2. Requiere --gaussian-source y --show-dh-periodic.; `--direction`: Filtrar por dirección: 'x','y','z','xy','xyz' o 'dx,dy,dz'; `--angle-tol`: Tolerancia angular en grados (default: 15); `--semilog`: Escala semilog en Y (default: lineal; usa |valor| para log); `--log-log`: Escala log-log en ambos ejes (usa |valor| automáticamente); `--log-log-error`: Escala log en eje Y del panel de error; `--show-error`: Mostrar panel de error relativo punto a punto vs DH inf (útil cuando se grafican todos los nodos sin --direction); `--error-ymax`: Límite superior del eje Y del panel de error en %% (útil para recortar outliers, ej: --error-ymax 200); `--xmin`; `--ymin`; `--output`: Nombre base del archivo de salida (default: qsi_ewald); `--max-points`: Máximo número de puntos a graficar (default: 10000); `--csv`: Exportar datos a CSV (nombre base; se añade -n_NNNNNN.csv); `--radial`: Graficar promedio radial (media ± std) en cascaras esféricas; `--radial-scatter`: Graficar todos los puntos individuales coloreados por cascara (en lugar de la media). Requiere --radial o se activa solo.; `--shell-width`: Grosor de cada cascara esférica en lu (default: 1.0); `--csv-radial`: Exportar promedio radial a CSV (nombre base); `--csv-scatter`: Exportar nodos individuales scatter a CSV (nombre base); `--angular-map`: Generar mapa de calor (θ,φ) en una cáscara esférica; `--angular-r`: Radio central de la cáscara para el mapa angular (lu); `--angular-dr`: Grosor de la cáscara para el mapa angular (default: 1.0 lu); `--angular-n-theta`: Resolución en θ del mapa angular (default: 90); `--angular-n-phi`: Resolución en φ del mapa angular (default: 180); `--angular-abs`: Usar |ρ_net| en lugar de ρ_net en el mapa angular; `--plane-map`: Generar mapa 2D de densidad de carga sobre un plano; `--plane`: Normal al plano (el plano pasa por la partícula): 'xy','xz','yz' o 'nx,ny,nz' (default: xy); `--plane-dz`: Semiancho de la rebanada en la dirección normal (default: 0.5 lu); `--plane-type`: Tipo de gráfico: 'pcolor' (sin interpolar) o 'contour' (default: pcolor); `--plane-contour-levels`: Número de niveles para contourf (default: 20); `--plane-abs`: Usar |ρ_net| en lugar de ρ_net en el mapa de plano

#### `microgel/scripts/plot_qsi_subgrid.py`
> plot_qsi_subgrid.py
>
> Grafica la distribución de carga del fluido (nodos qsi) junto con los
> subpuntos generados por subdivisión subgrid, comparando ambos con la
> curva teórica de Debye-Hückel en función de la distancia a la partícula.
>
> Muestra:
>   - Puntos de nodo (ρ_net por voxel) en gris
>   - Subpuntos (q_sub / volumen_sub) en color por voxel o por shell
>   - Curva teórica DH: ρ_net(r) = -κ²·q/(4π·r)·exp(-κr)
>
> Uso:
> … (6 líneas más en el archivo)

Uso (del archivo):
    python plot_qsi_subgrid.py -n 1000 -L 32 -q 1.0 -eps 1e4 \
    python plot_qsi_subgrid.py -n 1000 -L 32 -q 1.0 -eps 1e4 \

Argumentos (argparse): `-n / --nstep`; `-L / --grid-size`; `-q / --charge`; `-eps / --epsilon`; `--kappa`: Parámetro de Debye kappa (si no se da, se calcula desde --rho-el); `--rho-el`: Densidad de electrolito para calcular kappa; `-kt / --kt`; `--num-species`; `-d / --directory`; `--qsi-file`; `--cds-file`; `--rcut`: Radio de corte para voxels; `--rmin`: Radio mínimo para nodos (default: 0.5); `--rmax`: Radio máximo para el plot (default: rcut); `--rmin-sub`: Radio mínimo de exclusión para subpuntos (default: 0.0); `--rcut-sub`: Aplicar rcut también a subpuntos individuales; `--nsub`; `--gradient-dist`; `--quadratic-dist`; `--quartic-dist`; `--quadratic-full`; `--quartic-full`; `--no-nodes`: No mostrar los nodos qsi; `--no-subs`: No mostrar los subpuntos; `--no-dh`: No mostrar la curva DH; `--log-y`: Escala logarítmica en Y (usa |ρ|); `--output`: Nombre base del archivo de salida (default: qsi_subgrid)

#### `microgel/scripts/plot_velocities_grouped.py`
> Script para combinar plots de velocidad agrupados por un parámetro específico
> Lee los archivos datosfluid.csv y genera plots combinados
>
> Uso:
>     ./plot_velocities_grouped.py -d <directorio_padre> --group-by <parámetro> [opciones]
>
> Parámetros disponibles para agrupar:
>     - stencil: Agrupa por tipo de stencil (7, 19, 27, etc.)
>     - rhoel: Agrupa por valor de rhoel
>     - pos: Agrupa por posición
>     - q: Agrupa por carga
>     - eps: Agrupa por epsilon
> … (5 líneas más en el archivo)

Argumentos (argparse): `-d / --dir`: Directorio padre que contiene los datos; `--group-by`: Parámetro(s) para agrupar. Puede ser uno o varios separados por espacio (ej: stencil, stencil rhoel, stencil rhoel pos); `--plot-type`: Tipo de plot: both (partícula y fluido), particle (solo partícula), fluid (solo fluido); `--out_dir`: Directorio de salida para los plots (se usa si no se especifica -o); `-o / --output`: Ruta completa del archivo de salida (ejemplo: ./plots/mi_grafico.png). Si se especifica, se ignora --out_dir; `--component`: Componente de velocidad a plotear (all, x, y, z)

#### `microgel/scripts/plot_velocity_field.py`
> Script para visualizar campos de velocidad de archivos vel-* de Ludwig
> Permite graficar en planos xy, xz, yz o cualquier plano definido por el usuario

Argumentos (argparse): `-f / --file`: Archivo de velocidad específico; `-d / --directory`: Directorio con archivos de velocidad; `-t / --timestep`: Número de paso temporal (requerido si se usa -d); `-s / --size`: Tamaño de la malla (nx ny nz). Si no se especifica, asume malla cúbica; `--plane`: Plano a graficar (default: all); `--index`: Índice del plano (default: centro de la malla); `--skip`: Factor de submuestreo para flechas (default: 1); `--scale`: Escala para las flechas (default: auto); `--colormap`: Mapa de colores (default: viridis); `--magnitude-only`: Graficar solo la magnitud en plano sin vectores; `--single-plot`: Graficar solo un panel con magnitud y vectores (en lugar de 2 paneles); `-o / --output`: Directorio de salida (default: .); `--prefix`: Prefijo para archivos de salida (default: vel); `--show`: Mostrar gráficos interactivamente

#### `microgel/scripts/plot_velocity_field.sh`
> Run multiple plots to compare teoretical Electric field to simulations

#### `microgel/scripts/plot_verificar_simetria_radial.sh`
> Run multiple plots to compare teoretical Electric field to simulations

#### `microgel/scripts/plotdatos.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-i`: input file name; `-t / --tipo`: type of data to plot: 'posicion' or 'velocidad' (default: posicion); `-o / --output`: output file name (default: auto-generated based on type)

#### `microgel/scripts/plotvel.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-i`: input file name; `--out_dir`: Output directory; `--fluid-only`: plot only fluid average velocity

#### `microgel/scripts/qsi_subgrid_force.py`
> qsi_subgrid_force.py
>
> Calcula la fuerza electrostática que ejerce la distribución de carga del
> fluido sobre la partícula coloidal, usando un refinamiento subgrid:
>
>   - Cada voxel de la malla de Ludwig contiene una densidad de carga ρ_net.
>   - La carga total del voxel (ρ_net × 1 lu³) se redistribuye uniformemente
>     en n×n×n subpartículas equiespaciadas dentro del voxel.
>   - Las subpartículas se ubican en posiciones interiores al voxel: la
>     distancia desde el centro del voxel al borde es 0.5 lu, y los n puntos
>     se colocan de modo que la distancia de cada subpartícula a la pared más
>     cercana sea 1/(2*(n+1)) lu (división del arista en n+1 segmentos iguales).
> … (13 líneas más en el archivo)

Uso (del archivo):
    Uso:
    python qsi_subgrid_force.py -n 1000 -L 32 -q 1.0 -eps 1e4 \

Argumentos (argparse): `-n / --nstep`: Número de paso de tiempo; `-L / --grid-size`: Tamaño de la malla cúbica (NxNxN); `-q / --charge`: Carga de la partícula coloidal (default: 1.0); `-eps / --epsilon`: Permitividad dieléctrica (default: 1.0); `--rcut`: Radio de corte: solo voxels dentro de este radio; `--nsub`: Subdivisión por arista del voxel (default: 3 → 27 subpartículas); `--rmin-sub`: Radio mínimo desde la partícula: subpartículas más cercanas se excluyen y su carga se redistribuye entre las restantes (default: 0.0, sin exclusión); `--rcut-sub`: Aplicar rcut también a las subpartículas individuales. Por defecto solo se filtra el voxel por rcut; con esta opción las subpartículas que superen rcut también se excluyen y su carga se redistribuye entre las activas del mismo voxel.; `--gradient-dist`: Distribuir carga según gradiente lineal de ρ_net en el voxel (diferencias centradas). Por defecto: distribución uniforme.; `--quadratic-dist`: Distribuir carga según interpolación cuadrática de Lagrange usando los 3 nodos vecinos en cada eje (i-1, i, i+1). Separable en x, y, z; normalizado para conservar q_vox.; `--quartic-dist`: Distribuir carga según interpolación de Lagrange de grado 4 usando 5 nodos en cada eje (i-2, i-1, i, i+1, i+2). Separable en x, y, z; normalizado para conservar q_vox.; `--quadratic-full`: Lagrange grado 2 con producto tensorial completo (3x3x3=27 nodos). Captura correlaciones cruzadas diagonales. Más costoso que --quadratic-dist.; `--quartic-full`: Lagrange grado 4 con producto tensorial completo (5x5x5=125 nodos). Captura correlaciones cruzadas diagonales. Más costoso que --quartic-dist.; `--num-species`: Número de especies en el archivo qsi (default: 2); `--csv`: Nombre base del archivo CSV de salida; `--no-csv`: No exportar CSV aunque se especifique --csv; `--qsi-file`: Ruta explícita al archivo qsi (si no se usa el nombre estándar); `--cds-file`: Ruta explícita al archivo config.cds (si no se usa el nombre estándar)

#### `microgel/scripts/run_force_battery.sh`
> Corre baterias de fuerza entre 2 particulas (dirs con pos_*), N en paralelo. Reanudable:
> saltea las posiciones terminadas y rehace desde cero las que quedaron a medias.
> Uso: run_force_battery.sh NPAR dir1 [dir2 ...]   (las posiciones se intercalan entre dirs)

#### `microgel/scripts/run_pin_battery.sh`
> Corre baterias de pin-to-grid (dirs con x*), 4 en paralelo. Reanudable.
> Uso: run_pin_battery.sh dir1 [dir2 ...]

#### `microgel/scripts/run_template.sh`
> Plantilla de lanzamiento para corridas de microgel: corre Ludwig.exe y, en
> paralelo, va convirtiendo cada config.cds recien volcado a colloids-*.csv
> en colloid_data/ (via csv_watch.sh / extract_colloids), para poder abrir
> la corrida en Paraview mientras avanza -- igual que hace runbg.sh con
> runplot.sh/coloideacsv.sh para las corridas de una sola particula.

#### `microgel/scripts/runbg.sh`
> Run ludwig microgel simulation with MPI support
> Run './runbg.sh --help' for detailed usage information

#### `microgel/scripts/runmulti.sh`
> Run multiple Ludwig simulations with varying colloid parameters
> This script modifies config.cds.init.001-001 file to change charge, position,
> and al parameter for each simulation, then calls runbg.sh

#### `microgel/scripts/runplot.sh`
> Run ludwig microgel simulation plots
> Los .py de aca abajo se invocan via el shebang "#!/usr/bin/python3" si se
> corren como "./scripts/foo.py", que ignora PATH/venv activado. Para que
> usen numpy/scipy/matplotlib del venv del repo (y no el python del sistema,
> que no los tiene) se invocan explicitamente con PY_VENV.

#### `microgel/scripts/set_charges.py`
> Escribe un config.cds a partir de la estructura base equilibrada, aplicando un
> perfil de cargas. La geometria (posiciones, enlaces) nunca se toca: todas las
> corridas comparten exactamente la misma estructura, de modo que cualquier
> diferencia que se mida viene del perfil de carga o de la sal, no del punto de
> partida.
>
> Indices dentro del bloque de cada coloide en el formato config.cds:
>   0  = indice (id) del coloide      34 = posicion
>   2  = numero de enlaces            35 = velocidad lineal
>   9+ = ids de los coloides ligados  36 = velocidad angular
>                                     45 = q0
>                                     46 = q1
> … (7 líneas más en el archivo)

Uso (del archivo):
    Uso:

Argumentos (argparse): `--mode`; `--q`: carga por monomero cargado; `--frac`: fraccion cargada en modo random (solo interior); `--seed`; `--rmin`: radio interno en modo shell: se carga todo r >= rmin; `--nq`: modo outer: cuantos monomeros externos se cargan. Fija la carga total (y por lo tanto la concentracion de contraiones) sin depender de cuantos monomeros tenga la estructura.; `--base`; `-o / --out`

#### `microgel/scripts/verificar_campo_electrico.py`
> Script para verificar el campo eléctrico calculado por Ludwig contra el teórico.
>
> Este script:
> 1. Lee archivos efield de Ludwig
> 2. Extrae la posición de la partícula del archivo colloids
> 3. Calcula el módulo del campo eléctrico |E| para TODOS los puntos
> 4. Ordena por distancia y compara con el campo teórico de Debye-Hückel
> 5. Genera gráficos scatter de todos los puntos vs teoría
>
> Modelo teórico (campo NO reducido, como en Ludwig):
>     - kappa = 0 (sin iones): E = (q/4πε) · 1/r²
>     - kappa > 0 (con iones): E = (q/4πε) · exp(-κr) · (κ/r + 1/r²)
> … (5 líneas más en el archivo)

Argumentos (argparse): `-e / --efield`: Archivo efield de Ludwig (ej: efield-00000001.001-001); `-c / --colloid`: Archivo de coloides (ej: colloids-00000001.csv); `-L / --grid-size`: Tamaño de la malla (asume cubo LxLxL); `-q / --charge`: Carga de la partícula (en unidades de Ludwig); `-eps / --epsilon`: Permitividad relativa ε_r; `-kt / --kt`: Energía térmica k_B*T (solo para referencia, NO se usa en fórmula de E); `--kappa`: Parámetro de Debye κ = 1/λ_D (default: 0.0 = Coulomb); `-o / --output`: Archivo de salida para el gráfico (default: efield_comparison.png); `--rmin`: Radio mínimo de análisis (default: 0.5); `--rmax`: Radio máximo de análisis (default: L/2); `--max-points`: Máximo número de puntos a graficar (default: 10000); `--linear`: Usar escala lineal en los gráficos (default: logarítmica)

#### `microgel/scripts/verificar_simetria_psi.py`
> Script para verificar la simetría del campo eléctrico con respecto a la posición de la partícula.
>
> Este script:
> 1. Lee archivos psi de Ludwig
> 2. Extrae la posición de la partícula del archivo colloids
> 3. Calcula el campo eléctrico E = -∇ψ
> 4. Verifica la simetría radial y en diferentes planos
> 5. Compara con modelo teórico de Poisson-Boltzmann (Debye-Hückel)
> 6. Genera gráficos comparativos y métricas cuantitativas
>
> Autor: Script para análisis de simulaciones Ludwig

Argumentos (argparse): `-f / --file`: Archivo psi a analizar; `-s / --size`: Tamaño de la malla (nx ny nz); `-c / --colloid-file`: Archivo config.cds*.001-001 o colloids*.csv para leer posición de la partícula; `-p / --position`: Posición de la partícula (x y z); `--stencil`: Tipo de stencil para calcular el gradiente: 7 (D3Q7), 19 (D3Q19), 27 (D3Q27). Default: 19; `--min-radius`: Radio mínimo para análisis radial (default: 0.01 para evitar divergencia en r=0); `--max-radius`: Radio máximo para análisis radial (default: 15.0); `--num-radii`: Número de radios a analizar (default: 10); `--shell-thickness`: Grosor de las cáscaras esféricas para muestreo radial (default: 0.1, solo si --use-shells); `--use-all-points`: Usar TODOS los puntos de la grilla agrupados por bins de distancia (default: True, recomendado); `--use-shells`: Usar método de cáscaras esféricas en lugar de todos los puntos (menos preciso); `--distance-bin-size`: Tamaño del bin para agrupar distancias cuando se usan todos los puntos (default: 0.05); `--planes`: Planos a analizar (default: xy xz yz); `--log-scale`: Escala logarítmica: "y" (solo eje Y), "xy" (ambos ejes log-log), "symlog" (symlog en Y), "symlog-xy" (log-X y symlog-Y), "none" (lineal); `--symlog-linthresh`: Umbral lineal para escala symlog (región lineal: [-linthresh, +linthresh]). Default: 1e-6; `--radial-only`: Solo realizar análisis radial (salta planos, componentes y métricas). Por defecto solo guarda el gráfico radial.; `--efield-only`: En modo radial, solo graficar campo eléctrico y su error (omite potencial). Requiere --radial-only.; `--save-csv`: Guardar perfiles radiales en archivos CSV (además de los gráficos). Por defecto solo se guardan gráficos.; `--psi-offset-factor`: Factor para calcular offset del potencial en escala log: offset = -psi_min * factor (default: 1.1); `--compare-theory`: Comparar con modelo teórico de Poisson-Boltzmann (Debye-Hückel); `--charge`: Carga de la partícula (default: 1.0); `--epsilon`: Permitividad relativa del medio (default: 1.0); `--kt`: Energía térmica k_B*T (default: 1.0); `--kappa`: Parámetro de Debye κ = 1/λ_D. κ=0 usa Coulomb (vacío), κ>0 usa Debye-Hückel (electrolito). (default: 0.0); `--no-fit-scale`: No ajustar el factor de escala teórico a la simulación. Usa directamente q/(4πε*kt) para comparar magnitudes absolutas.; `--theory-start-distance`: Distancia mínima desde donde graficar la curva teórica. Si no se especifica, usa la distancia mínima de los datos de simulación.; `-o / --output-prefix`: Prefijo para archivos de salida (default: mostrar en pantalla)


### `microgel/scripts/EFIELD-1PART-LD8-gauss08eq/`

#### `microgel/scripts/EFIELD-1PART-LD8-gauss08eq/run.sh`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**


### `microgel/scripts/INGO-1PART-LBTERM-DIAG/`

#### `microgel/scripts/INGO-1PART-LBTERM-DIAG/cloud_profile.py`
> Radial counterion profile around the subgrid particle (particle frame).
>
> For each dump step s (psi/rho every 500): read counterion density n (qsi column 1),
> fluid density rho, particle position (config.cds at step s, colloid_io_freq 100).
> Distances: periodic minimum image in a cubic box L=27, nodes at integer coords 1..L.
> Outputs (text): profile n(r)/nbar in bins of 1 LU, averaged over a time window;
> time series of the mean excess in r<3 and the excess counterion number within r<R;
> far-field (r>9) ion vs fluid statistics.

#### `microgel/scripts/INGO-1PART-LBTERM-DIAG/cloud_profile_sal.py`
> Radial profiles of BOTH ion species around the subgrid particle, with salt, vs the discrete PB reference
> (LONG-SALT-LD2LB-CEN-3E4, 2026-11-16). Two-species copy of cloud_profile.py.
>
> For each dump s > 1 (qsi/rho/config.cds every 1000):
>   * n_+(r)/nbar_+ and n_-(r)/nbar_-, nbar_k = N_k/L^3 (totals from step 1; NP conserves them), bins of 1 LU,
>     periodic minimum image around the particle position of that dump (particle frame);
>   * PB reference (pb_reference_sal.py) solved at THAT dump's particle position, binned the same way
>     (so the reference sees the same particle-vs-grid offsets as the data);
>   * time series: excess of n_- and n_+ in r < 3, net ion charge within r < 3 and r < 6
>     (Q_ion(R) = sum_{r<R} (n_+ - n_-), screening of the +1 particle), same for PB;
>   * far field r > 9: std of dn_k/nbar_k, std of drho/rho and their correlation.
> Window average of the profiles over dumps with step >= WIN0 (env, default 15000 = second half of the run;
> … (4 líneas más en el archivo)

Uso (del archivo):
    Usage: [WIN0=15000] python3 cloud_profile_sal.py <run dir>

Lee argumentos de línea de comandos: sys.argv[1]

#### `microgel/scripts/INGO-1PART-LBTERM-DIAG/correr_serie_sal.sh`
> SAL-SERIE (2026-11-16): corre los casos de crear_serie_sal.py, de a NPAR en paralelo (default 3),
> dejando >= 2100 MiB libres de GPU. NO se corrió al crearlo (pedido del usuario: solo armar la serie).
> Costo estimado: casos de dx = 1 ~18 min c/u (0.035 s/paso x 3e4); DXH-* ~4.7 h c/u (8x nodos, 6e4 pasos).
> Uso: NPAR=3 bash microgel/scripts/INGO-1PART-LBTERM-DIAG/correr_serie_sal.sh [caso ...]   (sin casos: todos)

#### `microgel/scripts/INGO-1PART-LBTERM-DIAG/crear_sal.py`
> LONG-SALT-LD2LB-CEN-3E4 (2026-11-16, pedido del usuario): LONG-MASS-3E5 (parámetros de Ingo, una carga q = +1
> libre, LB térmico) pero CON SAL 1:1 tal que lambda_D = 2 l_B (elección del usuario entre l_B, 2 l_B y 3 l_B),
> advección centred (elección del usuario) y 3e4 pasos (prueba corta, elección del usuario).
>
> Sal: lambda_D = (8 pi l_B n_s)^(-1/2) para sal 1:1 (fórmula que usa Ludwig para el "Debye length" del log)
>      -> n_s = 1/(8 pi l_B lambda_D^2) = 1/(8 pi 2.84 5.68^2) = 4.3426e-4 por especie.
> Los contraiones del +1 (1/27^3 = 5.08e-5) los agrega Ludwig por electroneutralidad: con ellos, lambda_D efectivo
> (4 pi l_B sum z^2 n)^(-1/2) = 5.52 (no 5.68); Ludwig informa solo el de la sal.
> Todo lo demás idéntico a LONG-MASS-3E5 (config.cds.init, semilla, kernel Hann-8, D de Ingo, salida cada 1000).
> Ejecutable: el actual de src/ (269ec1b8: pesos advectivos + no-flujo en GPU).
>
> Uso:  python3 microgel/scripts/INGO-1PART-LBTERM-DIAG/crear_sal.py
> … (1 líneas más en el archivo)

#### `microgel/scripts/INGO-1PART-LBTERM-DIAG/crear_serie_sal.py`
> SAL-SERIE (2026-11-16, pedido del usuario): cómo cambia la nube de iones alrededor de la carga q = +1 con sal
> (lambda_D = 2 l_B) al variar D_ion, la temperatura y la resolución, y si la nube sigue a la partícula.
> Base: LONG-SALT-LD2LB-CEN-3E4 (parámetros de Ingo, 27^3, kT 1e-4, eta 0.01, l_B 2.84, n_s 4.3426e-4,
> lambda_D 5.68, Hann-8, advección centred, LB térmico, partícula libre, 3e4 pasos). Esa corrida ES el caso
> BASE-LIBRE y no se repite aquí.
>
> Decisiones del usuario (2026-11-16):
>   * D_ion x10 y x30; la partícula conserva su movilidad desnuda mu_b = 1/gamma (a_h, a_L sin cambio).
>   * kT/4 y kT/16 con epsilon ~ 1/kT: l_B = 1/(4 pi eps kT) = 2.84 fijo, misma sal y lambda_D; solo bajan el ruido
>     del LB, el empujón de la partícula y D_p = kT mu_b. D_ion no cambia (Ludwig usa D y el potencial en kT).
>   * dx/2, 1.5e4 pasos físicos: mismo sistema físico con escalado difusivo del tiempo (dt/4):
>         caja 54^3; Hann-16 (mismo ancho físico); D_LU y nu_LU iguales (d, viscosity sin cambio, tau_LB 0.53);
> … (12 líneas más en el archivo)

Uso (del archivo):
    Uso:  python3 microgel/scripts/INGO-1PART-LBTERM-DIAG/crear_serie_sal.py      (crea; NO corre)
    bash microgel/scripts/INGO-1PART-LBTERM-DIAG/correr_serie_sal.sh [caso ...]   (para correr, cuando se decida)

#### `microgel/scripts/INGO-1PART-LBTERM-DIAG/particle_motion.py`
> Mean velocity / position of the free subgrid particle under thermal LB (2026-11-12).
> Reads config.cds* (line 36: r, line 37: v = displacement per step). Unwraps periodic jumps.

#### `microgel/scripts/INGO-1PART-LBTERM-DIAG/pb_reference.py`
> Discrete Poisson-Boltzmann equilibrium reference for INGO-1PART-LBTERM-DIAG
> (2026-11-12). Counterions only, canonical (fixed total), same discretisation as
> the Ludwig run, so it is the exact fixed point of Ludwig's NP + Poisson:
>
>   * Particle charge q = +1 spread with the separable Hann kernel (subgrid.c, p = 2):
>         w(r) = (1/n) [1 + cos(2 pi r/n)],  |r| <= n/2,  n = 8,
>     on nodes at integer coordinates 1..L (periodic), particle at r0.
>   * Poisson with the FFT solver's discrete 27-point Laplacian (psi_fft.c,
>     stencil_d3q27.c after the 20251203 fix): lambda(k) = 6 sum_p w_p (1 - cos k.c_p),
>     and phi = beta e psi (psi_fft.c multiplies rho by beta):
>         lambda(k) phi_hat(k) = 4 pi lB rho_hat(k),   k != 0;  phi_hat(0) = 0.
>   * NP zero-flux condition of Ludwig's link flux (Capuani form) is exactly
> … (5 líneas más en el archivo)

#### `microgel/scripts/INGO-1PART-LBTERM-DIAG/pb_reference_sal.py`
> Discrete Poisson-Boltzmann reference WITH SALT (two species) for LONG-SALT-LD2LB-CEN-3E4 (2026-11-16).
>
> Copy of pb_reference.py (counterions only) extended to a 1:1 electrolyte; same discretisation as the Ludwig run,
> so it is the fixed point of Ludwig's NP + Poisson for a particle held at r0:
>   * particle charge q = +1, separable Hann kernel of order 8 (subgrid.c), nodes at integer coords 1..L (periodic);
>   * Poisson with the FFT solver's 27-point Laplacian, lambda(k) = 6 sum_p w_p (1 - cos k.c_p), phi = beta e psi:
>         lambda(k) phi_hat(k) = 4 pi lB rho_hat(k),  k != 0,  rho = q_p + n_+ - n_-;
>   * zero Capuani flux <=> n_k,i proportional to exp(-z_k phi_i): n_+ ~ exp(-phi), n_- ~ exp(+phi);
>   * CANONICAL, one constraint per species (NP conserves each total): N_+ and N_- read from the run's step-1 qsi
>         n_+,i = N_+ exp(-phi_i)/sum_j exp(-phi_j),   n_-,i = N_- exp(phi_i)/sum_j exp(phi_j).
> Damped Picard iteration (n <- 0.5 n + 0.5 n_new), stop when max relative change < 1e-12 (same as pb_reference.py).
>
> … (3 líneas más en el archivo)

Uso (del archivo):
    Usage (functions are imported by cloud_profile_sal.py):
    python3 pb_reference_sal.py <run dir>      -> solves at the INITIAL particle position (step 1) and prints

Lee argumentos de línea de comandos: sys.argv[1]


### `microgel/scripts/MGEL-RED-SURF-NOSALT/`

#### `microgel/scripts/MGEL-RED-SURF-NOSALT/monitor_equilibrio.sh`
> Actualiza estado_equilibrio.txt cada 30 min hasta que Ludwig termine.

#### `microgel/scripts/MGEL-RED-SURF-NOSALT/run.sh`
> Plantilla de lanzamiento para corridas de microgel: corre Ludwig.exe y, en
> paralelo, va convirtiendo cada config.cds recien volcado a colloids-*.csv
> en colloid_data/ (via csv_watch.sh / extract_colloids), para poder abrir
> la corrida en Paraview mientras avanza -- igual que hace runbg.sh con
> runplot.sh/coloideacsv.sh para las corridas de una sola particula.


### `microgel/scripts/MGEL-RED-neutro/`

#### `microgel/scripts/MGEL-RED-neutro/run.sh`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**


### `microgel/scripts/MGEL-SURF-NOSALT-hann06/`

#### `microgel/scripts/MGEL-SURF-NOSALT-hann06/run.sh`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**


### `microgel/scripts/MICROGEL-RED/`

#### `microgel/scripts/MICROGEL-RED/gen.sh`
> Red esferica entrecruzada (microgel): colocacion por rechazo con separacion
> minima 2.3 y enlace de TODO par a distancia [2.5, 3.5], que es lo que da la
> conectividad alta (z=4.72) que el generador de cadena+entrecruzamientos no
> alcanza (z=1.92). Radio 8 y densidad de saturacion para igualar el tamano
> del microgel anterior (Rg 6.67 vs 6.41; diametro 16.0 vs 15.9).


### `microgel/scripts/NBODY-additivity/`

#### `microgel/scripts/NBODY-additivity/analyze.py`
> Aditividad: F_i(N cuerpos) vs suma de F_i medidas en los pares aislados.
> Fuerzas del paso 1 (nube ionica aun uniforme => electrostatica pura).

#### `microgel/scripts/NBODY-additivity/mkrun.py`
> Genera un directorio de corrida con N particulas en posiciones dadas.

Lee argumentos de línea de comandos: sys.argv[1, 2, 3]


### `microgel/scripts/PIN2GRID-gauss08/`

#### `microgel/scripts/PIN2GRID-gauss08/run.sh`
> Pin-to-grid con kernel gaussiano (support 8, sigma 1.44612), 4 en paralelo. Reanudable:
> saltea las posiciones terminadas y rehace desde cero las que quedaron a medias.


### `microgel/scripts/SLIDES-LD8-hann04-corr_no-dx2/`

#### `microgel/scripts/SLIDES-LD8-hann04-corr_no-dx2/correr.sh`
> Corre todas las posiciones r_* (2026-11-09), N en paralelo. Saltea las terminadas (no borra nada).


### `microgel/scripts/SLIDES-LD8-hann04-corr_no-dx2-r4-shift/`

#### `microgel/scripts/SLIDES-LD8-hann04-corr_no-dx2-r4-shift/correr.sh`
> Corre todas las posiciones d_* (2026-11-09), N en paralelo. Saltea las terminadas (no borra nada).


### `microgel/scripts/SLIDES-LD8-hann08-corr_no-N10000/`

#### `microgel/scripts/SLIDES-LD8-hann08-corr_no-N10000/correr.sh`
> Corre todas las posiciones pos_* (2026-11-09), N en paralelo. Saltea las terminadas (no borra nada).


### `microgel/scripts/SLIDES-LD8-hann08-corr_no-N10000-d4a5/`

#### `microgel/scripts/SLIDES-LD8-hann08-corr_no-N10000-d4a5/correr.sh`
> Corre todas las posiciones pos_* (2026-11-09), N en paralelo. Saltea las terminadas (no borra nada).


### `microgel/scripts/SLIDES-LD8-hann16-corr_no-dxhalf/`

#### `microgel/scripts/SLIDES-LD8-hann16-corr_no-dxhalf/correr.sh`
> Corre todas las posiciones r_* (2026-11-09), N en paralelo. Saltea las terminadas (no borra nada).


### `microgel/scripts/SLIDES-LD8-hann16-corr_no-dxhalf-cont/`

#### `microgel/scripts/SLIDES-LD8-hann16-corr_no-dxhalf-cont/correr_cont.sh`
> Resumes every r_* from its step N0 (N in parallel). Deletes nothing; log in run_cont.log.


### `microgel/scripts/SS-SCAN-repulsion/`

#### `microgel/scripts/SS-SCAN-repulsion/freepos.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Lee argumentos de línea de comandos: sys.argv[1, 2, 3]


### `microgel/scripts/TEST-HIDRO-dx1/`

#### `microgel/scripts/TEST-HIDRO-dx1/grafico_psi_E_hidro_no.py`
> Potential and electric field of the dx = 1, Hann-8, r = 4 run WITHOUT hydrodynamics (2026-11-11, user's request).
>
> Run: ../hidro_no_exe45dc94bc (box 48^3, charges q = +1 at x = 8 and 12, y = z = 24; exe 45dc94bc), last dump (step 30000).
> Units: psi dump in e psi/kT; efield dump / kT -> e E/kT per lattice unit (dx = 1: per physical length).
> Reference (user's choice): linear Debye-Hueckel, periodic box, continuous Hann-8 sources = dh_nodes of
> error_capas_2part.py (zero node mean, E_DH = -grad phi_DH), evaluated on the nodes, M = 3N modes per direction.
> Top: profiles on the line y = z = 24 (through both charges), Ludwig (points) and reference (line).
> Middle: pointwise error on that line, Ludwig - reference on the nodes (signed, absolute).
> Bottom: maps on the plane z = 24 (contains both charges): psi, and |E| with arrows of (E_x, E_y) every 3 nodes.
> Node at array index j <-> Ludwig coordinate j + 1.
> Usage: python3 grafico_psi_E_hidro_no.py   -> ../graficos/psi_E_hidro_no_paso30000.png


### `microgel/scripts/TEST-HIDRO-dx2/`

#### `microgel/scripts/TEST-HIDRO-dx2/analizar_visc.py`
> Viscosity test (2026-11-10, user's request): does the chemical-potential floor come from the steady flow?
>
> Runs (dx = 2, r = 4, current exe 4f0c20e7): ../hidro_si (viscosity 0.05, as SLIDES-LD8-hann04-corr_no-dx2/r_04) and
> ../hidro_si_eta0.5 (viscosity and bulk viscosity 0.5). If the floor of Delta mu comes from the flow, u ~ 1/eta should drop
> ~10x and Delta mu with it.
> Per dump (definitions of microgel/local/scripts/energia_2part.py): F/kT and dF between dumps, Delta mu_a = max - min of
> mu_a/kT = ln rho_a + Z_a psi, max|u| and RMS|u| (vel dump, lattice units).
> 2026-11-10: extended to the 7-viscosity sweep (0.02 ... 2) + summary table of the last dump.
> Last dump: RMS of delta mu = mu - <mu> per shell of physical distance to the nearest charge (width 1, as mu_capas_2part.py).
> Usage: python3 analizar_visc.py   (from this directory)

#### `microgel/scripts/TEST-HIDRO-dx2/comprobar_flla_descomposicion.py`
> Exact decomposition of the flla (mode 1) uniform-field source S, and its correlations (2026-11-16, user's request).
>
> Follows comprobar_flla_fuerza.py (same reconstruction of the post-collision f* from dist-<step>, reused from it).
>     S_i = D_i/rho_i + R_i,
>     D_i = sum_a [f*_abar(x+c_a) - f*_a(x)]                 net LB mass flux into x (= rho(t+1) - rho*)
>     R_i = sum_a f*_abar(x+c_a) (1/rho(x+c_a) - 1/rho(x))   part due to LB density differences between neighbours
> Prints per viscosity: RMS of S, D/rho, R, corr(S, R), RMS and max of rho - <rho>; and for eta = 0.02 and 2 the
> node-by-node correlation of S with -div(F/2), -div(rho u), -div(j*) (central differences; F/2 = j* - rho u).
> Outcome 2026-11-16: D ~ 0, S = R exactly -> the hypothesis "flla advects with u + F/2rho" was REJECTED; the source
> comes from the normalisation by the local LB density (documented limitation of flla in src/nernst_planck.c).
> Reads files only. Usage: python3 scripts/comprobar_flla_descomposicion.py

#### `microgel/scripts/TEST-HIDRO-dx2/comprobar_flla_fuerza.py`
> Check of the hypothesis "flla advects the ions with u + F/(2 rho), not with u" (2026-11-16, user's request).
>
> Background: with electrokinetics_advection flla the Delta mu floor does not decrease with viscosity
> (grafico_dmu_u_visc_3adv.py). flla takes the link mass flux from the POST-collision distributions f*
> (nernst_planck_lattice_flux_compute, src/nernst_planck.c), whose momentum is
>     j* = sum_a f*_a c_a = j_pre + F = rho u + F/2      (collision.c: mode j += F; u = (j_pre + F/2)/rho),
> while centred advection uses u. F/2 does not depend on the viscosity.
>
> Exact reconstruction, no steady-state assumption: dist-<step> is written after propagation, so
> f*_a(x) = f_dump_a(x + c_a) (propagation moves f*_a(x) to x + c_a). vel-<step> is u of the collision of the
> same step. Then, per node:
>     F/2 = j* - rho u                       (half the body force of that step)
> … (7 líneas más en el archivo)

Uso (del archivo):
    Reads files only. Usage: python3 scripts/comprobar_flla_fuerza.py [step]   (default 3000)

Lee argumentos de línea de comandos: sys.argv[1]

#### `microgel/scripts/TEST-HIDRO-dx2/correr_advfix.sh`
> Corre advfix/* (crear_advfix.py, 2026-11-16) en paralelo dejando >= 2100 MiB libres en la GPU
> (regla del usuario); cada corrida 24^3 usa ~904 MiB (medido 2026-11-10). No borra nada.

#### `microgel/scripts/TEST-HIDRO-dx2/correr_advfix_flla.sh`
> Corre advfix_flla/* (crear_advfix_flla.py, 2026-11-16) en paralelo dejando >= 2100 MiB libres en la GPU
> (regla del usuario); cada corrida 24^3 usa ~904 MiB (medido 2026-11-10). No borra nada.

#### `microgel/scripts/TEST-HIDRO-dx2/correr_advfix_flla_mass.sh`
> Corre advfix_flla_mass/* (crear_advfix_flla_mass.py, 2026-11-16) en paralelo dejando >= 2100 MiB libres en la GPU
> (regla del usuario); cada corrida 24^3 usa ~904 MiB (medido 2026-11-10). No borra nada.

#### `microgel/scripts/TEST-HIDRO-dx2/correr_cont_visc.sh`
> Corre las continuaciones cont<NEND>/<variante>/<corrida> (crear_cont_visc.py, 2026-11-16) en paralelo dejando
> >= 2100 MiB libres en la GPU (regla del usuario); ~904 MiB por corrida 24^3. No borra nada.
> [OK] exige "finished normally" y "Re-starting simulation at step N0".
> Uso: scripts/correr_cont_visc.sh NEND   (corre todas las carpetas de cont<NEND> que no hayan terminado)

#### `microgel/scripts/TEST-HIDRO-dx2/correr_exe45dc.sh`
> Prueba del ejecutable 45dc94bc (2026-11-11): hidro_si (debe reproducir hidro_si exacto) y hidro_no, en paralelo.

#### `microgel/scripts/TEST-HIDRO-dx2/correr_visc.sh`
> Barrido de viscosidad (2026-11-10). Lanza en PARALELO las viscosidades que entren en la GPU dejando siempre libre
> RESERVA MiB (pedido del usuario: "dejar una memoria libre como la que hay ahora"; medido ~2100-2300 MiB con las tres
> continuaciones de 96^3 corriendo). Cada corrida 24^3 usa ~904 MiB (medido): lanza una si libre >= RESERVA + 904,
> luego espera 60 s a que asigne su memoria antes de volver a medir. No borra nada; log de cada una en run.log.

#### `microgel/scripts/TEST-HIDRO-dx2/crear_advfix.py`
> Rehace el barrido de viscosidad de TEST-HIDRO-dx2 (7 corridas, 2026-11-10) con el ejecutable ACTUAL
> (src/Ludwig.exe, md5 221f3229, compilado 2026-11-12 19:35), que incluye el arreglo del peso de la
> advección de NP (CHANGE 20261112: W_p = 2 wgradients[p]; antes los iones se advectaban 9x con stencil 27).
> Pedido del usuario 2026-11-16: comparar el potencial químico con hidrodinámica en función de la
> viscosidad antes y después del arreglo.
>
> Cada corrida nueva es una copia del input y de config.cds.init de la original, sin otro cambio físico
> (misma advección "centred" por defecto, mismos 3000 pasos, mismo solver fft); solo se agregan líneas de
> comentario en la cabecera.
>
> Uso:  python3 scripts/crear_advfix.py   (desde TEST-HIDRO-dx2). Crea advfix/<corrida>/ (no sobrescribe).

#### `microgel/scripts/TEST-HIDRO-dx2/crear_advfix_flla.py`
> Barrido de viscosidad de TEST-HIDRO-dx2 (7 corridas, 2026-11-10) con el ejecutable actual (md5 221f3229)
> y advección de NP "flla" (DeVincentis & Thomas 2023, limitador superbee; CHANGE 20261112 en
> src/nernst_planck.c) en lugar de "centred". Pedido del usuario 2026-11-16: comparar Delta mu vs viscosidad
> con la advección vieja, la centrada corregida (advfix/) y flla.
>
> Cada corrida es una copia del input y de config.cds.init de la original; cambian el ejecutable y la
> línea agregada "electrokinetics_advection flla" (requisitos: nhalo >= 2 (aquí 3), un solo proceso MPI,
> enlaces D3Q19 del LB contenidos en el stencil 27 de NP).
>
> Uso:  python3 scripts/crear_advfix_flla.py   (desde TEST-HIDRO-dx2). Crea advfix_flla/<corrida>/ (no sobrescribe).

#### `microgel/scripts/TEST-HIDRO-dx2/crear_advfix_flla_mass.py`
> Barrido de viscosidad de TEST-HIDRO-dx2 (7 corridas, 2026-11-10) con el ejecutable actual (md5 221f3229)
> y advección de NP "flla_mass" (variante de flla de DeVincentis & Thomas 2023, limitador superbee; CHANGE 20261112 en
> src/nernst_planck.c) en lugar de "centred". Pedido del usuario 2026-11-16: comparar Delta mu vs viscosidad
> con la advección vieja, la centrada corregida (advfix/), flla (advfix_flla/) y flla_mass.
>
> Cada corrida es una copia del input y de config.cds.init de la original; cambian el ejecutable y la
> línea agregada "electrokinetics_advection flla_mass" (requisitos: nhalo >= 2 (aquí 3), un solo proceso MPI,
> enlaces D3Q19 del LB contenidos en el stencil 27 de NP).
>
> Uso:  python3 scripts/crear_advfix_flla_mass.py   (desde TEST-HIDRO-dx2). Crea advfix_flla_mass/<corrida>/ (no sobrescribe).

#### `microgel/scripts/TEST-HIDRO-dx2/crear_cont_visc.py`
> Continuation of the 4 x 7 viscosity runs of TEST-HIDRO-dx2 (2026-11-16, user's request: run 3000 more steps, check
> for steady state, and if not reached run 3000 more; then redo the Delta mu curve with the steady state).
>
> Variants (source of step N0):   old       ./hidro_si*              (exe 4f0c20e7, unweighted centred advection)
>                                 centred   ./advfix/hidro_si*       (exe 221f3229)
>                                 flla      ./advfix_flla/hidro_si*
>                                 flla_mass ./advfix_flla_mass/hidro_si*
> Destination: cont<NEND>/<variant>/<run>/. Recipe of crear_continuacion_dxhalf.py: copies of input, config.cds.init,
> the SAME Ludwig.exe as the source run, the *-metadata files, the files of step N0 only, and proceced_data/*.csv
> truncated to step <= N0. Input: N_start N0, N_cycles NEND - N0 (nothing else changed). Ludwig restarts from
> dist/psi/qsi/rho/vel of step N0 ("Re-starting simulation at step N0").
> Note: flla/flla_mass store the LB link flux after each collision; right after a restart there is none yet, so the
> … (6 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 scripts/crear_cont_visc.py N0 NEND [run ...]
    e.g. python3 scripts/crear_cont_visc.py 3000 6000            (all 28)
    python3 scripts/crear_cont_visc.py 6000 9000 flla/hidro_si_eta2.0 ...   (sources taken from cont6000/)

Lee argumentos de línea de comandos: sys.argv[1, 2]

#### `microgel/scripts/TEST-HIDRO-dx2/estacionario_visc.py`
> Steady-state check of the 4 x 7 viscosity runs (2026-11-16, user's request).
>
> For each variant (old, centred, flla, flla_mass) and viscosity, the time series of Delta mu_a = max - min of
> mu_a/kT = ln rho_a + Z_a psi (as analizar_visc.py) is built from the dumps every 250 steps of the original run
> (0..3000) and of its continuations cont6000/, cont9000/ (if present).
>
> Steady-state criterion (HEURISTIC, proposed by Claude, not a standard; stated to the user):
>     r_a = |Delta mu_a(N) - Delta mu_a(N - 1000)| / Delta mu_a(N) < TOL = 1 %   for both species,
> with N the last available dump. Printed: Delta mu_a at 3000, at N-1000 and at N, r_a, and STEADY yes/no.
> Writes nothing. Usage: python3 scripts/estacionario_visc.py

#### `microgel/scripts/TEST-HIDRO-dx2/grafico_dmu_u_visc.py`
> Two-panel plot vs viscosity: chemical-potential residual (top) and fluid speed (bottom) (2026-11-10, user's request).
>
> Runs: the 7-viscosity sweep of analizar_visc.py (dx = 2, r = 4, 24^3), last dump (step 3000).
> Top: Delta mu_a = max_n - min_n of mu_a/kT = ln rho_a + Z_a psi (as grafico_dmu_visc.py).
> Bottom: max_n |u_n| and RMS_n |u_n| (vel dump, lattice units), as analizar_visc.py.
> Shared log eta axis; the points are joined by straight segments only as a guide (no fit).
> Usage: python3 grafico_dmu_u_visc.py   (from this directory) -> ../graficos/dmu_u_vs_viscosidad.png
>
> 2026-11-16 (user's request): added to the top panel
>   * Delta mu with the weighted centred advection (advfix/, exe 221f3229, W_p = 2 wgradients[p]), same step 3000;
>     open markers + dashed = old exe (unweighted advection, the original curves), filled + solid = weighted;
>   * Delta mu WITHOUT hydrodynamics (hidro_no_exe45dc94bc, exe 45dc94bc, hydrodynamics no, eta irrelevant) at step
> … (11 líneas más en el archivo)

#### `microgel/scripts/TEST-HIDRO-dx2/grafico_dmu_u_visc_3adv.py`
> Delta mu and |u| vs viscosity for three NP advection variants (2026-11-16, user's request).
>
>   old     ../hidro_si*              exe 4f0c20e7, centred advection WITHOUT link weights (ions 9x too fast)
>   centred ../advfix/hidro_si*       exe 221f3229, centred advection with W_p = 2 wgradients[p] (crear_advfix.py)
>   flla    ../advfix_flla/hidro_si*  exe 221f3229, electrokinetics_advection flla (crear_advfix_flla.py)
> Same quantities and dump as grafico_dmu_u_visc.py: Delta mu_a = max_n - min_n of mu_a/kT = ln rho_a + Z_a psi
> (top), max_n|u_n| and RMS_n|u_n| (bottom), step 3000. Straight segments only guide the eye (no fit).
> Prints the table of values and of the ratios old/variant.
> Usage: python3 scripts/grafico_dmu_u_visc_3adv.py  -> graficos/dmu_u_vs_viscosidad_3adv.png

#### `microgel/scripts/TEST-HIDRO-dx2/grafico_dmu_u_visc_4adv.py`
> Delta mu and |u| vs viscosity for four NP advection variants (2026-11-16, user's request).
>
>   old     ../hidro_si*              exe 4f0c20e7, centred advection WITHOUT link weights (ions 9x too fast)
>   centred ../advfix/hidro_si*       exe 221f3229, centred advection with W_p = 2 wgradients[p] (crear_advfix.py)
>   flla    ../advfix_flla/hidro_si*  exe 221f3229, electrokinetics_advection flla (crear_advfix_flla.py)
>   flla_mass ../advfix_flla_mass/hidro_si*  exe 221f3229, electrokinetics_advection flla_mass
>           (crear_advfix_flla_mass.py)
> Same quantities and dump as grafico_dmu_u_visc.py: Delta mu_a = max_n - min_n of mu_a/kT = ln rho_a + Z_a psi
> (top), max_n|u_n| and RMS_n|u_n| (bottom), step 3000. Straight segments only guide the eye (no fit).
> Prints the table of values and of the ratios centred/variant.
> Usage: python3 scripts/grafico_dmu_u_visc_4adv.py  -> graficos/dmu_u_vs_viscosidad_4adv.png

#### `microgel/scripts/TEST-HIDRO-dx2/grafico_dmu_u_visc_4adv_estacionario.py`
> Delta mu and |u| vs viscosity for four NP advection variants, STEADY STATE (step 6000) (2026-11-16, user's request).
>
> Same as grafico_dmu_u_visc_4adv.py (step 3000) but read from the continuations cont6000/<variant>/<run>
> (crear_cont_visc.py; steady by estacionario_visc.py: relative change of Delta mu over steps 5000-6000
> <= 1.7e-6 in all 28 runs). Variants as listed below (the cont6000 copies use the same exe as the source).
>
>   old     ../hidro_si*              exe 4f0c20e7, centred advection WITHOUT link weights (ions 9x too fast)
>   centred ../advfix/hidro_si*       exe 221f3229, centred advection with W_p = 2 wgradients[p] (crear_advfix.py)
>   flla    ../advfix_flla/hidro_si*  exe 221f3229, electrokinetics_advection flla (crear_advfix_flla.py)
>   flla_mass ../advfix_flla_mass/hidro_si*  exe 221f3229, electrokinetics_advection flla_mass
>           (crear_advfix_flla_mass.py)
> Same quantities and dump as grafico_dmu_u_visc.py: Delta mu_a = max_n - min_n of mu_a/kT = ln rho_a + Z_a psi
> … (3 líneas más en el archivo)

Uso (del archivo):
    Usage: python3 scripts/grafico_dmu_u_visc_4adv_estacionario.py  -> graficos/dmu_u_vs_viscosidad_4adv_estacionario.png

#### `microgel/scripts/TEST-HIDRO-dx2/grafico_dmu_u_visc_advfix.py`
> Delta mu and |u| vs viscosity, BEFORE vs AFTER the NP advection-weight fix (2026-11-16, user's request).
>
> Before: the 7-viscosity sweep of 2026-11-10 (../hidro_si*, exe 4f0c20e7: unweighted centred advection, ions advected
> 9x too fast with stencil 27). After: the same 7 inputs rerun with the current exe (../advfix/hidro_si*, md5 221f3229,
> CHANGE 20261112 W_p = 2 wgradients[p]); see crear_advfix.py. Same quantities and dump as grafico_dmu_u_visc.py:
> Delta mu_a = max_n - min_n of mu_a/kT = ln rho_a + Z_a psi (top), max_n|u_n| and RMS_n|u_n| (bottom), step 3000.
> Straight segments only guide the eye (no fit). Prints a table with the before/after ratios.
> Usage: python3 scripts/grafico_dmu_u_visc_advfix.py  -> graficos/dmu_u_vs_viscosidad_advfix.png

#### `microgel/scripts/TEST-HIDRO-dx2/grafico_dmu_visc.py`
> Plot of the chemical-potential residual vs viscosity (2026-11-10, user's request).
>
> Runs: the 7-viscosity sweep of analizar_visc.py (dx = 2, r = 4, 24^3, 3000 steps), last dump (step 3000).
> Delta mu_a = max_n - min_n of mu_a/kT = ln rho_a + Z_a psi (same definition as energia_2part.py / analizar_visc.py).
> Log-log axes; the points are joined by straight segments only as a guide (no fit).
> Usage: python3 grafico_dmu_visc.py   (from this directory) -> ../graficos/dmu_vs_viscosidad.png


### `microgel/scripts/aaa-otro/`

#### `microgel/scripts/aaa-otro/analizar_campo_electrico.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

#### `microgel/scripts/aaa-otro/batch_plot_electric_field.py`
> Script para procesar múltiples archivos psi y generar gráficos del campo eléctrico.
>
> Este script busca todos los archivos psi en un directorio y genera gráficos
> para cada uno de ellos de forma automática.
>
> Autor: Script generado para análisis de simulaciones Ludwig

Argumentos (argparse): `-d / --directory`: Directorio donde buscar archivos psi (por defecto: directorio actual); `-s / --size`: Tamaño de la malla (nx ny nz); `-m / --mode`: Modo de graficación: plane (plano 2D), line (línea 1D), o plane3d (superficie 3D); `--pattern`: Patrón de búsqueda de archivos (por defecto: psi-*.001-001); `--output-dir`: Directorio de salida (por defecto: electric_field_plots); `-p / --plane`: Plano a graficar (obligatorio para mode=plane); `-pos / --position`: Posición del plano (índice). Por defecto: centro; `-v / --vectors`: Mostrar vectores del campo (solo para mode=plane); `--vector-stride`: Espaciado entre vectores (por defecto: 2); `--start`: Punto inicial de la línea (obligatorio para mode=line); `--end`: Punto final de la línea (obligatorio para mode=line); `--num-points`: Número de puntos a lo largo de la línea (por defecto: 100); `-c / --component`: Componente a graficar. Plane/Plane3D: magnitude/x/y/z/psi. Line: all/magnitude/x/y/z/psi. (default: magnitude); `--colormap`: Mapa de colores para modo plane3d (default: viridis); `--elevation`: Ángulo de elevación para modo plane3d (default: 30); `--azimuth`: Ángulo azimutal para modo plane3d (default: -60); `--max-files`: Procesar solo los primeros N archivos; `--skip`: Procesar cada N archivos (por defecto: 1, procesar todos)

#### `microgel/scripts/aaa-otro/calculos.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-nciclos`: number of cycles; `-npaso`: cycles per step; `-o`: output file name

#### `microgel/scripts/aaa-otro/calculosvel.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-nciclos`: number of cycles; `-npaso`: cycles per step; `-o`: output file name

#### `microgel/scripts/aaa-otro/calculosvelfluid.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-nciclos`: number of cycles; `-ninicio`: Number of cycle to start from; `-npaso`: cycles per step; `-o`: output file name; `--idir`: input directory

#### `microgel/scripts/aaa-otro/calculosvelfluidonly.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-nciclos`: number of cycles; `-ninicio`: Number of cycle to start from; `-npaso`: cycles per step; `-o`: output file name

#### `microgel/scripts/aaa-otro/coloideacsv.sh`
> Convert config colloid files by batch to .csv deleting vtk files
> Parameters

#### `microgel/scripts/aaa-otro/compare_field_theory.py`
> Script para comparar el campo eléctrico simulado con el campo teórico
> de una carga puntual.
>
> Para una carga puntual q en el vacío (o medio uniforme), el campo eléctrico es:
> E = (q / (4πε)) * (r / |r|³)
>
> donde r es el vector desde la carga al punto de observación.
>
> En Ludwig, el sistema está discretizado en una malla y se usa el sistema de
> unidades de lattice Boltzmann.

Argumentos (argparse): `-f / --file`: Archivo psi a leer; `-s / --size`: Tamaño de la malla; `--charge-pos`: Posición de la carga puntual; `--charge`: Magnitud de la carga (default: 1.0); `--epsilon`: Permitividad del medio (default: 1.0); `--kt`: Energía térmica k_B*T (default: 1.0); `--kappa`: Parámetro de Debye κ (inverso de longitud de Debye) para screening iónico. kappa=0 usa Coulomb simple (vacío). kappa>0 usa Debye-Hückel (electrolito). En unidades de lattice. (default: 0.0); `--ionic-strength`: Fuerza iónica en mol/L. Si se proporciona, calcula kappa automáticamente. Requiere conocer la escala de longitud del sistema.; `-m / --mode`: Modo de comparación: line o plane; `--start`: Punto inicial de la línea (obligatorio para mode=line); `--end`: Punto final de la línea (obligatorio para mode=line); `--num-points`: Número de puntos a lo largo de la línea (default: 100); `--show-nodes`: Mostrar solo los nodos de Ludwig con marcadores (sin interpolación); `-p / --plane`: Plano a comparar (obligatorio para mode=plane); `--position`: Posición del plano (índice). Default: centro; `-o / --output`: Archivo de salida (default: mostrar en pantalla); `--external-field`: Componentes del campo externo constante a restar (Ex, Ey, Ez). Ejemplo: --external-field 0.1 0.0 0.0

#### `microgel/scripts/aaa-otro/compare_field_theory_peskin.py`
> Script para comparar el campo eléctrico simulado (con interpolación Peskin)
> con el campo teórico de una carga puntual.
>
> DIFERENCIA con compare_field_theory.py:
> - Este script usa interpolación con delta de Peskin (igual que Ludwig)
> - compare_field_theory.py usa diferencias finitas estándar
>
> Para una carga puntual q en el vacío (o medio uniforme), el campo eléctrico es:
> E = (q / (4πε)) * (r / |r|³)
>
> donde r es el vector desde la carga al punto de observación.
>
> … (2 líneas más en el archivo)

Argumentos (argparse): `-f / --file`: Archivo psi a leer; `-s / --size`: Tamaño de la malla; `--charge-pos`: Posición de la carga puntual; `--charge`: Magnitud de la carga (default: 1.0); `--epsilon`: Permitividad del medio (default: 1.0); `-m / --mode`: Modo de comparación: line o plane; `--start`: Punto inicial de la línea (obligatorio para mode=line); `--end`: Punto final de la línea (obligatorio para mode=line); `--num-points`: Número de puntos a lo largo de la línea (default: 100); `-p / --plane`: Plano a comparar (obligatorio para mode=plane); `--position`: Posición del plano (índice). Default: centro; `-o / --output`: Archivo de salida (default: mostrar en pantalla)

#### `microgel/scripts/aaa-otro/compare_potential_theory.py`
> Script para comparar el potencial eléctrico simulado con el potencial teórico
> de una carga puntual.
>
> Para una carga puntual q, el potencial eléctrico es:
> ψ = q / (4πε|r|)
>
> donde |r| es la distancia desde la carga al punto de observación.

Argumentos (argparse): `-f / --file`: Archivo psi a leer; `-s / --size`: Tamaño de la malla; `--charge-pos`: Posición de la carga puntual; `--charge`: Magnitud de la carga (default: 1.0); `--epsilon`: Permitividad del medio (default: 1.0); `-m / --mode`: Modo de comparación: line o plane; `--start`: Punto inicial de la línea (obligatorio para mode=line); `--end`: Punto final de la línea (obligatorio para mode=line); `--num-points`: Número de puntos a lo largo de la línea (default: 100); `-p / --plane`: Plano a comparar (obligatorio para mode=plane); `--position`: Posición del plano (índice). Default: centro; `-o / --output`: Archivo de salida (default: mostrar en pantalla)

#### `microgel/scripts/aaa-otro/convert_to_pdf.py`
> Simple script to convert markdown to PDF using available tools

#### `microgel/scripts/aaa-otro/copy_plots.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-origen`: Carpeta origen; `-destino`: Carpeta destino

#### `microgel/scripts/aaa-otro/delete_files_dist_rho.sh`
> Delete unnecesary files.
> Leave last file so that restart is possible.

#### `microgel/scripts/aaa-otro/ejemplo_runmulti.sh`
> Ejemplo de uso de runmulti.sh
> Este script muestra diferentes ejemplos de cómo usar runmulti.sh

#### `microgel/scripts/aaa-otro/ejemplos_carga.sh`
> Ejemplos de uso del script de distribución de carga
> Instrucciones: Copia y pega los comandos que necesites, ajustando los parámetros

#### `microgel/scripts/aaa-otro/ejemplos_uso.sh`
> Ejemplos de uso de los scripts de graficación del campo eléctrico
> Instrucciones: Copia y pega los comandos que necesites, ajustando los parámetros

#### `microgel/scripts/aaa-otro/extraer_posicion.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-nciclos`: number of cycles; `-ninicio`: Number of cycle to start from; `-npaso`: cycles per step; `-o`: output file name

#### `microgel/scripts/aaa-otro/field_potencial.sh`
> Run ludwig microgel simulation with MPI support
> Run './runbg.sh --help' for detailed usage information

#### `microgel/scripts/aaa-otro/ludwig_redo.sh`
> Compile ludwig and copy program to current dir

#### `microgel/scripts/aaa-otro/microgel_redo.sh`
> Compile microgel_poly_cross_density from ../util and copy program to current dir

#### `microgel/scripts/aaa-otro/microgel_run.sh`
> Run microgel_poly_cross_density

#### `microgel/scripts/aaa-otro/organize_by_parameter.sh`
> Script para organizar carpetas de simulaciones por parámetro
> Uso: ./organize_by_parameter.sh <directorio_base> <parametro> [prefijo]
> Ejemplo: ./organize_by_parameter.sh /path/to/simulations rhoel PB_Sn

#### `microgel/scripts/aaa-otro/organize_by_parameter_dryrun.sh`
> Script para VISUALIZAR cómo se organizarían las carpetas de simulaciones por parámetro
> (DRY RUN - no mueve nada, solo muestra qué haría)
> Uso: ./organize_by_parameter_dryrun.sh <directorio_base> <parametro> [prefijo]

#### `microgel/scripts/aaa-otro/plot.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

#### `microgel/scripts/aaa-otro/plot_all_velocities.py`
> Script para combinar plots de velocidad de múltiples subcarpetas stencil
> Lee los archivos datosfluid.csv de cada subcarpeta y genera un plot combinado
>
> Uso:
>     ./plot_all_velocities.py -d <directorio_padre> [--out_dir <directorio_salida>] [--fluid-only]
>
> Ejemplo:
>     ./plot_all_velocities.py -d /path/to/parent/folder --out_dir ./plots
>
> El script busca automáticamente todas las subcarpetas que contengan 'stencil_' en su nombre
> y combina los plots de velocidad de todas ellas en un solo gráfico.

Argumentos (argparse): `-d / --dir`: Directorio padre que contiene las subcarpetas stencil; `--out_dir`: Directorio de salida para los plots; `--fluid-only`: Plotear solo velocidad promedio del fluido

#### `microgel/scripts/aaa-otro/plot_campo_electrico_rhoel.py`
> Script para graficar campo_electrico_medio.csv de múltiples simulaciones
> variando el parámetro rhoel

Lee argumentos de línea de comandos: sys.argv[1]

#### `microgel/scripts/aaa-otro/plot_charge_distribution.py`
> Script para graficar la distribución de carga desde archivos qsi de Ludwig.
>
> Los archivos qsi contienen la densidad de carga en cada punto de la malla.
> Formato típico (2 columnas):
>   - Columna 0: Densidad de carga positiva ρ₊ (cationes, carga +1)
>   - Columna 1: Densidad de carga negativa ρ₋ (aniones, carga -1)
>
> IMPORTANTE: Las densidades son siempre ≥ 0 (representan concentraciones).
> La carga neta en cada nodo es: ρ_net = ρ₊ - ρ₋
>
> Modos de visualización:
> 1. Plano 2D: Muestra la distribución de carga en un plano específico (XY, XZ, o YZ)
> … (11 líneas más en el archivo)

Argumentos (argparse): `-f / --file`: Archivo qsi a leer; `-s / --size`: Tamaño de la malla (nx ny nz); `-m / --mode`: Modo: plane (2D), line (1D), plane3d (superficie 3D), stats (solo estadísticas); `-p / --plane`: Plano a graficar (obligatorio para mode=plane/plane3d); `-pos / --position`: Posición del plano (índice). Por defecto: centro; `--start`: Punto inicial de la línea (obligatorio para mode=line); `--end`: Punto final de la línea (obligatorio para mode=line); `--num-points`: Número de puntos a lo largo de la línea (por defecto: 100); `-c / --component`: Componente a graficar: net (carga neta), species0/1 (especie individual), total (densidad total), all (todas en modo line). Default: net; `--colormap`: Mapa de colores para modo plane3d (default: RdBu_r); `--elevation`: Ángulo de elevación de la vista 3D (default: 30); `--azimuth`: Ángulo azimutal de la vista 3D (default: -60); `--num-species`: Número de especies iónicas (default: 2); `-o / --output`: Archivo de salida (por defecto: mostrar en pantalla)

#### `microgel/scripts/aaa-otro/plot_electric_field.py`
> Script para graficar el campo eléctrico desde archivos psi de Ludwig.
>
> El campo eléctrico se calcula como E = -∇ψ, donde ψ es el potencial eléctrico.
>
> Modos de visualización:
> 1. Plano 2D: Muestra el campo eléctrico en un plano específico (XY, XZ, o YZ)
> 2. Línea 1D: Muestra el campo eléctrico a lo largo de una línea en el espacio
> 3. Plano 3D: Muestra superficie 3D del potencial o campo eléctrico en un plano
>
> SISTEMA DE COORDENADAS:
> - Los nodos de la malla están centrados en: (0.5, 1.5, 2.5, ..., nx-0.5)
> - Los coloides en colloids-*.csv tienen coordenadas continuas (ej: 16.0, 16.0, 16.0)
> … (4 líneas más en el archivo)

Argumentos (argparse): `-f / --file`: Archivo psi a leer; `-s / --size`: Tamaño de la malla (nx ny nz); `-m / --mode`: Modo de graficación: plane (plano 2D), line (línea 1D), o plane3d (superficie 3D); `-p / --plane`: Plano a graficar (obligatorio para mode=plane); `-pos / --position`: Posición del plano (índice). Por defecto: centro; `-v / --vectors`: Mostrar vectores del campo (solo para mode=plane); `--vector-stride`: Espaciado entre vectores (por defecto: 2); `--start`: Punto inicial de la línea (obligatorio para mode=line); `--end`: Punto final de la línea (obligatorio para mode=line); `--num-points`: Número de puntos a lo largo de la línea (por defecto: 100); `-c / --component`: Componente a graficar. Plane/Plane3D: magnitude/x/y/z/psi (default: magnitude). Line: all/magnitude/x/y/z/psi (default: magnitude). Usar "all" en modo line muestra todas las componentes juntas; `--colormap`: Mapa de colores para modo plane3d (default: viridis). Opciones: viridis, plasma, inferno, magma, RdBu_r, coolwarm, seismic; `--elevation`: Ángulo de elevación de la vista 3D en grados (default: 30); `--azimuth`: Ángulo azimutal de la vista 3D en grados (default: -60); `--external-field`: Campo eléctrico externo constante a restar (Ex Ey Ez). Útil para visualizar solo el campo generado por cargas, sin el campo aplicado externamente.; `-o / --output`: Archivo de salida (por defecto: mostrar en pantalla)

#### `microgel/scripts/aaa-otro/plot_electric_field_peskin.py`
> Script para graficar el campo eléctrico usando interpolación con delta de Peskin.
>
> Este script usa el MISMO esquema de interpolación que Ludwig (delta de Peskin)
> en lugar de diferencias finitas estándar.
>
> El delta de Peskin es una función de interpolación suave que se usa en el método
> de frontera inmersa (Immersed Boundary Method). En Ludwig se usa para interpolar
> campos entre la posición de las partículas y la malla de lattice.
>
> Función delta de Peskin:
> - Para |r| <= 1.0:
>   δ(r) = (1/8) * [3 - 2|r| + sqrt(1 + 4|r| - 4|r|²)]
> … (10 líneas más en el archivo)

Argumentos (argparse): `-f / --file`: Archivo psi a leer; `-s / --size`: Tamaño de la malla (nx ny nz); `-m / --mode`: Modo de graficación: plane (plano 2D), line (línea 1D), o plane3d (superficie 3D); `-p / --plane`: Plano a graficar (obligatorio para mode=plane); `-pos / --position`: Posición del plano (índice). Por defecto: centro; `-v / --vectors`: Mostrar vectores del campo (solo para mode=plane); `--vector-stride`: Espaciado entre vectores (por defecto: 2); `--start`: Punto inicial de la línea (obligatorio para mode=line); `--end`: Punto final de la línea (obligatorio para mode=line); `--num-points`: Número de puntos a lo largo de la línea (por defecto: 100); `-c / --component`: Componente a graficar. Plane/Plane3D: magnitude/x/y/z/psi (default: magnitude). Line: all/magnitude/x/y/z/psi (default: magnitude). Usar "all" en modo line muestra todas las componentes juntas; `--colormap`: Mapa de colores para modo plane3d (default: viridis). Opciones: viridis, plasma, inferno, magma, RdBu_r, coolwarm, seismic; `--elevation`: Ángulo de elevación de la vista 3D en grados (default: 30); `--azimuth`: Ángulo azimutal de la vista 3D en grados (default: -60); `--external-field`: Campo eléctrico externo constante a restar (Ex Ey Ez). Útil para visualizar solo el campo generado por cargas, sin el campo aplicado externamente.; `-o / --output`: Archivo de salida (por defecto: mostrar en pantalla)

#### `microgel/scripts/aaa-otro/plot_field_comparizon.sh`
> Run multiple plots to compare teoretical Electric field to simulations

#### `microgel/scripts/aaa-otro/plot_isdw.py`
> Plot the d_isdw (Inverse Squared Distance Weighting) interpolation weights
> This script visualizes how the ISDW function distributes weights among neighboring nodes

#### `microgel/scripts/aaa-otro/plot_peskin.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

#### `microgel/scripts/aaa-otro/plot_velocity_field.py`
> Script para visualizar campos de velocidad de archivos vel-* de Ludwig
> Permite graficar en planos xy, xz, yz o cualquier plano definido por el usuario

Argumentos (argparse): `-f / --file`: Archivo de velocidad específico; `-d / --directory`: Directorio con archivos de velocidad; `-t / --timestep`: Número de paso temporal (requerido si se usa -d); `--plane`: Plano a graficar (default: all); `--index`: Índice del plano (default: centro de la malla); `--skip`: Factor de submuestreo para flechas (default: 1); `--scale`: Escala para las flechas (default: auto); `--colormap`: Mapa de colores (default: viridis); `--magnitude-only`: Graficar solo la magnitud en plano sin vectores; `--single-plot`: Graficar solo un panel con magnitud y vectores (en lugar de 2 paneles); `-o / --output`: Directorio de salida (default: .); `--prefix`: Prefijo para archivos de salida (default: vel); `--show`: Mostrar gráficos interactivamente

#### `microgel/scripts/aaa-otro/plotdatos.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-i`: input file name; `-t / --tipo`: type of data to plot: 'posicion' or 'velocidad' (default: posicion); `-o / --output`: output file name (default: auto-generated based on type)

#### `microgel/scripts/aaa-otro/plotvariacioncarga.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-i`: output file name

#### `microgel/scripts/aaa-otro/plotvel.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

Argumentos (argparse): `-i`: input file name; `--out_dir`: Output directory; `--fluid-only`: plot only fluid average velocity

#### `microgel/scripts/aaa-otro/plotvelfluid.py`
**SIN DOCUMENTACIÓN EN EL ARCHIVO.**

#### `microgel/scripts/aaa-otro/runbg copy.sh`
> Run ludwig microgel simulation with MPI support
> Run './runbg.sh --help' for detailed usage information

#### `microgel/scripts/aaa-otro/runbg.sh`
> Run ludwig microgel simulation with MPI support
> Run './runbg.sh --help' for detailed usage information

#### `microgel/scripts/aaa-otro/runbg_old.sh`
> Run ludwig microgel simulation with MPI support
> Input parameters:
> -n  Nsteps       : Number of steps to calculate
> -i  Ninicio      : Initial step number
> -p  paso         : Step interval for output
> -l  ladox        : Frame size in the X direction
> -y  ladoyz       : Frame size in the Y and Z directions
> -d  del          : Flag to run script that deletes files
> -v  viscosidad   : Viscosity value
> -e  energy       : Free energy model to use: 'fe_electro' or 'none'
> -f  electric_e0  : External electric field: Ex_Ey_Ez
> -m  mpi_procs    : Number of MPI processes (default: 1)
> … (5 líneas más en el archivo)

#### `microgel/scripts/aaa-otro/runmulti copy.sh`
> Run multiple Ludwig simulations with varying colloid parameters
> This script modifies config.cds.init.001-001 file to change charge, position,
> and al parameter for each simulation, then calls runbg.sh

#### `microgel/scripts/aaa-otro/runmulti.sh`
> Run multiple Ludwig simulations with varying colloid parameters
> This script modifies config.cds.init.001-001 file to change charge, position,
> and al parameter for each simulation, then calls runbg.sh

#### `microgel/scripts/aaa-otro/runplot.sh`
> Run ludwig microgel simulation plots
> Input parameters:
> -n  Nsteps       : Number of steps to calculate
> -i  Ninicio      : Initial step number
> -p  paso         : Step interval for output
> -l  ladox        : Frame size in the X direction
> -y  ladoyz       : Frame size in the Y and Z directions
> -d  del          : Flag to run script that deletes files
> -v  viscosidad   : Viscosity value
> -e  energy       : Free energy model to use: 'fe_electro' or 'none'
> -f  electric_e0  : Excternal electric field: Ex_Ey_Ez
> -q  fluid_only   : Plot only fluid average velocity (y/n)
> … (5 líneas más en el archivo)

#### `microgel/scripts/aaa-otro/verificar_simetria_psi.py`
> Script para verificar la simetría del campo eléctrico con respecto a la posición de la partícula.
>
> Este script:
> 1. Lee archivos psi de Ludwig
> 2. Extrae la posición de la partícula del archivo colloids
> 3. Calcula el campo eléctrico E = -∇ψ
> 4. Verifica la simetría radial y en diferentes planos
> 5. Compara con modelo teórico de Poisson-Boltzmann (Debye-Hückel)
> 6. Genera gráficos comparativos y métricas cuantitativas
>
> Autor: Script para análisis de simulaciones Ludwig

Argumentos (argparse): `-f / --file`: Archivo psi a analizar; `-s / --size`: Tamaño de la malla (nx ny nz); `-c / --colloid-file`: Archivo config.cds*.001-001 o colloids*.csv para leer posición de la partícula; `-p / --position`: Posición de la partícula (x y z); `--max-radius`: Radio máximo para análisis radial (default: 15.0); `--num-radii`: Número de radios a analizar (default: 10); `--shell-thickness`: Grosor de las cáscaras esféricas para muestreo radial (default: 0.1); `--planes`: Planos a analizar (default: xy xz yz); `--log-scale`: Escala logarítmica: "y" (solo eje Y), "xy" (ambos ejes log-log), "none" (lineal); `--radial-only`: Solo realizar análisis radial (salta planos, componentes y métricas). Solo guarda el gráfico radial.; `--save-csv`: Guardar perfiles radiales en archivos CSV (solo se aplica si se especifica -o); `--psi-offset-factor`: Factor para calcular offset del potencial en escala log: offset = -psi_min * factor (default: 1.1); `--compare-theory`: Comparar con modelo teórico de Poisson-Boltzmann (Debye-Hückel); `--charge`: Carga de la partícula (default: 1.0); `--epsilon`: Permitividad relativa del medio (default: 1.0); `--kt`: Energía térmica k_B*T (default: 1.0); `--kappa`: Parámetro de Debye κ = 1/λ_D. κ=0 usa Coulomb (vacío), κ>0 usa Debye-Hückel (electrolito). (default: 0.0); `-o / --output-prefix`: Prefijo para archivos de salida (default: mostrar en pantalla)


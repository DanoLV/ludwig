# !/bin/bash
# ------------------------------------------------------------------------------------
# Run ludwig microgel simulation
# ------------------------------------------------------------------------------------
# Input parameters:
#   -n  Nsteps       : Number of steps to calculate
#   -i  Ninicio      : Initial step number
#   -p  paso         : Step interval for output
#   -l  ladox        : Frame size in the X direction
#   -y  ladoyz       : Frame size in the Y and Z directions
#   -d  del          : Flag to run script that deletes files
#   -v  viscosidad   : Viscosity value
#   -e  energy       : Free energy model to use: 'fe_electro' or 'none'
#   -f  electric_e0  : Excternal electric field: Ex_Ey_Ez
# ------------------------------------------------------------------------------------
# Example:
# ./runbg.sh -n 1000000 -i 0 -p 500 -t 300 -c 10 -l 26 -y 26 -v 0.5 -e fe_electro -g petsc -f 0.0_0.0_0.0 -s y -o outdir
# ./coloideacsv.sh -n 372000 -i 000000 -p 500
# ./calculosvel.py -nciclos 372000 -npaso 500 -o datos.csv
# ./runplot.sh  -n 610000 -i 316750 -p 250  -s y >> outputplot.txt 2>&1 
# ------------------------------------------------------------------------------------

# PRUEBA DE CLAUDE---------------------------------------------------------------------------------
# PRUEBA DE CLAUDE---------------------------------------------------------------------------------
# PRUEBA DE CLAUDE---------------------------------------------------------------------------------
# PRUEBA DE CLAUDE---------------------------------------------------------------------------------

#!/bin/bash
#------------------------------------------------------------------------------------
# Run ludwig microgel simulation with MPI support
#------------------------------------------------------------------------------------
# Input parameters:
#   -n  Nsteps       : Number of steps to calculate
#   -i  Ninicio      : Initial step number
#   -p  paso         : Step interval for output
#   -l  ladox        : Frame size in the X direction
#   -y  ladoyz       : Frame size in the Y and Z directions
#   -d  del          : Flag to run script that deletes files
#   -v  viscosidad   : Viscosity value
#   -e  energy       : Free energy model to use: 'fe_electro' or 'none'
#   -f  electric_e0  : External electric field: Ex_Ey_Ez
#   -m  mpi_procs    : Number of MPI processes (default: 1)
#   -r  mpi_grid     : MPI grid decomposition: NX_NY_NZ (e.g., 2_2_1)
#------------------------------------------------------------------------------------
# Example with MPI:
# ./runbg.sh -n 1000000 -i 0 -p 500 -t 300 -c 10 -l 26 -y 26 -v 0.5 -e fe_electro \
#            -g petsc -f 0.0_0.0_0.0 -s y -o outdir -m 4 -r 2_2_1
#------------------------------------------------------------------------------------

clear

# Default values for MPI
mpi_procs=1
mpi_grid=""

while getopts "a:n:i:p:l:d:y:v:e:f:g:s:t:c:o:m:r:" flag
do
    case "${flag}" in
        a) fconfig=${OPTARG};;      # Config out frequency
        n) Nsteps=${OPTARG};;       # Number of steps to calculate
        i) Ninicio=${OPTARG};;      # Initial step number
        p) paso=${OPTARG};;         # Delta steps for output
        t) deltat=${OPTARG};;       # time interval
        l) ladox=${OPTARG};;        # Frame size in X direction
        y) ladoyz=${OPTARG};;       # Frame size in Y and Z directions
        d) del=${OPTARG};;          # Run script to delete files
        v) viscosidad=${OPTARG};;   # Viscosity
        e) energy=${OPTARG};;       # free_energy: fe_electro/none
        f) electric_e0=${OPTARG};;  # electric_e0 
        g) solver=${OPTARG};;       # electrokinetics_solver_type: petsc / sor
        s) single=${OPTARG};;       # Single monomer
        c) cores=${OPTARG};;        # Num of cores (threads per process)
        o) dir=${OPTARG};;          # Output dir
        m) mpi_procs=${OPTARG};;    # Number of MPI processes
        r) mpi_grid=${OPTARG};;     # MPI grid decomposition
    esac
done

# Check output dir is specified
if [ -z "$dir" ]; then
        echo 'Missing output directory -o' >&2
        exit 1
fi

# Check mandatory options
if [[ -z "$Ninicio" || -z "$Nsteps" || -z "$paso" ]]; then
        echo 'Missing mandatory input line parameters' >&2
        exit 1
fi

# outputdir
RESULTS_DIR="$dir"

# Check if simulation starts from a previous one
if [ "$Ninicio" -gt 0 ]; then
    if [ ! -d "$RESULTS_DIR" ]; then # Check output dir exists
        echo "Error: 'output' directory does not exist in the current working directory."
        exit 1
    fi

elif [ "$Ninicio" -eq 0  ]; then 

    if [ -d "$RESULTS_DIR" ]; then # Check if there is already output dir
    echo "Error: 'output' directory already exists in the current working directory."
    exit 1
    fi

    # Create outputdir
    mkdir -p $RESULTS_DIR
    cp config.cds.init.001-001 $RESULTS_DIR
    
    # If using MPI, copy initial config for each process
    if [ "$mpi_procs" -gt 1 ]; then
        for ((proc=1; proc<=mpi_procs; proc++)); do
            proc_str=$(printf "%03d" $proc)
            cp config.cds.init.001-001 $RESULTS_DIR/config.cds.init.001-${proc_str}
        done
    fi

else 
    echo "Error: 'Ninicio' is not a valid positive integer number"
    exit 1
fi

# Copy files to execute simulation
cp input $RESULTS_DIR
cp Ludwig.exe $RESULTS_DIR
cp del.sh $RESULTS_DIR

# Copy files to execute plot
cp runplot.sh $RESULTS_DIR
cp extract_colloids $RESULTS_DIR
cp coloideacsv.sh $RESULTS_DIR
cp calculosvel.py $RESULTS_DIR
cp calculosvelfluid.py $RESULTS_DIR
cp plotvel.py $RESULTS_DIR

# Provide access to plot
chmod +x $RESULTS_DIR/runplot.sh
chmod +x $RESULTS_DIR/extract_colloids
chmod +x $RESULTS_DIR/calculosvel.py
chmod +x $RESULTS_DIR/calculosvelfluid.py
chmod +x $RESULTS_DIR/plotvel.py

# Change to result dir
cd $RESULTS_DIR/

#Delete files from previous runs
if [ "$del" == "y" ]; then
    ./del.sh
fi

# Config for step 0
cp config.cds.init.001-001 config.cds00000000.001-001

echo "Inicia simulacion con MPI ($mpi_procs procesos):"

# Change parameters in input file
sed -i -e "/freq_config/c\freq_config $fconfig" input
sed -i -e "/N_start/c\N_start $Ninicio" input
sed -i -e "/N_cycles/c\N_cycles $Nsteps" input
sed -i -e "/^viscosity /c\viscosity $viscosidad" input
sed -i -e "/^viscosity_bulk/c\viscosity_bulk $viscosidad" input
sed -i -e "/^free_energy/c\free_energy $energy" input
sed -i -e "/^electrokinetics_solver_type/c\electrokinetics_solver_type $solver" input 
sed -i -e "/^electric_e0/c\electric_e0 $electric_e0" input
sed -i -e "/colloid_io_freq/c\colloid_io_freq $paso" input
sed -i -e "/vel_io_freq/c\vel_io_freq $paso" input
# sed -i -e "/psi_io_freq/c\psi_io_freq $paso" input
sed -i -e "/size/c\size $ladox\_$ladoyz\_$ladoyz" input

# Configure MPI grid if specified
if [ -n "$mpi_grid" ]; then
    # Check if grid_decomposition line exists in input file
    if grep -q "^grid" input; then
        sed -i -e "/^grid/c\grid $mpi_grid" input
    else
        echo "grid $mpi_grid" >> input
    fi
fi

# Total steps of simulation
NT=$((Nsteps + Ninicio))

# Set number of OpenMP threads per MPI process
if [ -z "$cores" ]; then
    cores=1
fi
export OMP_NUM_THREADS=$cores

# --- 1. Ejecutar tarea principal en segundo plano ---
echo "Ejecutando tarea de fondo con $mpi_procs procesos MPI y $cores threads por proceso..."

# Determine MPI command (mpirun or mpiexec)
if command -v mpirun &> /dev/null; then
    MPI_CMD="mpirun"
elif command -v mpiexec &> /dev/null; then
    MPI_CMD="mpiexec"
else
    echo "Error: No MPI command found (mpirun or mpiexec)"
    exit 1
fi

# Ludwig run with MPI on background
if [ "$mpi_procs" -gt 1 ]; then
    $MPI_CMD -np $mpi_procs ./Ludwig.exe >> output.txt 2>&1 &
else
    # Single process without MPI
    ./Ludwig.exe >> output.txt 2>&1 &
fi

PID_BG=$!
echo "PID = $PID_BG"

# --- 2. Ejecutar tarea periódica mientras la principal corre ---
count=0
ni=$Ninicio

while kill -0 $PID_BG 2>/dev/null; do

    echo "[INFO] Plot inicia en $(date)"

    # Plot
    count=$(find . -maxdepth 1 -type f -name "config.cds*" | wc -l)
    
    # Adjust count for MPI processes
    if [ "$mpi_procs" -gt 1 ]; then
        # count=$((count / mpi_procs))
        count=$((count - mpi_procs + 1))
    fi
    
    count=$((paso*(count-2)-ni))
    
    echo "Ni=$ni"

    if [ "$count" -lt 0 ]; then
        sleep $deltat
        continue
    fi

    ./runplot.sh   \
                -n "$count"  \
                -i "$ni"   \
                -p "$paso"   \
                -s "$single" >> outputplot.txt 2>&1 

    ni=$((ni+count+paso))

    sleep $deltat

done

echo "[INFO] Tarea periódica en $(date)"

# Plot final
echo "Plot inicia"
count=$(find . -maxdepth 1 -type f -name "config.cds*" | wc -l)

# Adjust count for MPI processes
if [ "$mpi_procs" -gt 1 ]; then
    count=$((count / mpi_procs))
fi

count=$((paso*(count-2)-ni))

./runplot.sh   \
                -n "$count"  \
                -i "$ni"   \
                -p "$paso"   \
                -s "$single" >> outputplot.txt 2>&1 

echo "Plot termino"
echo "Simulacion finalizada con $mpi_procs procesos MPI"

# clear

# # set number of CPU cores to use
# export OMP_NUM_THREADS=1

# while getopts "a:n:i:p:l:d:y:v:e:f:g:s:t:c:o:" flag
# do
#     case "${flag}" in
#         a) fconfig=${OPTARG};;      # Config out frequency
#         n) Nsteps=${OPTARG};;       # Number of steps to calculate
#         i) Ninicio=${OPTARG};;      # Initial step number
#         p) paso=${OPTARG};;         # Delta steps for output
#         t) deltat=${OPTARG};;       # time interval
#         l) ladox=${OPTARG};;        # Frame size in X direction
#         y) ladoyz=${OPTARG};;       # Frame size in Y and Z directions
#         d) del=${OPTARG};;          # Run script to delete files
#         v) viscosidad=${OPTARG};;   # Viscosity
#         e) energy=${OPTARG};;       # free_energy: fe_electro/none
#         f) electric_e0=${OPTARG};;  # electric_e0 
#         g) solver=${OPTARG};;       # electrokinetics_solver_type: petsc / sor
#         s) single=${OPTARG};;       # Single monomer
#         c) cores=${OPTARG};;        # Num of cores
#         o) dir=${OPTARG};;          # Output dir
#     esac
# done

# # Check output dir is specified
# if [ -z "$dir" ]; then
#         echo 'Missing output directory -o' >&2
#         exit 1
# fi

# # Check mandatory options
# if [[ -z "$Ninicio" || -z "$Nsteps" || -z "$paso" ]]; then
#         echo 'Missing mandatory input line parameters' >&2
#         exit 1
# fi

# # outputdir
# RESULTS_DIR="$dir"

# # Check if simulation startr from a previous one
# if [ "$Ninicio" -gt 0 ]; then
#     if [ ! -d "$RESULTS_DIR" ]; then # Check output dir exists
#         echo "Error: 'output' directory does not exists in the current working directory."
#         exit 1
#     fi

# elif [ "$Ninicio" -eq 0  ]; then 

#     if [ -d "$RESULTS_DIR" ]; then # Check if there is already output dir
#     echo "Error: 'output' directory already exists in the current working directory."
#     exit 1
#     fi

#     # Create outputdir
#     mkdir -p $RESULTS_DIR
#     cp config.cds.init.001-001 $RESULTS_DIR
#     # # Config for step 0
#     # cp config.cds.init.001-001 $RESULTS_DIR/config.cds00000000.001-001

# else 
#     echo "Error: 'Ninicio' is not a valid positive integer number"
#     exit 1
# fi
# # Copy files to execute simulation
# cp input $RESULTS_DIR
# cp Ludwig.exe $RESULTS_DIR
# cp del.sh $RESULTS_DIR

# # Copy files to execute plot
# cp runplot.sh $RESULTS_DIR
# cp extract_colloids $RESULTS_DIR
# cp coloideacsv.sh $RESULTS_DIR
# cp calculosvel.py $RESULTS_DIR
# cp calculosvelfluid.py $RESULTS_DIR
# cp plotvel.py $RESULTS_DIR

# # Provide access tu plot
# chmod +x $RESULTS_DIR/runplot.sh
# chmod +x $RESULTS_DIR/extract_colloids
# chmod +x $RESULTS_DIR/calculosvel.py
# chmod +x $RESULTS_DIR/calculosvelfluid.py
# chmod +x $RESULTS_DIR/plotvel.py

# # Change to result dir
# cd $RESULTS_DIR/

# #Delete files from previus runs
# if [ "$del" == "y" ]; then
#     ./del.sh
# fi

# # Config for step 0
# cp config.cds.init.001-001 config.cds00000000.001-001

# echo "Inicia simulacion:"

# # Change parameters in input file

# sed -i -e "/freq_config/c\freq_config $fconfig" input
# sed -i -e "/N_start/c\N_start $Ninicio" input
# sed -i -e "/N_cycles/c\N_cycles $Nsteps" input
# sed -i -e "/^viscosity /c\viscosity $viscosidad" input
# sed -i -e "/^viscosity_bulk/c\viscosity_bulk $viscosidad" input
# sed -i -e "/^free_energy/c\free_energy $energy" input
# sed -i -e "/^electrokinetics_solver_type/c\electrokinetics_solver_type $solver" input 
# sed -i -e "/^electric_e0/c\electric_e0 $electric_e0" input
# sed -i -e "/colloid_io_freq/c\colloid_io_freq $paso" input
# sed -i -e "/vel_io_freq/c\vel_io_freq $paso" input
# sed -i -e "/psi_io_freq/c\psi_io_freq $paso" input
# sed -i -e "/size/c\size $ladox\_$ladoyz\_$ladoyz" input

# # Total steps of simulation
# NT=$((Nsteps + Ninicio))

# # set number of CPU cores to use
# export OMP_NUM_THREADS=$cores

# # --- 1. Ejecutar tarea principal en segundo plano ---
# echo "Ejecutando tarea de fondo..."

# # Ludwing run on background
# ./Ludwig.exe >> output.txt 2>&1 &
# PID_BG=$!
# echo "PID = $PID_BG"

# # --- 2. Ejecutar tarea periódica mientras la principal corre ---
# count=0
# ni=$Ninicio

# while kill -0 $PID_BG 2>/dev/null; do

#     echo "[INFO] Plot inicia en $(date)"

#     # Plot
#     count=$(find . -maxdepth 1 -type f -name "config.cds*" | wc -l)
#     count=$((paso*(count-2)-ni))
#     # echo "count=$count"
#     echo "Ni=$ni"

#     if [ "$count" -lt 0 ]; then
#         sleep $deltat
#         continue
#     fi

#     ./runplot.sh   \
#                 -n "$count"  \
#                 -i "$ni"   \
#                 -p "$paso"   \
#                 -s "$single" >> outputplot.txt 2>&1 

#     ni=$((ni+count+paso))

#     sleep $deltat

# done

# echo "[INFO] Tarea periódica en $(date)"

# # Plot
# echo "Plot inicia"
# count=$(find . -maxdepth 1 -type f -name "config.cds*" | wc -l)
# count=$((paso*(count-2)-ni))
# # echo "count=$count"

# ./runplot.sh   \
#                 -n "$count"  \
#                 -i "$ni"   \
#                 -p "$paso"   \
#                 -s "$single" >> outputplot.txt 2>&1 

# echo "Plot termino"

# # # Run Ludwig
# # ./Ludwig.exe

# # # Postprocesing - convert data to .cvs files
# # cp config.cds.init.001-001 config.cds00000000.001-001
# # ./coloideacsv.sh -n $NT -i $Ninicio -p $paso

# # if [ "$single" == "y" ]; then
# # # Calculates and plots velocity for a single subgrid monomer
# # ./calculosvel.py -nciclos $NT -npaso $paso -o datos.csv
# # ./plotvel.py
# # # cp config.cds.init.001-001 config.cds00000000.001-001
# # # ./coloideacsv.sh -n 567000 -i 0 -p 500
# # # ./calculosvel.py -nciclos 567000 -npaso 500 -o datos.csv
# # # ./plotvel.py
# # else 
# # # Calculates and plots density, medium bond length and inertia moments for a microgel
# # ./calculos.py -nciclos $NT -npaso $paso -o datos.csv
# # ./plot.py
# # fi
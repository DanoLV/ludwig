#!/bin/bash 

#SBATCH --job-name=Ludwig_run_Dano
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# Printing date and time
date

# Simulation parameters
while getopts "n:i:p:t:l:d:y:v:e:f:s:o:" flag
do
    case "${flag}" in
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
        s) single=${OPTARG};;       # Single monomer
        o) dir=${OPTARG};;          # Output dir
    esac
done

# Check output dir is specified
if [ -z "$dir" ]; then
        echo 'Missing output directory -o' >&2
        exit 1
fi

# Check mandatory options
if [[ -z "$Ninicio" || -z "$Nsteps" || -z "$paso" || -z "$single" ]]; then
        echo 'Missing mandatory input line parameters' >&2
        exit 1
fi

# outputdir
RESULTS_DIR="$SLURM_SUBMIT_DIR/$dir"

# Check if simulation starts from a previous one
if [ "$Ninicio" -gt 0 ]; then
    if [ ! -d "$RESULTS_DIR" ]; then # Check output dir exists
        echo "Error: 'output' directory does not exists in the current working directory."
        exit 1
    fi

    # Copy necesary files to continue simulation from Ninicio
    # Colloid info
    file=$(printf "$RESULTS_DIR/config.cds%08d.001-001" "$Ninicio")  
    cp $file $TMPDIR/
    # Distribution info
    file=$(printf "$RESULTS_DIR/dist-%09d.001-001" "$Ninicio")
    cp $file $TMPDIR/
    # Electric field info
    file=$(printf "$RESULTS_DIR/psi-%09d.001-001" "$Ninicio") 
    cp $file $TMPDIR/
    # electric charge density  info
    file=$(printf "$RESULTS_DIR/qsi-%09d.001-001" "$Ninicio") 
    cp $file $TMPDIR/
    # density info
    file=$(printf "$RESULTS_DIR/rho-%09d.001-001" "$Ninicio") 
    cp $file $TMPDIR/
    # velicity field info
    file=$(printf "$RESULTS_DIR/vel-%09d.001-001" "$Ninicio") 
    cp $file $TMPDIR/

elif [ "$Ninicio" -eq 0  ]; then 

    if [ -d "$RESULTS_DIR" ]; then # Check if there is already output dir
    echo "Error: 'output' directory already exists in the current working directory."
    exit 1
    fi

    # Create outputdir
    mkdir -p $RESULTS_DIR
    cp config.cds.init.001-001 $TMPDIR
    # Config for step 0
    cp config.cds.init.001-001 $TMPDIR/config.cds00000000.001-001

else 
    echo "Error: 'Ninicio' is not a valid positive integer number"
    exit 1
fi

# Load modules
module load OpenMPI/5.0.7-GCC-14.2.0 #OpenMPI/4.1.6-GCC-13.2.0  # OpenMPI/4.1.5-GCC-12.3.0
#module load GCC # /13.2.0
#module load git 

echo "Running on " `hostname`

# Now comes the commands to be executed
# export OMP_NUM_THREADS=20
#export LD_LIBRARY_PATH=/home/gdallava/petsc/arch-linux-c-debug/lib/ 

# Copy files to execute simulation
cp input $TMPDIR
# cp config.cds.init.001-001 $TMPDIR
cp Ludwig.exe $TMPDIR
cp run-local.sh $TMPDIR

# Copy files to execute plot
cp runplot.sh $RESULTS_DIR
cp extract_colloids $RESULTS_DIR
cp coloideacsv.sh $RESULTS_DIR
cp calculosvel.py $RESULTS_DIR
cp plotvel.py $RESULTS_DIR

# Provide access tu plot
chmod +x $RESULTS_DIR/runplot.sh
chmod +x $RESULTS_DIR/extract_colloids
chmod +x $RESULTS_DIR/calculosvel.py
chmod +x $RESULTS_DIR/plotvel.py

# Change to temporal dir
cd $TMPDIR/
echo $TMPDIR
ls -l $TMPDIR

# Provide access tu run
chmod +x Ludwig.exe
chmod +x run-local.sh

echo "Inicia simulacion:"

#-------------------------------------------------------------------------------------------------------------
# Run simulation in background and copying file every time deltat
#-------------------------------------------------------------------------------------------------------------
# Total steps of simulation
NT=$((Nsteps + Ninicio))

# --- 1. Ejecutar tarea principal en segundo plano ---
echo "Ejecutando tarea de fondo..."

# stdbuf -oL bash ./run-local.sh  \
nice -n 0 srun --exclusive -n1 \
        # --output=outputludwig.txt \
        stdbuf -oL -eL \
        bash ./run-local.sh  \
        -n "$Nsteps"  \
        -i "$Ninicio"  \
        -p "$paso"  \
        -l "$ladox"  \
        -y "$ladoyz"  \
        -v "$viscosidad"  \
        -e "$energy"  \
        -f "$electric_e0" >> outputludwig.txt 2>&1 &
PID_BG=$!

# --- 2. Ejecutar tarea periódica mientras la principal corre ---
count=0
ni=$Ninicio

while kill -0 $PID_BG 2>/dev/null; do
    echo "[INFO] Tarea periódica en $(date)"

    # Copy files to output dir and delete from temp folder
    cp -n $TMPDIR/config.cds* $RESULTS_DIR
    ls -t $TMPDIR/config.cds0* | tail -n +2 | xargs -d '\n' rm --

    cp -n $TMPDIR/rho-* $RESULTS_DIR
    ls -t $TMPDIR/rho-* | tail -n +2 | xargs -d '\n' rm --

    cp -n $TMPDIR/qsi-* $RESULTS_DIR
    ls -t $TMPDIR/qsi-* | tail -n +2 | xargs -d '\n' rm --

    cp -n $TMPDIR/psi-* $RESULTS_DIR
    ls -t -n $TMPDIR/psi-* | tail -n +2 | xargs -d '\n' rm --

    cp -n $TMPDIR/dist-* $RESULTS_DIR
    ls -t $TMPDIR/dist-* | tail -n +2 | xargs -d '\n' rm --

    cp -n $TMPDIR/vel-* $RESULTS_DIR
    ls -t $TMPDIR/vel-* | tail -n +2 | xargs -d '\n' rm --

    cp $TMPDIR/output*.txt $RESULTS_DIR

    cp $TMPDIR/input $RESULTS_DIR

    # Plot
    cd $RESULTS_DIR/
    echo "Plot inicia"
    count=$(find . -maxdepth 1 -type f -name "config.cds*" | wc -l)
    count=$((paso*(count-2)-ni))
    echo "count=$count"
    if [ "$count" -lt 0 ]; then
        sleep $deltat
        continue
    fi
    
    stdbuf -oL ./runplot.sh   \
                -n "$count"  \
                -i "$ni"   \
                -p "$paso"   \
                -s "$single" > outputplot.txt 2>&1 

    echo "Plot termino"

    ni=$((ni+count))

    cd $TMPDIR/
    ls -l

    sleep $deltat
done

# --- 3. Ejecutar tarea final ---
echo "[INFO] Tarea principal finalizada. Ejecutando tarea final..."

 # Copy files to output dir and delete from temp folder
cp -n $TMPDIR/config.cds* $RESULTS_DIR
ls -t $TMPDIR/config.cds* | tail -n +2 | xargs -d '\n' rm --

cp -n $TMPDIR/rho-* $RESULTS_DIR
ls -t $TMPDIR/rho-* | tail -n +2 | xargs -d '\n' rm --

cp -n $TMPDIR/qsi-* $RESULTS_DIR
ls -t $TMPDIR/qsi-* | tail -n +2 | xargs -d '\n' rm --

cp -n $TMPDIR/psi-* $RESULTS_DIR
ls -t -n $TMPDIR/psi-* | tail -n +2 | xargs -d '\n' rm --

cp -n $TMPDIR/dist-* $RESULTS_DIR
ls -t $TMPDIR/dist-* | tail -n +2 | xargs -d '\n' rm --

cp -n $TMPDIR/vel-* $RESULTS_DIR
ls -t $TMPDIR/vel-* | tail -n +2 | xargs -d '\n' rm --

cp $TMPDIR/output*.txt $RESULTS_DIR

cp $TMPDIR/input $RESULTS_DIR

# Plot
cd $RESULTS_DIR/
stdbuf -oL ./runplot.sh -n "$Nsteps" -i "$Ninicio" -p "$paso" -s "$single" > outputplot.txt 2>&1
cd $TMPDIR/

cp -n $TMPDIR/input $RESULTS_DIR

module purge

date

echo "Job finished successfully."
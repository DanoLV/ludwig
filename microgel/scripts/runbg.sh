#!/bin/bash
#------------------------------------------------------------------------------------
# Run ludwig microgel simulation with MPI support
#------------------------------------------------------------------------------------
# Run './runbg.sh --help' for detailed usage information
#------------------------------------------------------------------------------------

show_help() {
    cat << EOF
Usage: ./runbg.sh [OPTIONS]

Run Ludwig microgel simulation with MPI support.

MANDATORY OPTIONS:
  -i, --initial-step NUMBER       Initial step number (0 for new simulation)
  -n, --nsteps NUMBER             Number of simulation steps to calculate
  -o, --output-dir DIR            Output directory for results
  -s, --step-interval NUMBER      Step interval for output files

SIMULATION OPTIONS (alphabetically sorted):
  -a, --angle-harmonic PARAMS     Angle harmonic: on_k_theta0 (e.g., 1_1.0e-04_2.0944)
  -b, --bond-harmonic PARAMS      Bond harmonic: on_k_r0 (e.g., 1_5.0e-5_0.5)
  -c, --cores NUMBER              Number of OpenMP threads per MPI process (default: 1)
      --config-file FILE          Path to config file (default: config.cds.init.001-001)
  -d, --delete-files y/n          Delete previous files before starting (default: n)
  -e, --electric-field Ex_Ey_Ez   External electric field components (e.g., 0.0_0.0_0.0)
  -f, --freq-config NUMBER        Configuration output frequency
  -w, --fluctuations 0/1          LB fluctuations: 0=off, 1=on (default: 0)
  -g, --grid-mpi NX_NY_NZ         MPI grid decomposition (e.g., 2_2_1 for 2x2x1 grid)
  -j, --gravity Gx_Gy_Gz          Colloid gravity vector (e.g., 0.0_0.0_-0.000005)
      --input-file FILE           Path to input file (default: input)
  -k, --temperature NUMBER        Temperature (kT) for thermal fluctuations
  -l, --fluid-only y/n            Plot only fluid average velocity (y/n, default: n)
  -m, --mpi-procs NUMBER          Number of MPI processes (default: 1)
  -p, --plot-interval SECONDS     Time interval between periodic plots (default: 300)
  -q, --relaxation-scheme TYPE    LB relaxation scheme: 'M10' or 'BGK' (default: M10)
  -r, --rho NUMBER                Fluid density (rho0)
  -t, --free-energy TYPE          Free energy model: 'fe_electro' or 'none' (default: none)
  -u, --single-monomer y/n        Single monomer mode (y/n, default: n)
  -v, --viscosity NUMBER          Fluid viscosity value
  -x, --size-x NUMBER             Frame size in X direction (grid units)
  -y, --size-y NUMBER             Frame size in Y direction (grid units)
  -z, --size-z NUMBER             Frame size in Z direction (grid units)
      --solver TYPE               Electrokinetics solver: 'petsc' or 'sor' (default: sor)
      --stencil NUMBER            Electrokinetics solver stencil: 7, 19, or 27 (default: 7)
      --epsilon NUMBER            Electrokinetics epsilon value
      --electrokinetics-init TYPE Electrokinetics init mode: 'point_charges', 'uniform', or 'none'
                                  (default: none - line stays commented)
      --rho-el NUMBER             Electrokinetics init rho el value

HELP:
  -h, --help                      Display this help message and exit

EXAMPLES:
  # New simulation with single MPI process:
  ./runbg_new.sh --initial-step 0 --nsteps 1000000 --step-interval 500 \\
                 --output-dir results --size-x 26 --size-y 26 --size-z 26 --viscosity 0.5 \\
                 --free-energy fe_electro --electric-field 0.0_0.0_0.1

  # Simulation with MPI (4 processes on 2x2x1 grid):
  ./runbg_new.sh -i 0 -n 1000000 -s 500 -o results -x 26 -y 26 -z 26 -v 0.5 \\
                 -t fe_electro -e 0.0_0.0_0.1 -m 4 -g 2_2_1 -c 10

  # Continue previous simulation:
  ./runbg_new.sh -i 1000000 -n 500000 -s 500 -o results

  # New simulation with custom input and config files:
  ./runbg_new.sh -i 0 -n 1000000 -s 500 -o results \\
                 --input-file input.custom --config-file config.cds.init.001-001-L32

EOF
}

# Default values
cores=1
del="n"
deltat=60
energy="none"
fluctuations=0
fluid_only="n"
mpi_procs=1
mpi_grid=""
single="n"
solver="petsc"
stencil="7"
relaxation_scheme="M10"
input_file="input"
config_file="config.cds.init.001-001"

# Parse command line arguments
OPTS=$(getopt -o a:b:c:d:e:f:w:g:hj:k:i:l:m:n:o:p:q:r:s:t:u:v:x:y:z: \
              --long angle-harmonic:,bond-harmonic:,cores:,config-file:,delete-files:,electric-field:,freq-config:,fluctuations:,grid-mpi:,help,gravity:,input-file:,temperature:,initial-step:,fluid-only:,mpi-procs:,nsteps:,output-dir:,plot-interval:,relaxation-scheme:,rho:,step-interval:,free-energy:,single-monomer:,viscosity:,size-x:,size-y:,size-z:,solver:,stencil:,epsilon:,electrokinetics-init:,rho-el: \
              -n 'runbg_new.sh' -- "$@")

if [ $? != 0 ]; then
    echo "Error parsing options. Try './runbg_new.sh --help' for more information." >&2
    exit 1
fi

eval set -- "$OPTS"

# Process options (alphabetically sorted by long name)
while true; do
    case "$1" in
        -a|--angle-harmonic)    angle_harmonic="$2"; shift 2 ;;
        -b|--bond-harmonic)     bond_harmonic="$2"; shift 2 ;;
        -c|--cores)             cores="$2"; shift 2 ;;
        --config-file)          config_file="$2"; shift 2 ;;
        -d|--delete-files)      del="$2"; shift 2 ;;
        -e|--electric-field)    electric_e0="$2"; shift 2 ;;
        -w|--fluctuations)      fluctuations="$2"; shift 2 ;;
        -f|--freq-config)       fconfig="$2"; shift 2 ;;
        -t|--free-energy)       energy="$2"; shift 2 ;;
        -l|--fluid-only)        fluid_only="$2"; shift 2 ;;
        -j|--gravity)           gravity="$2"; shift 2 ;;
        -g|--grid-mpi)          mpi_grid="$2"; shift 2 ;;
        -h|--help)              show_help; exit 0 ;;
        -i|--initial-step)      Ninicio="$2"; shift 2 ;;
        --input-file)           input_file="$2"; shift 2 ;;
        -m|--mpi-procs)         mpi_procs="$2"; shift 2 ;;
        -n|--nsteps)            Nsteps="$2"; shift 2 ;;
        -o|--output-dir)        dir="$2"; shift 2 ;;
        -p|--plot-interval)     deltat="$2"; shift 2 ;;
        -q|--relaxation-scheme) relaxation_scheme="$2"; shift 2 ;;
        -r|--rho)               rho="$2"; shift 2 ;;
        -u|--single-monomer)    single="$2"; shift 2 ;;
        -s|--step-interval)     paso="$2"; shift 2 ;;
        -k|--temperature)       temperature="$2"; shift 2 ;;
        -v|--viscosity)         viscosidad="$2"; shift 2 ;;
        -x|--size-x)            ladox="$2"; shift 2 ;;
        -y|--size-y)            ladoy="$2"; shift 2 ;;
        -z|--size-z)            ladoz="$2"; shift 2 ;;
        --solver)               solver="$2"; shift 2 ;;
        --stencil)              stencil="$2"; shift 2 ;;
        --epsilon)              epsilon="$2"; shift 2 ;;
        --electrokinetics-init) electrokinetics_init="$2"; shift 2 ;;
        --rho-el)               rho_el="$2"; shift 2 ;;
        --)                     shift; break ;;
        *)                      echo "Internal error!"; exit 1 ;;
    esac
done

# Check output dir is specified
if [ -z "$dir" ]; then
        echo 'ERROR: Missing output directory. Use -o or --output-dir' >&2
        echo "Run './runbg.sh --help' for usage information." >&2
        exit 1
fi

# Check mandatory options
if [[ -z "$Ninicio" || -z "$Nsteps" || -z "$paso" ]]; then
        echo 'ERROR: Missing mandatory parameters:' >&2
        [ -z "$Ninicio" ] && echo '  - Initial step (-i or --initial-step)' >&2
        [ -z "$Nsteps" ] && echo '  - Number of steps (-n or --nsteps)' >&2
        [ -z "$paso" ] && echo '  - Step interval (-s or --step-interval)' >&2
        echo "Run './runbg_new.sh --help' for usage information." >&2
        exit 1
fi

# Verify input and config files exist
if [ ! -f "$input_file" ]; then
    echo "ERROR: Input file '$input_file' not found" >&2
    exit 1
fi

if [ ! -f "$config_file" ]; then
    echo "ERROR: Config file '$config_file' not found" >&2
    exit 1
fi

# outputdir
RESULTS_DIR="$dir"

# Check if simulation starts from a previous one
if [ "$Ninicio" -gt 0 ]; then
    if [ ! -d "$RESULTS_DIR" ]; then # Check output dir exists
        echo "ERROR: Output directory '$RESULTS_DIR' does not exist." >&2
        echo "       For continuing simulations (initial-step > 0), the directory must exist." >&2
        exit 1
    fi

elif [ "$Ninicio" -eq 0  ]; then

    if [ -d "$RESULTS_DIR" ]; then # Check if there is already output dir
        echo "ERROR: Output directory '$RESULTS_DIR' already exists." >&2
        echo "       For new simulations (initial-step = 0), the directory must not exist." >&2
        exit 1
    fi

    # Create outputdir
    mkdir -p $RESULTS_DIR
    cp "$config_file" $RESULTS_DIR/config.cds.init.001-001

    # If using MPI, copy initial config for each process
    if [ "$mpi_procs" -gt 1 ]; then
        for ((proc=1; proc<=mpi_procs; proc++)); do
            proc_str=$(printf "%03d" $proc)
            cp "$config_file" $RESULTS_DIR/config.cds.init.001-${proc_str}
        done
    fi

else
    echo "ERROR: Initial step must be a valid non-negative integer" >&2
    exit 1
fi

# Copy files to execute simulation
cp "$input_file" $RESULTS_DIR/input
cp Ludwig.exe $RESULTS_DIR
# Copy PETSc options file if present
[ -f ".petscrc" ] && cp .petscrc $RESULTS_DIR/.petscrc
# cp efield_self_*.bin $RESULTS_DIR/
# cp del.sh $RESULTS_DIR

# Create graphics subdirectory
mkdir -p $RESULTS_DIR/scripts
mkdir -p $RESULTS_DIR/plots
mkdir -p $RESULTS_DIR/proceced_data
mkdir -p $RESULTS_DIR/colloid_data
mkdir -p $RESULTS_DIR/logs

# Copy files to execute plot
cp runplot.sh $RESULTS_DIR/scripts
cp extract_colloids $RESULTS_DIR/scripts
cp coloideacsv.sh $RESULTS_DIR/scripts

# Scripts that process raw data (need access to colloids-*.csv and vel-*) stay in root
cp calculosvel.py $RESULTS_DIR/scripts
cp calculosvelfluid.py $RESULTS_DIR/scripts
cp calculosvelfluidonly.py $RESULTS_DIR/scripts
cp extraer_posicion.py $RESULTS_DIR/scripts
cp calculos.py $RESULTS_DIR/scripts
cp calc_mui.py $RESULTS_DIR/scripts

# Plotting scripts (only read processed CSVs) go to plot subdirectory
cp plotvel.py $RESULTS_DIR/scripts
cp plot.py $RESULTS_DIR/scripts
cp plot_velocity_field.py $RESULTS_DIR/scripts
cp plot_charge_distribution.py $RESULTS_DIR/scripts
cp plot_electric_field.py $RESULTS_DIR/scripts
# cp plot_electric_field_peskin.py $RESULTS_DIR/scripts
cp plotdatos.py $RESULTS_DIR/scripts
# cp batch_plot_electric_field.py $RESULTS_DIR/scripts
cp compare_field_theory.py $RESULTS_DIR/scripts
# cp compare_field_theory_peskin.py $RESULTS_DIR/scripts

cp plot_*ewald.py $RESULTS_DIR/scripts
cp plot_mui*.py $RESULTS_DIR/scripts

# Provide access to plot scripts
chmod +x $RESULTS_DIR/scripts/runplot.sh
chmod +x $RESULTS_DIR/scripts/extract_colloids
chmod +x $RESULTS_DIR/scripts/calculosvel.py
chmod +x $RESULTS_DIR/scripts/calculosvelfluid.py
chmod +x $RESULTS_DIR/scripts/extraer_posicion.py
chmod +x $RESULTS_DIR/scripts/calculos.py
chmod +x $RESULTS_DIR/scripts/plotvel.py
chmod +x $RESULTS_DIR/scripts/plot.py
chmod +x $RESULTS_DIR/scripts/plot_velocity_field.py
chmod +x $RESULTS_DIR/scripts/plot_charge_distribution.py
chmod +x $RESULTS_DIR/scripts/plot_electric_field.py
chmod +x $RESULTS_DIR/scripts/plotdatos.py
chmod +x $RESULTS_DIR/scripts/compare_field_theory.py

chmod +x $RESULTS_DIR/scripts/plot_*ewald.py

# Change to result dir
cd $RESULTS_DIR/

# Delete files from previous runs
if [ "$del" == "y" ]; then
    ./del.sh
fi

# Config for step 0
cp config.cds.init.001-001 config.cds00000000.001-001

echo "=========================================="
echo "Starting Ludwig simulation with MPI"
echo "=========================================="
echo "MPI processes:     $mpi_procs"
echo "OpenMP threads:    $cores"
echo "Output directory:  $RESULTS_DIR"
echo "Input file:        $input_file"
echo "Config file:       $config_file"
echo "Initial step:      $Ninicio"
echo "Total steps:       $Nsteps"
echo "Step interval:     $paso"
echo "=========================================="

# Change parameters in input file
[ -n "$fconfig" ] && sed -i -e "/freq_config/c\freq_config $fconfig" input
sed -i -e "/N_start/c\N_start $Ninicio" input
sed -i -e "/N_cycles/c\N_cycles $Nsteps" input
[ -n "$rho" ] && sed -i -e "/^fluid_rho0 /c\fluid_rho0 $rho" input
[ -n "$viscosidad" ] && sed -i -e "/^viscosity /c\viscosity $viscosidad" input
[ -n "$viscosidad" ] && sed -i -e "/^viscosity_bulk/c\viscosity_bulk $viscosidad" input
sed -i -e "/^free_energy/c\free_energy $energy" input
[ -n "$solver" ] && sed -i -e "/^electrokinetics_solver_type/c\electrokinetics_solver_type $solver" input
[ -n "$stencil" ] && sed -i -e "/^electrokinetics_solver_stencil/c\electrokinetics_solver_stencil $stencil" input
[ -n "$electric_e0" ] && sed -i -e "/^electric_e0/c\electric_e0 $electric_e0" input
sed -i -e "/colloid_io_freq/c\colloid_io_freq $paso" input
sed -i -e "/vel_io_freq/c\vel_io_freq $paso" input
[ -n "$ladox" ] && [ -n "$ladoy" ] && [ -n "$ladoz" ] && sed -i -e "/size/c\size $ladox\_$ladoy\_$ladoz" input

# Electrokinetics epsilon parameter
if [ -n "$epsilon" ]; then
    if grep -q "^electrokinetics_epsilon" input; then
        sed -i -e "/^electrokinetics_epsilon/c\electrokinetics_epsilon $epsilon" input
    else
        sed -i -e "/^# electrokinetics_epsilon/a electrokinetics_epsilon $epsilon" input
    fi
fi

# Electrokinetics init parameter
if [ -n "$electrokinetics_init" ]; then
    if [ "$electrokinetics_init" = "none" ]; then
        # Comment out the line if it exists uncommented
        if grep -q "^electrokinetics_init " input; then
            sed -i -e "s/^electrokinetics_init /# electrokinetics_init /" input
        fi
    else
        # Uncomment and set the value
        if grep -q "^electrokinetics_init " input; then
            sed -i -e "/^electrokinetics_init /c\electrokinetics_init $electrokinetics_init" input
        elif grep -q "^# electrokinetics_init " input; then
            sed -i -e "0,/^# electrokinetics_init /s/^# electrokinetics_init.*/electrokinetics_init $electrokinetics_init/" input
        else
            # If line doesn't exist at all, add it after the commented lines
            sed -i -e "/^# electric_e0 /a electrokinetics_init $electrokinetics_init" input
        fi
    fi
else
    # If not specified, comment out the line if it exists uncommented
    if grep -q "^electrokinetics_init " input; then
        sed -i -e "s/^electrokinetics_init /# electrokinetics_init /" input
    fi
fi

# Electrokinetics init rho el parameter
if [ -n "$rho_el" ]; then
    if grep -q "^electrokinetics_init_rho_el" input; then
        sed -i -e "/^electrokinetics_init_rho_el/c\electrokinetics_init_rho_el $rho_el" input
    else
        sed -i -e "/^# electrokinetics_init_rho_el/a electrokinetics_init_rho_el $rho_el" input
    fi
fi

# Temperature parameter
if [ -n "$temperature" ]; then
    if grep -q "^temperature" input; then
        sed -i -e "/^temperature/c\temperature $temperature" input
    else
        sed -i -e "/^# temperature/a temperature $temperature" input
    fi
fi

# LB fluctuations
if grep -q "^lb_fluctuations" input; then
    sed -i -e "/^lb_fluctuations/c\lb_fluctuations $fluctuations" input
else
    sed -i -e "/^# lb_fluctuations/a lb_fluctuations $fluctuations" input
fi

# LB relaxation scheme
if grep -q "^lb_relaxation_scheme" input; then
    sed -i -e "/^lb_relaxation_scheme/c\lb_relaxation_scheme $relaxation_scheme" input
else
    sed -i -e "/^# lb_relaxation_scheme/a lb_relaxation_scheme $relaxation_scheme" input
fi

# Colloid gravity
if [ -n "$gravity" ]; then
    if grep -q "^colloid_gravity" input; then
        sed -i -e "/^colloid_gravity/c\colloid_gravity $gravity" input
    else
        sed -i -e "/^# colloid_gravity/a colloid_gravity $gravity" input
    fi
fi

# Bond harmonic parameters
if [ -n "$bond_harmonic" ]; then
    IFS='_' read -r bond_on bond_k bond_r0 <<< "$bond_harmonic"

    if grep -q "^bond_harmonic_on" input; then
        sed -i -e "/^bond_harmonic_on/c\bond_harmonic_on $bond_on" input
    else
        sed -i -e "/^# bond_harmonic_on/a bond_harmonic_on $bond_on" input
    fi

    if grep -q "^bond_harmonic_k" input; then
        sed -i -e "/^bond_harmonic_k/c\bond_harmonic_k $bond_k" input
    else
        sed -i -e "/^# bond_harmonic_k/a bond_harmonic_k $bond_k" input
    fi

    if grep -q "^bond_harmonic_r0" input; then
        sed -i -e "/^bond_harmonic_r0/c\bond_harmonic_r0 $bond_r0" input
    else
        sed -i -e "/^# bond_harmonic_r0/a bond_harmonic_r0 $bond_r0" input
    fi
fi

# Angle harmonic parameters
if [ -n "$angle_harmonic" ]; then
    IFS='_' read -r angle_on angle_k angle_theta0 <<< "$angle_harmonic"

    if grep -q "^angle_harmonic_on" input; then
        sed -i -e "/^angle_harmonic_on/c\angle_harmonic_on $angle_on" input
    else
        sed -i -e "/^# angle_harmonic_on/a angle_harmonic_on $angle_on" input
    fi

    if grep -q "^angle_harmonic_k" input; then
        sed -i -e "/^angle_harmonic_k/c\angle_harmonic_k $angle_k" input
    else
        sed -i -e "/^# angle_harmonic_k/a angle_harmonic_k $angle_k" input
    fi

    if grep -q "^angle_harmonic_theta0" input; then
        sed -i -e "/^angle_harmonic_theta0/c\angle_harmonic_theta0 $angle_theta0" input
    else
        sed -i -e "/^# angle_harmonic_theta0/a angle_harmonic_theta0 $angle_theta0" input
    fi
fi

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
export OMP_NUM_THREADS=$cores

# Set GPU options for OpenCL (Ludwig uses OpenCL, not CUDA directly)
export GPU_DEVICE_ORDINAL=0
export OPENCL_VISIBLE_DEVICES=0
# Allow concurrent kernel execution
export CUDA_VISIBLE_DEVICES=0

# --- 1. Execute main task in background ---
echo "Starting background simulation task..."

# Determine MPI command (mpirun or mpiexec)
if command -v mpirun &> /dev/null; then
    MPI_CMD="mpirun"
elif command -v mpiexec &> /dev/null; then
    MPI_CMD="mpiexec"
else
    echo "ERROR: No MPI command found (mpirun or mpiexec)" >&2
    exit 1
fi

# Ludwig run with MPI on background
if [ "$mpi_procs" -gt 1 ]; then
    # Use process binding for better performance and consistency
    $MPI_CMD -np $mpi_procs --bind-to core --map-by core ./Ludwig.exe >> logs/output.txt 2>&1 &
else
    # Single process without MPI
    ./Ludwig.exe >> logs/output.txt 2>&1 &
fi

PID_BG=$!
echo "Simulation PID: $PID_BG"

# --- 2. Execute periodic plotting task while main simulation runs ---
count=0
ni=$Ninicio

# sleep 10
cp $RESULTS_DIR/vel-000000001.001-001 $RESULTS_DIR/vel-000000000.001-001

while kill -0 $PID_BG 2>/dev/null; do

    echo "[INFO] Starting periodic plot at $(date)"

    # Count output files
    file_count=$(find . -maxdepth 1 -type f -name "config.cds*" | wc -l)
    file_count=$((file_count - 2))

    if  [ "$fluid_only" == "y" ]; then
        file_count=$(find . -maxdepth 1 -type f -name "vel-*" | wc -l)
    fi

    # Adjust file_count for MPI processes
    if [ "$mpi_procs" -gt 1 ]; then
        file_count=$((file_count - mpi_procs + 1))
    fi

    # Calculate the current step number from file count
    current_step=$((paso*(file_count-1)))

    # Number of new steps to process (from ni to current_step)
    count=$((current_step - ni))

    echo "  Files found:     $file_count"
    echo "  Current step:    $current_step"
    echo "  Steps to process: $count"

    if [ "$count" -lt 0 ]; then
        sleep $deltat
        continue
    fi

    ./scripts/runplot.sh \
                -n "$count"  \
                -i "$ni"   \
                -p "$paso"   \
                -s "$single" \
                -q "$fluid_only" >> logs/outputplot.txt 2>&1

    ni=$((ni+count+paso))
    sleep $deltat

done

echo "[INFO] Periodic plotting task completed at $(date)"

# Wait to do final plot
sleep 10

# Final plot
echo "=========================================="
echo "Starting final plot..."
echo "=========================================="

# Count final files
file_count=$(find . -maxdepth 1 -type f -name "config.cds*" | wc -l)

if  [ "$fluid_only" == "y" ]; then
    file_count=$(find . -maxdepth 1 -type f -name "vel-*" | wc -l)
fi

# Adjust file_count for MPI processes
if [ "$mpi_procs" -gt 1 ]; then
    file_count=$((file_count - mpi_procs + 1))
fi

# Calculate the final step number from file count
current_step=$((paso*(file_count-1)))

# Number of steps to process for final plot
count=$((current_step - ni))

echo "Final plot - Files found:        $file_count"
echo "Final plot - Current step:       $current_step"
echo "Final plot - Steps to process:   $count"

# Execute final plot
./scripts/runplot.sh   \
                -n "$count"  \
                -i "$ni"   \
                -p "$paso"   \
                -s "$single" \
                -q "$fluid_only" >> logs/outputplot.txt 2>&1

echo "=========================================="
echo "Simulation completed successfully!"
echo "MPI processes: $mpi_procs"
echo "Output directory: $RESULTS_DIR"
echo "=========================================="

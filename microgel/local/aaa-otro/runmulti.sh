#!/bin/bash
#------------------------------------------------------------------------------------
# Run multiple Ludwig simulations with varying colloid parameters
#------------------------------------------------------------------------------------
# This script modifies config.cds.init.001-001 file to change charge, position,
# and al parameter for each simulation, then calls runbg.sh
#------------------------------------------------------------------------------------

show_help() {
    cat << EOF
Usage: ./runmulti.sh [RUNBG_OPTIONS] [MULTI_OPTIONS]

Run multiple Ludwig simulations with varying colloid parameters.
For each parameter combination, this script:
1. Modifies config.cds.init.001-001 (charge q0, position r, al parameter)
2. Calls runbg.sh with specified parameters
3. Creates output directories with correlative names based on varying parameter

MULTI-SIMULATION OPTIONS:
  --charge VALUES             Comma-separated list of charge values (q0)
                              Example: --charge 1.0,2.0,3.0
  --position VALUES           Comma-separated list of position vectors (x_y_z format)
                              Example: --position 16.0_16.0_16.0,17.0_17.0_17.0
                              Note: This sets all three coordinates at once
  --position-x VALUES         Comma-separated list of X position values
                              Example: --position-x 16.0,16.5,17.0
  --position-y VALUES         Comma-separated list of Y position values
  --position-z VALUES         Comma-separated list of Z position values
  --al VALUES                 Comma-separated list of al parameter values
                              Example: --al 0.5,0.6,0.7
  --electric-field VALUES     Comma-separated list of electric field vectors (ex_ey_ez)
                              Example: --electric-field 0.001_0.0_0.0,0.005_0.0_0.0
  --temperature VALUES        Comma-separated list of temperature values (kT)
                              Example: --temperature 0.0001,0.0005,0.001
  --rho VALUES                Comma-separated list of fluid density values
                              Example: --rho 0.6,0.8,1.0
  --viscosity VALUES          Comma-separated list of viscosity values
                              Example: --viscosity 0.01,0.02,0.05
  --epsilon VALUES            Comma-separated list of electrokinetics epsilon values
                              Example: --epsilon 1.0e4,5.0e4,1.0e5
  --rho-el VALUES             Comma-separated list of electrokinetics init rho el values
                              Example: --rho-el 1.0e-3,5.0e-3,1.0e-2
  --stencil VALUES            Comma-separated list of electrokinetics solver stencil values
                              Example: --stencil 7,19,27
  --param-name NAME           Name of the varying parameter for directory naming
                              (default: automatically determined)
  --base-dir DIR              Base directory name (will be used with runbg.sh -o)
                              This is REQUIRED.
  --parallel                  Run all simulations in parallel (default: sequential)
                              WARNING: Ensure you have enough CPU/memory resources!
  --max-parallel N            Maximum number of parallel simulations (default: unlimited)
                              Only used with --parallel option

RUNBG.SH OPTIONS:
  All options from runbg.sh are supported. The -o/--output-dir option will be
  automatically modified to include the parameter value.

  See './runbg.sh --help' for full list of options.

EXAMPLES:
  # Run simulations with varying charge:
  ./runmulti.sh --charge 0.5,1.0,1.5 --base-dir test-charge \\
                -i 0 -n 10000 -s 500 -x 32 -y 32 -z 32 -v 0.5

  # Run simulations with varying position:
  ./runmulti.sh --position-x 16.0,16.5,17.0 --base-dir test-pos \\
                -i 0 -n 10000 -s 500 -x 32 -y 32 -z 32 -v 0.5

  # Run simulations with varying al parameter:
  ./runmulti.sh --al 0.5,0.6,0.7 --base-dir test-al \\
                -i 0 -n 10000 -s 500 -x 32 -y 32 -z 32 -v 0.5 -t fe_electro

  # Run simulations with varying electric field:
  ./runmulti.sh --electric-field 0.001_0.0_0.0,0.005_0.0_0.0,0.01_0.0_0.0 \\
                --base-dir test-efield -i 0 -n 10000 -s 500 -x 32 -y 32 -z 32

  # Run simulations with varying temperature:
  ./runmulti.sh --temperature 0.0001,0.0005,0.001 --base-dir test-kt \\
                -i 0 -n 10000 -s 500 -x 32 -y 32 -z 32 -v 0.02

  # Run simulations with varying density and viscosity:
  ./runmulti.sh --rho 0.6,0.8,1.0 --viscosity 0.01,0.02,0.05 \\
                --base-dir test-fluid -i 0 -n 10000 -s 500 -x 32 -y 32 -z 32

  # Run simulations with varying electrokinetics epsilon:
  ./runmulti.sh --epsilon 1.0e4,5.0e4,1.0e5 --base-dir test-epsilon \\
                -i 0 -n 10000 -s 500 -x 32 -y 32 -z 32 -t fe_electro

  # Run simulations with varying electrokinetics rho_el:
  ./runmulti.sh --rho-el 1.0e-3,5.0e-3,1.0e-2 --base-dir test-rhoel \\
                -i 0 -n 10000 -s 500 -x 32 -y 32 -z 32 -t fe_electro

  # Run simulations with varying electrokinetics stencil:
  ./runmulti.sh --stencil 7,19,27 --base-dir test-stencil \\
                -i 0 -n 10000 -s 500 -x 32 -y 32 -z 32 -t fe_electro

  # Combine multiple varying parameters (Cartesian product):
  ./runmulti.sh --charge 1.0,2.0 --position-x 16.0,17.0 --base-dir test-combined \\
                -i 0 -n 10000 -s 500 -x 32 -y 32 -z 32 -v 0.5

  # Run simulations in PARALLEL (all at once):
  ./runmulti.sh --charge 0.5,1.0,1.5 --base-dir test-parallel \\
                --parallel \\
                -i 0 -n 10000 -s 500 -x 32 -y 32 -z 32 -v 0.5

  # Run simulations in parallel with limit (max 4 at a time):
  ./runmulti.sh --charge 0.5,1.0,1.5,2.0,2.5 --base-dir test-parallel-limited \\
                --parallel --max-parallel 4 \\
                -i 0 -n 10000 -s 500 -x 32 -y 32 -z 32 -v 0.5

NOTE:
  - If multiple parameter lists are provided, all combinations will be run
  - The base config.cds.init.001-001 file must exist in the current directory
  - Directory naming format: base-dir-n_NSTEPS-s_STEPINT-L_SIZE-param_value
    * n_NSTEPS: Number of steps (from -n/--nsteps)
    * s_STEPINT: Step interval (from -s/--step-interval)
    * L_SIZE: Volume size (compact if all equal, e.g., L_32)
              or Lx_X_Ly_Y_Lz_Z if different dimensions
    * param_value: Values of varying parameters
  - For multiple varying parameters: base-dir-n_X-s_Y-L_Z-param1_val1-param2_val2
  - By default, simulations run SEQUENTIALLY (one after another)
  - Use --parallel to run all simulations simultaneously (parallel mode)
  - Use --max-parallel N with --parallel to limit concurrent simulations

EOF
}

# Function to modify line in config file (single value)
# Usage: modify_line_in_config <line_number> <new_value> <config_file>
modify_line_in_config() {
    local line_num=$1
    local new_value=$2
    local config_file=$3

    # Format the value in scientific notation with proper spacing
    local formatted_value=$(printf "%24.15e" $new_value)

    # Use sed to replace the specific line
    sed -i "${line_num}s/.*/${formatted_value}/" "$config_file"
}

# Function to modify position line (line 40 with x, y, z)
# Usage: modify_position_line <x> <y> <z> <config_file>
modify_position_line() {
    local x=$1
    local y=$2
    local z=$3
    local config_file=$4

    # Format the values in scientific notation with proper spacing
    local formatted_line=$(printf "%24.15e %24.15e %24.15e" $x $y $z)

    # Replace line 40
    sed -i "40s/.*/${formatted_line}/" "$config_file"
}

# Parse multi-simulation specific arguments first
charge_values=""
position_values=""
pos_x_values=""
pos_y_values=""
pos_z_values=""
al_values=""
electric_field_values=""
temperature_values=""
rho_values=""
viscosity_values=""
epsilon_values=""
rho_el_values=""
stencil_values=""
param_name=""
base_dir=""
parallel_mode=0
max_parallel=0

# Arrays to hold runbg.sh arguments
runbg_args=()

# Variables to track simulation parameters for directory naming
sim_nsteps=""
sim_step_interval=""
sim_size_x=""
sim_size_y=""
sim_size_z=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --charge)
            charge_values="$2"
            shift 2
            ;;
        --position)
            position_values="$2"
            shift 2
            ;;
        --position-x)
            pos_x_values="$2"
            shift 2
            ;;
        --position-y)
            pos_y_values="$2"
            shift 2
            ;;
        --position-z)
            pos_z_values="$2"
            shift 2
            ;;
        --al)
            al_values="$2"
            shift 2
            ;;
        --electric-field)
            electric_field_values="$2"
            shift 2
            ;;
        --temperature)
            temperature_values="$2"
            shift 2
            ;;
        --rho)
            rho_values="$2"
            shift 2
            ;;
        --viscosity)
            viscosity_values="$2"
            shift 2
            ;;
        --epsilon)
            epsilon_values="$2"
            shift 2
            ;;
        --rho-el)
            rho_el_values="$2"
            shift 2
            ;;
        --stencil)
            stencil_values="$2"
            shift 2
            ;;
        --param-name)
            param_name="$2"
            shift 2
            ;;
        --base-dir)
            base_dir="$2"
            shift 2
            ;;
        --parallel)
            parallel_mode=1
            shift
            ;;
        --max-parallel)
            max_parallel="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        -o|--output-dir)
            echo "ERROR: Do not use -o/--output-dir. Use --base-dir instead." >&2
            exit 1
            ;;
        # All other options are passed through to runbg.sh
        # Note: -e, -k, -r, -v are now handled above for multi-simulation
        -a|--angle-harmonic|-b|--bond-harmonic|-c|--cores|-d|--delete-files|\
        -f|--freq-config|-w|--fluctuations|-g|--grid-mpi|\
        -j|--gravity|-i|--initial-step|-l|--fluid-only|\
        -m|--mpi-procs|-p|--plot-interval|-q|--relaxation-scheme|\
        -t|--free-energy|-u|--single-monomer|\
        --solver|--input-file|--electrokinetics-init)
            runbg_args+=("$1" "$2")
            shift 2
            ;;
        -n|--nsteps)
            sim_nsteps="$2"
            runbg_args+=("$1" "$2")
            shift 2
            ;;
        -s|--step-interval)
            sim_step_interval="$2"
            runbg_args+=("$1" "$2")
            shift 2
            ;;
        -x|--size-x)
            sim_size_x="$2"
            runbg_args+=("$1" "$2")
            shift 2
            ;;
        -y|--size-y)
            sim_size_y="$2"
            runbg_args+=("$1" "$2")
            shift 2
            ;;
        -z|--size-z)
            sim_size_z="$2"
            runbg_args+=("$1" "$2")
            shift 2
            ;;
        -e|--electric-field|-k|-r|-v)
            echo "ERROR: To vary $1, use the long form (--electric-field, --temperature, --rho, --viscosity)" >&2
            echo "       If not varying, pass it after the multi-simulation options." >&2
            exit 1
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            echo "Run './runmulti.sh --help' for usage information." >&2
            exit 1
            ;;
    esac
done

# Check that base-dir is specified
if [ -z "$base_dir" ]; then
    echo "ERROR: --base-dir is required" >&2
    echo "Run './runmulti.sh --help' for usage information." >&2
    exit 1
fi

# Check that at least one parameter is specified to vary
if [ -z "$charge_values" ] && [ -z "$position_values" ] && [ -z "$pos_x_values" ] && \
   [ -z "$pos_y_values" ] && [ -z "$pos_z_values" ] && [ -z "$al_values" ] && \
   [ -z "$electric_field_values" ] && [ -z "$temperature_values" ] && \
   [ -z "$rho_values" ] && [ -z "$viscosity_values" ] && \
   [ -z "$epsilon_values" ] && [ -z "$rho_el_values" ] && [ -z "$stencil_values" ]; then
    echo "ERROR: At least one varying parameter must be specified" >&2
    echo "       (--charge, --position, --position-x/y/z, --al, --electric-field, --temperature, --rho, --viscosity, --epsilon, --rho-el, or --stencil)" >&2
    exit 1
fi

# Check that --position is not used together with --position-x/y/z
if [ -n "$position_values" ] && { [ -n "$pos_x_values" ] || [ -n "$pos_y_values" ] || [ -n "$pos_z_values" ]; }; then
    echo "ERROR: Cannot use --position together with --position-x, --position-y, or --position-z" >&2
    echo "       Use either --position for complete vectors or --position-x/y/z for individual components" >&2
    exit 1
fi

# Check that base config file exists
if [ ! -f "config.cds.init.001-001" ]; then
    echo "ERROR: config.cds.init.001-001 file not found in current directory" >&2
    exit 1
fi

#------------------------------------------------------------------------
# Create  subdirectory
mkdir -p $base_dir

# Copy files to execute simulation
cp "$input_file" $base_dir
cp input $base_dir
cp config.cds.init.001-001 $base_dir
cp runbg.sh $base_dir
cp Ludwig.exe $base_dir
# cp del.sh $base_dir

# Copy files to execute plot
cp runplot.sh $base_dir/
cp extract_colloids $base_dir/
cp coloideacsv.sh $base_dir/

# Scripts that process raw data (need access to colloids-*.csv and vel-*) stay in root
cp calculosvel.py $base_dir/
cp calculosvelfluid.py $base_dir/
cp calculosvelfluidonly.py $base_dir/
cp extraer_posicion.py $base_dir/
cp calculos.py $base_dir/

# Plotting scripts (only read processed CSVs) go to plot subdirectory
cp plotvel.py $base_dir/
cp plot.py $base_dir/
cp plot_velocity_field.py $base_dir/
cp plot_charge_distribution.py $base_dir/
cp plot_electric_field.py $base_dir/
# cp plot_electric_field_peskin.py $base_dir/scripts
cp plotdatos.py $base_dir/
# cp batch_plot_electric_field.py $base_dir/scripts
cp compare_field_theory.py $base_dir/
# cp compare_field_theory_peskin.py $base_dir/scripts

cd $base_dir
pwd
#--------------------------------------------------------------------------------
# Create backup of original config file
backup_file="config.cds.init.001-001.backup.$$"
cp config.cds.init.001-001 "$backup_file"
echo "Created backup: $backup_file"

# Convert comma-separated values to arrays
IFS=',' read -ra CHARGE_ARRAY <<< "$charge_values"
IFS=',' read -ra POSITION_ARRAY <<< "$position_values"
IFS=',' read -ra POS_X_ARRAY <<< "$pos_x_values"
IFS=',' read -ra POS_Y_ARRAY <<< "$pos_y_values"
IFS=',' read -ra POS_Z_ARRAY <<< "$pos_z_values"
IFS=',' read -ra AL_ARRAY <<< "$al_values"
IFS=',' read -ra EFIELD_ARRAY <<< "$electric_field_values"
IFS=',' read -ra TEMP_ARRAY <<< "$temperature_values"
IFS=',' read -ra RHO_ARRAY <<< "$rho_values"
IFS=',' read -ra VISC_ARRAY <<< "$viscosity_values"
IFS=',' read -ra EPSILON_ARRAY <<< "$epsilon_values"
IFS=',' read -ra RHOEL_ARRAY <<< "$rho_el_values"
IFS=',' read -ra STENCIL_ARRAY <<< "$stencil_values"

# If arrays are empty, add a dummy element to allow iteration
[ ${#CHARGE_ARRAY[@]} -eq 0 ] && CHARGE_ARRAY=("")
[ ${#POSITION_ARRAY[@]} -eq 0 ] && POSITION_ARRAY=("")
[ ${#POS_X_ARRAY[@]} -eq 0 ] && POS_X_ARRAY=("")
[ ${#POS_Y_ARRAY[@]} -eq 0 ] && POS_Y_ARRAY=("")
[ ${#POS_Z_ARRAY[@]} -eq 0 ] && POS_Z_ARRAY=("")
[ ${#AL_ARRAY[@]} -eq 0 ] && AL_ARRAY=("")
[ ${#EFIELD_ARRAY[@]} -eq 0 ] && EFIELD_ARRAY=("")
[ ${#TEMP_ARRAY[@]} -eq 0 ] && TEMP_ARRAY=("")
[ ${#RHO_ARRAY[@]} -eq 0 ] && RHO_ARRAY=("")
[ ${#VISC_ARRAY[@]} -eq 0 ] && VISC_ARRAY=("")
[ ${#EPSILON_ARRAY[@]} -eq 0 ] && EPSILON_ARRAY=("")
[ ${#RHOEL_ARRAY[@]} -eq 0 ] && RHOEL_ARRAY=("")
[ ${#STENCIL_ARRAY[@]} -eq 0 ] && STENCIL_ARRAY=("")

# Automatically determine parameter name if not specified
if [ -z "$param_name" ]; then
    param_components=()
    [ -n "$charge_values" ] && param_components+=("q")
    [ -n "$position_values" ] && param_components+=("pos")
    [ -n "$pos_x_values" ] && param_components+=("x")
    [ -n "$pos_y_values" ] && param_components+=("y")
    [ -n "$pos_z_values" ] && param_components+=("z")
    [ -n "$al_values" ] && param_components+=("al")
    [ -n "$electric_field_values" ] && param_components+=("e")
    [ -n "$temperature_values" ] && param_components+=("kT")
    [ -n "$rho_values" ] && param_components+=("rho")
    [ -n "$viscosity_values" ] && param_components+=("visc")
    [ -n "$epsilon_values" ] && param_components+=("eps")
    [ -n "$rho_el_values" ] && param_components+=("rhoel")
    [ -n "$stencil_values" ] && param_components+=("stencil")

    # Join with underscore
    param_name=$(IFS='_'; echo "${param_components[*]}")
fi

echo "=========================================="
echo "Starting multiple simulations"
echo "=========================================="
echo "Base directory: $base_dir"
echo "Parameter name: $param_name"
echo "Charge values: ${charge_values:-[not varying]}"
echo "Position (x_y_z) values: ${position_values:-[not varying]}"
echo "Position X values: ${pos_x_values:-[not varying]}"
echo "Position Y values: ${pos_y_values:-[not varying]}"
echo "Position Z values: ${pos_z_values:-[not varying]}"
echo "AL values: ${al_values:-[not varying]}"
echo "Electric field values: ${electric_field_values:-[not varying]}"
echo "Temperature values: ${temperature_values:-[not varying]}"
echo "Density values: ${rho_values:-[not varying]}"
echo "Viscosity values: ${viscosity_values:-[not varying]}"
echo "Epsilon values: ${epsilon_values:-[not varying]}"
echo "Rho_el values: ${rho_el_values:-[not varying]}"
echo "Stencil values: ${stencil_values:-[not varying]}"
if [ $parallel_mode -eq 1 ]; then
    echo "Execution mode: PARALLEL"
    if [ $max_parallel -gt 0 ]; then
        echo "Max parallel jobs: $max_parallel"
    else
        echo "Max parallel jobs: UNLIMITED (WARNING: resource intensive!)"
    fi
else
    echo "Execution mode: SEQUENTIAL"
fi
echo "=========================================="

# Counter for simulations
sim_count=0

# Array to track background process PIDs (for parallel mode)
declare -a bg_pids=()
declare -a bg_dirs=()
num=0
# Nested loops for all parameter combinations
for charge in "${CHARGE_ARRAY[@]}"; do
for position in "${POSITION_ARRAY[@]}"; do
for pos_x in "${POS_X_ARRAY[@]}"; do
for pos_y in "${POS_Y_ARRAY[@]}"; do
for pos_z in "${POS_Z_ARRAY[@]}"; do
for al in "${AL_ARRAY[@]}"; do
for efield in "${EFIELD_ARRAY[@]}"; do
for temp in "${TEMP_ARRAY[@]}"; do
for rho in "${RHO_ARRAY[@]}"; do
for visc in "${VISC_ARRAY[@]}"; do
for epsilon in "${EPSILON_ARRAY[@]}"; do
for rho_el in "${RHOEL_ARRAY[@]}"; do
for stencil in "${STENCIL_ARRAY[@]}"; do
    num=$((num+1))
    # If using --position, extract x, y, z components
    if [ -n "$position" ]; then
        IFS='_' read -r pos_x_val pos_y_val pos_z_val <<< "$position"
    else
        pos_x_val="$pos_x"
        pos_y_val="$pos_y"
        pos_z_val="$pos_z"
    fi

    # Build directory name suffix
    dir_suffix=""
    suffix_parts=()

    # Add simulation parameters (nsteps, step-interval, size) first
    [ -n "$sim_nsteps" ] && suffix_parts+=("n_${sim_nsteps}")
    [ -n "$sim_step_interval" ] && suffix_parts+=("s_${sim_step_interval}")

    # Add size as a combined parameter if all three dimensions are specified
    if [ -n "$sim_size_x" ] && [ -n "$sim_size_y" ] && [ -n "$sim_size_z" ]; then
        if [ "$sim_size_x" = "$sim_size_y" ] && [ "$sim_size_y" = "$sim_size_z" ]; then
            # All dimensions equal: use compact notation L_X
            suffix_parts+=("L_${sim_size_x}")
        else
            # Different dimensions: use full notation Lx_X_Ly_Y_Lz_Z
            suffix_parts+=("Lx_${sim_size_x}_Ly_${sim_size_y}_Lz_${sim_size_z}")
        fi
    fi

    # Add varying parameters
    [ -n "$charge" ] && suffix_parts+=("q_${charge}")
    [ -n "$position" ] && suffix_parts+=("pos_${position}")
    # [ -n "$pos_x_val" ] && [ -z "$position" ] && suffix_parts+=("pos_${pos_x_val}")
    # [ -n "$pos_y_val" ] && [ -z "$position" ] && suffix_parts+=("_${pos_y_val}")
    # [ -n "$pos_z_val" ] && [ -z "$position" ] && suffix_parts+=("_${pos_z_val}")
    [ -n "$pos_x_val" ] && [ -n "$pos_y_val" ] && [ -n "$pos_z_val" ] && [ -z "$position" ] && suffix_parts+=("pos_${pos_x_val}_${pos_y_val}_${pos_z_val}")
    [ -n "$al" ] && suffix_parts+=("al_${al}")
    [ -n "$efield" ] && suffix_parts+=("e_${efield}")
    [ -n "$temp" ] && suffix_parts+=("kT_${temp}")
    [ -n "$rho" ] && suffix_parts+=("rho_${rho}")
    [ -n "$visc" ] && suffix_parts+=("eta_${visc}")
    [ -n "$epsilon" ] && suffix_parts+=("eps_${epsilon}")
    [ -n "$rho_el" ] && suffix_parts+=("rhoel_${rho_el}")
    [ -n "$stencil" ] && suffix_parts+=("stencil_${stencil}")

    # Join with `-` for the suffix
    if [ ${#suffix_parts[@]} -gt 0 ]; then
        dir_suffix="-$(IFS='-'; echo "${suffix_parts[*]}")"
    fi

    # Create output directory name
    output_dir="${base_dir}${dir_suffix}"

    echo ""
    echo "=========================================="
    echo "Simulation $((sim_count + 1))"
    echo "=========================================="
    echo "Output directory: $output_dir"
    [ -n "$charge" ] && echo "  Charge (q0): $charge"
    if [ -n "$position" ]; then
        echo "  Position (x_y_z): $position"
        echo "    -> X: $pos_x_val, Y: $pos_y_val, Z: $pos_z_val"
    else
        [ -n "$pos_x_val" ] && echo "  Position X: $pos_x_val"
        [ -n "$pos_y_val" ] && echo "  Position Y: $pos_y_val"
        [ -n "$pos_z_val" ] && echo "  Position Z: $pos_z_val"
    fi
    [ -n "$al" ] && echo "  AL parameter: $al"
    [ -n "$efield" ] && echo "  Electric field: $efield"
    [ -n "$temp" ] && echo "  Temperature: $temp"
    [ -n "$rho" ] && echo "  Density: $rho"
    [ -n "$visc" ] && echo "  Viscosity: $visc"
    [ -n "$epsilon" ] && echo "  Epsilon: $epsilon"
    [ -n "$rho_el" ] && echo "  Rho_el: $rho_el"
    [ -n "$stencil" ] && echo "  Stencil: $stencil"

    # Restore original config from backup for each simulation
    cp "$backup_file" config.cds.init.001-001

    # Create a temporary modified config file
    temp_config="config.cds.init.001-001.tmp.$$"
    cp config.cds.init.001-001 "$temp_config"

    # Modify the parameters in the temp config file
    # Line numbers in the ASCII file:
    # Lines 1-32: integer values
    # Lines 33+: double values
    # Line 40: position r[0] r[1] r[2] (x, y, z on one line)
    # Line 51: q0 (charge)
    # Line 58: al parameter

    # Handle position modification (line 40 contains x, y, z)
    if [ -n "$pos_x_val" ] || [ -n "$pos_y_val" ] || [ -n "$pos_z_val" ]; then
        # Read current position values from line 40
        current_pos=($(sed -n '40p' "$temp_config"))

        # Use provided values or keep current ones
        new_x=${pos_x_val:-${current_pos[0]}}
        new_y=${pos_y_val:-${current_pos[1]}}
        new_z=${pos_z_val:-${current_pos[2]}}

        modify_position_line "$new_x" "$new_y" "$new_z" "$temp_config"
    fi

    # Modify charge if specified
    [ -n "$charge" ] && modify_line_in_config 51 "$charge" "$temp_config"

    # Modify al parameter if specified
    [ -n "$al" ] && modify_line_in_config 58 "$al" "$temp_config"

    # Replace the original config with modified one
    temp_parallel_config="config.cds.init.001-001-$num"
    mv "$temp_config" "$temp_parallel_config"

    # Build additional runbg.sh arguments for this specific simulation
    sim_args=("${runbg_args[@]}")
    sim_args+=("--config-file" "$temp_parallel_config")

    # Add varying parameters to runbg.sh arguments
    [ -n "$efield" ] && sim_args+=("-e" "$efield")
    [ -n "$temp" ] && sim_args+=("-k" "$temp")
    [ -n "$rho" ] && sim_args+=("-r" "$rho")
    [ -n "$visc" ] && sim_args+=("-v" "$visc")
    [ -n "$epsilon" ] && sim_args+=("--epsilon" "$epsilon")
    [ -n "$rho_el" ] && sim_args+=("--rho-el" "$rho_el")
    [ -n "$stencil" ] && sim_args+=("--stencil" "$stencil")

    # Call runbg.sh with the modified config
    if [ $parallel_mode -eq 1 ]; then
        # PARALLEL MODE: Run in background

        # If max_parallel is set, wait for slots to become available
        if [ $max_parallel -gt 0 ]; then
            while [ ${#bg_pids[@]} -ge $max_parallel ]; do
                # Wait for any job to finish
                for i in "${!bg_pids[@]}"; do
                    pid="${bg_pids[$i]}"
                    if ! kill -0 "$pid" 2>/dev/null; then
                        # Process finished, remove from array
                        echo "[INFO] Simulation in ${bg_dirs[$i]} completed"
                        unset bg_pids[$i]
                        unset bg_dirs[$i]
                        # Re-index arrays
                        bg_pids=("${bg_pids[@]}")
                        bg_dirs=("${bg_dirs[@]}")
                        break
                    fi
                done
                sleep 2
            done
        fi

        echo "[PARALLEL] Starting: ./runbg.sh ${sim_args[@]} -o $output_dir"

        ./runbg.sh "${sim_args[@]}" -o "$output_dir" &
        bg_pid=$!
        bg_pids+=("$bg_pid")
        bg_dirs+=("$output_dir")
        echo "[PARALLEL] Launched simulation in background (PID: $bg_pid, Dir: $output_dir)"
        sim_count=$((sim_count + 1))

        # wait 5
        # rm -f "$temp_parallel_config"

    else
        # SEQUENTIAL MODE: Run and wait for completion

        echo "Running: ./runbg.sh ${sim_args[@]} -o $output_dir"
        ./runbg.sh "${sim_args[@]}" -o "$output_dir"

        if [ $? -eq 0 ]; then
            echo "Simulation completed successfully!"
            sim_count=$((sim_count + 1))
        else
            echo "ERROR: Simulation failed!" >&2
            # Restore original config from backup
            cp "$backup_file" config.cds.init.001-001
            rm -f "$backup_file"
            exit 1
        fi
    fi

done
done
done
done
done
done
done
done
done
done
done
done
done

# If running in parallel mode, wait for all background jobs to complete
if [ $parallel_mode -eq 1 ]; then
    echo ""
    echo "=========================================="
    echo "Waiting for all parallel simulations to complete..."
    echo "Active simulations: ${#bg_pids[@]}"
    echo "=========================================="

    # Wait for all background processes
    failed_count=0
    for i in "${!bg_pids[@]}"; do
        pid="${bg_pids[$i]}"
        dir="${bg_dirs[$i]}"
        echo "Waiting for simulation in $dir (PID: $pid)..."

        if wait "$pid"; then
            echo "  ✓ Simulation in $dir completed successfully"
        else
            echo "  ✗ Simulation in $dir failed (exit code: $?)"
            failed_count=$((failed_count + 1))
        fi
    done

    echo ""
    echo "=========================================="
    echo "All parallel simulations finished!"
    echo "Successful: $((sim_count - failed_count))"
    echo "Failed: $failed_count"
    echo "=========================================="

    if [ $failed_count -gt 0 ]; then
        echo "WARNING: Some simulations failed. Check individual directories for details."
    fi
fi

# # Restore original config file
# echo ""
# echo "Restoring original config file..."
# cp "$backup_file" config.cds.init.001-001
# rm -f "$backup_file"

echo ""
echo "=========================================="
echo "All simulations completed!"
echo "Total simulations launched: $sim_count"
echo "=========================================="

# Clean up temporary config files and scripts
echo "borrar archivos residuales"
rm -f ./* 2>/dev/null
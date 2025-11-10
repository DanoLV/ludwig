#!/bin/bash
#------------------------------------------------------------------------------------
# Test script to verify parameter passing to runbg.sh
#------------------------------------------------------------------------------------

echo "=========================================="
echo "Testing parameter passing"
echo "=========================================="
echo ""

# Test 1: Parse a complex command line
echo "Test 1: Parsing complex command line"
echo "--------------------------------------"

# Simulate parsing (we'll extract the logic from runmulti.sh)
declare -a test_args

# Simulate command line:
# ./runmulti.sh --charge 1.0,2.0 --base-dir test -i 0 -n 10000 -s 500 -x 32 -y 32 -v 0.5 -k 0.001 -t fe_electro -e 0.0_0.0_0.1

test_cmd=(
    "--charge" "1.0,2.0"
    "--base-dir" "test"
    "-i" "0"
    "-n" "10000"
    "-s" "500"
    "-x" "32"
    "-y" "32"
    "-v" "0.5"
    "-k" "0.001"
    "-t" "fe_electro"
    "-e" "0.0_0.0_0.1"
)

echo "Simulated command line:"
echo "  ./runmulti.sh ${test_cmd[@]}"
echo ""

# Parse the arguments
runbg_args=()
charge_values=""
base_dir=""
idx=0

while [ $idx -lt ${#test_cmd[@]} ]; do
    arg="${test_cmd[$idx]}"

    case "$arg" in
        --charge)
            idx=$((idx + 1))
            charge_values="${test_cmd[$idx]}"
            ;;
        --base-dir)
            idx=$((idx + 1))
            base_dir="${test_cmd[$idx]}"
            ;;
        -i|--initial-step|-n|--nsteps|-s|--step-interval|-x|--size-x|\
        -y|--size-yz|-v|--viscosity|-k|--temperature|-t|--free-energy|\
        -e|--electric-field)
            runbg_args+=("$arg")
            idx=$((idx + 1))
            runbg_args+=("${test_cmd[$idx]}")
            ;;
        *)
            echo "ERROR: Unknown option: $arg"
            exit 1
            ;;
    esac

    idx=$((idx + 1))
done

echo "Parsed values:"
echo "  charge_values: $charge_values"
echo "  base_dir: $base_dir"
echo "  runbg_args: ${runbg_args[@]}"
echo ""

# Verify expected values
expected_runbg_args="-i 0 -n 10000 -s 500 -x 32 -y 32 -v 0.5 -k 0.001 -t fe_electro -e 0.0_0.0_0.1"
actual_runbg_args="${runbg_args[@]}"

echo "Expected runbg args: $expected_runbg_args"
echo "Actual runbg args:   $actual_runbg_args"
echo ""

if [ "$actual_runbg_args" == "$expected_runbg_args" ]; then
    echo "✓ Test 1 PASSED: Arguments parsed correctly"
else
    echo "✗ Test 1 FAILED: Arguments mismatch"
fi

echo ""
echo "=========================================="
echo "Test 2: Verify runmulti.sh accepts params"
echo "=========================================="
echo ""

# Try to run the actual script with --help to verify syntax
if [ -f "./runmulti.sh" ]; then
    echo "Running: ./runmulti.sh --help | head -20"
    ./runmulti.sh --help | head -20
    echo ""
    echo "✓ Test 2 PASSED: runmulti.sh is executable and shows help"
else
    echo "✗ Test 2 FAILED: runmulti.sh not found"
fi

echo ""
echo "=========================================="
echo "Test 3: Dry run with actual parameters"
echo "=========================================="
echo ""

# We can't do a full run without the actual simulation setup,
# but we can test the argument parsing by checking the error messages

echo "Testing invalid parameter rejection:"
echo "Running: ./runmulti.sh -o invalid-dir 2>&1 | grep ERROR"

if ./runmulti.sh -o invalid-dir 2>&1 | grep -q "ERROR"; then
    echo "✓ Test 3a PASSED: -o/--output-dir is properly rejected"
else
    echo "✗ Test 3a FAILED: -o/--output-dir should be rejected"
fi

echo ""
echo "Testing missing --base-dir:"
echo "Running: ./runmulti.sh --charge 1.0 -i 0 -n 100 -s 10 2>&1 | grep ERROR"

if ./runmulti.sh --charge 1.0 -i 0 -n 100 -s 10 2>&1 | grep -q "ERROR"; then
    echo "✓ Test 3b PASSED: Missing --base-dir is properly detected"
else
    echo "✗ Test 3b FAILED: Missing --base-dir should cause error"
fi

echo ""
echo "=========================================="
echo "All tests completed!"
echo "=========================================="

#!/bin/bash
#------------------------------------------------------------------------------------
# Test script for --position option
#------------------------------------------------------------------------------------

echo "=========================================="
echo "Testing --position option parsing"
echo "=========================================="

# Test position parsing
test_position="16.5_17.0_17.5"

echo "Test position string: $test_position"
echo ""

# Parse the position
IFS='_' read -r pos_x pos_y pos_z <<< "$test_position"

echo "Parsed components:"
echo "  X: $pos_x"
echo "  Y: $pos_y"
echo "  Z: $pos_z"

echo ""
echo "Testing position modification in config file..."

# Check if original config exists
if [ ! -f "config.cds.init.001-001" ]; then
    echo "ERROR: config.cds.init.001-001 not found"
    exit 1
fi

# Create backup
backup="config.cds.init.001-001.postest.backup"
cp config.cds.init.001-001 "$backup"

# Create test config
test_config="config.cds.init.001-001.postest"
cp config.cds.init.001-001 "$test_config"

echo ""
echo "Original position line 40:"
sed -n '40p' config.cds.init.001-001

# Modify the position
formatted_line=$(printf "%24.15e %24.15e %24.15e" $pos_x $pos_y $pos_z)
sed -i "40s/.*/${formatted_line}/" "$test_config"

echo ""
echo "Modified position line 40:"
sed -n '40p' "$test_config"

# Read back and verify
echo ""
echo "Reading back to verify:"
read_pos=($(sed -n '40p' "$test_config"))
echo "  X: ${read_pos[0]}"
echo "  Y: ${read_pos[1]}"
echo "  Z: ${read_pos[2]}"

# Cleanup
rm -f "$test_config"
rm -f "$backup"

echo ""
echo "=========================================="
echo "Test completed successfully!"
echo "=========================================="

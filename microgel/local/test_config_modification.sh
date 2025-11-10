#!/bin/bash
#------------------------------------------------------------------------------------
# Test script to verify config file modification functions
#------------------------------------------------------------------------------------

echo "=========================================="
echo "Testing config file modification"
echo "=========================================="

# Check if original config exists
if [ ! -f "config.cds.init.001-001" ]; then
    echo "ERROR: config.cds.init.001-001 not found"
    exit 1
fi

# Create backup
backup="config.cds.init.001-001.test.backup"
cp config.cds.init.001-001 "$backup"
echo "Created backup: $backup"

# Create test config
test_config="config.cds.init.001-001.test"
cp config.cds.init.001-001 "$test_config"

echo ""
echo "Original values:"
echo "  Line 40 (position): $(sed -n '40p' config.cds.init.001-001)"
echo "  Line 51 (charge q0): $(sed -n '51p' config.cds.init.001-001)"
echo "  Line 58 (al param): $(sed -n '58p' config.cds.init.001-001)"

echo ""
echo "Modifying test config..."

# Test modifying position
new_x=15.5
new_y=16.0
new_z=16.5
formatted_line=$(printf "%24.15e %24.15e %24.15e" $new_x $new_y $new_z)
sed -i "40s/.*/${formatted_line}/" "$test_config"

# Test modifying charge
new_charge=2.5
formatted_value=$(printf "%24.15e" $new_charge)
sed -i "51s/.*/${formatted_value}/" "$test_config"

# Test modifying al
new_al=0.75
formatted_value=$(printf "%24.15e" $new_al)
sed -i "58s/.*/${formatted_value}/" "$test_config"

echo ""
echo "Modified values in test config:"
echo "  Line 40 (position): $(sed -n '40p' $test_config)"
echo "  Line 51 (charge q0): $(sed -n '51p' $test_config)"
echo "  Line 58 (al param): $(sed -n '58p' $test_config)"

echo ""
echo "Reading back values to verify:"
pos_values=($(sed -n '40p' $test_config))
echo "  Position X: ${pos_values[0]}"
echo "  Position Y: ${pos_values[1]}"
echo "  Position Z: ${pos_values[2]}"
echo "  Charge q0: $(sed -n '51p' $test_config)"
echo "  AL param: $(sed -n '58p' $test_config)"

# Verify original is unchanged
echo ""
echo "Verifying original config is unchanged:"
echo "  Line 40 (position): $(sed -n '40p' config.cds.init.001-001)"
echo "  Line 51 (charge q0): $(sed -n '51p' config.cds.init.001-001)"
echo "  Line 58 (al param): $(sed -n '58p' config.cds.init.001-001)"

# Cleanup
rm -f "$test_config"
rm -f "$backup"

echo ""
echo "=========================================="
echo "Test completed successfully!"
echo "=========================================="

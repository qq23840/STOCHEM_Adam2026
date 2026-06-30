#!/bin/bash
# Script to replace local file system paths with generic relative paths
# for publication

echo "Cleaning up local paths in workspace..."

# Replace paths in .for files
echo "Processing .for files..."
for file in sensitivity_*.for; do
    if [ -f "$file" ]; then
        sed -i "s|/user/work/qq23840/stochem_met/|./met_data/|g" "$file"
        sed -i "s|/user/home/qq23840/stochem/model_files/|./model_files/|g" "$file"
        sed -i "s|/user/work/qq23840/stochem_output/|./output/|g" "$file"
        echo "  ✓ $file"
    fi
done

# Replace paths in .SOA2 files
echo "Processing .SOA2 files..."
for file in model_files/*.SOA2; do
    if [ -f "$file" ]; then
        sed -i "s|/user/work/qq23840/stochem_met/|./met_data/|g" "$file"
        sed -i "s|/user/home/qq23840/stochem/model_files/|./model_files/|g" "$file"
        sed -i "s|/user/work/qq23840/stochem_output/|./output/|g" "$file"
        echo "  ✓ $file"
    fi
done

echo ""
echo "Cleanup complete!"
echo ""
echo "Summary of replacements:"
echo "  /user/work/qq23840/stochem_met/       → ./met_data/"
echo "  /user/home/qq23840/stochem/model_files/ → ./model_files/"
echo "  /user/work/qq23840/stochem_output/    → ./output/"
echo ""
echo "Verifying no local paths remain..."
if grep -r "/user/home/qq23840\|/user/work/qq23840" *.for model_files/*.SOA2 2>/dev/null; then
    echo "WARNING: Some local paths still found!"
else
    echo "✓ All local paths successfully replaced!"
fi

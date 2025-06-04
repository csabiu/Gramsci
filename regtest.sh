#!/bin/bash
set -e

# Build executables
make -C src all >/dev/null

# Run gramsci to compute 2pcf
bin/gramsci -gal example/test.gal -ran example/test.ran -rmin 1.0 -rmax 30.0 -nbins 10 -wgt -out example/test.2pcf -2pcf >/dev/null

# Capture summary statistics
TEST_2PCF_LINES=$(wc -l < example/test.2pcf)
TEST_2PCF_MD5=$(md5sum example/test.2pcf | awk '{print $1}')

# Prepare input for domain_decomposition
awk '{print $0,1,0,1}' example/test.gal > allpoints
awk '{print $0,0,0,1}' example/test.ran >> allpoints

# Run domain_decomposition with rmax 30
bin/domain_decomposition 30.0 >/dev/null

# Capture summary statistics for domain_decomposition output
LOADNODES_LINES=$(wc -l < 1.loadnodes)
LOADNODES_MD5=$(md5sum 1.loadnodes | awk '{print $1}')

# Cleanup large files
rm 1.loadnodes allpoints

# Load expected values
source regtest.expected

# Compare results
status=0
if [ "$TEST_2PCF_LINES" != "$EXP_TEST_2PCF_LINES" ]; then
  echo "Mismatch in 2pcf line count" >&2
  status=1
fi
if [ "$TEST_2PCF_MD5" != "$EXP_TEST_2PCF_MD5" ]; then
  echo "Mismatch in 2pcf md5" >&2
  status=1
fi
if [ "$LOADNODES_LINES" != "$EXP_LOADNODES_LINES" ]; then
  echo "Mismatch in loadnodes line count" >&2
  status=1
fi
if [ "$LOADNODES_MD5" != "$EXP_LOADNODES_MD5" ]; then
  echo "Mismatch in loadnodes md5" >&2
  status=1
fi

# Output summary
echo "$TEST_2PCF_LINES"
echo "$TEST_2PCF_MD5  example/test.2pcf"
echo "$LOADNODES_LINES"
echo "$LOADNODES_MD5  1.loadnodes"

exit $status

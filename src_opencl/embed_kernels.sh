#!/bin/sh
# embed_kernels.sh kernels.cl cl_kernels_module.F90
#
# Generate a Fortran module that returns the OpenCL kernel source as a string,
# so the built binary has no runtime dependency on the .cl file's location.
# One Fortran statement per kernel line (not continuation lines) keeps us clear
# of the standard 255-continuation-line limit for large kernels.
SRC="$1"
OUT="$2"

{
  echo "! AUTO-GENERATED from kernels.cl by embed_kernels.sh -- DO NOT EDIT."
  echo "module cl_kernels_module"
  echo "  implicit none"
  echo "  character(len=1), parameter :: NL = char(10)"
  echo "contains"
  echo "  function kernel_src() result(s)"
  echo "    character(len=:), allocatable :: s"
  echo "    s = ''"
  # Escape single quotes (double them) so each line is a valid Fortran literal.
  awk '{
    line = $0
    gsub(/\x27/, "\x27\x27", line)
    printf "    s = s // \x27%s\x27 // NL\n", line
  }' "$SRC"
  echo "  end function kernel_src"
  echo "end module cl_kernels_module"
} > "$OUT"

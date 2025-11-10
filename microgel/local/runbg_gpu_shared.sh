#!/bin/bash
# Version of runbg.sh optimized for sharing GPU between processes
# This limits GPU memory usage per process to allow 2 processes on 1 GPU

# Source the original runbg.sh but with GPU memory limit
export CUDA_VISIBLE_DEVICES=0
export GPU_MAX_HEAP_SIZE=50  # Limit to 50% of GPU memory
export GPU_FORCE_64BIT_PTR=0

# Execute original runbg.sh with all arguments passed through
exec "$(dirname "$0")/runbg.sh" "$@"

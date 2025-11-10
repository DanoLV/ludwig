#!/bin/bash
# Enable NVIDIA MPS for GPU sharing between multiple Ludwig processes

echo "=== Enabling NVIDIA Multi-Process Service (MPS) ==="
echo "This allows multiple processes to share the GPU efficiently"
echo ""

# Check if running as root
if [ "$EUID" -ne 0 ]; then
    echo "Please run as root or with sudo"
    exit 1
fi

# Create directories for MPS
mkdir -p /tmp/nvidia-mps
mkdir -p /tmp/nvidia-log

# Set permissions
chmod 777 /tmp/nvidia-mps
chmod 777 /tmp/nvidia-log

# Export environment variables
export CUDA_VISIBLE_DEVICES=0
export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps
export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-log

# Stop any existing MPS daemon
echo "Stopping any existing MPS daemon..."
nvidia-cuda-mps-control -d 2>/dev/null

# Start MPS daemon
echo "Starting MPS daemon..."
nvidia-cuda-mps-control -d

# Check status
sleep 1
if pgrep -x "nvidia-cuda-mps" > /dev/null; then
    echo "✓ MPS daemon started successfully"
    echo ""
    echo "Now you can run multiple Ludwig simulations in parallel"
    echo "They will share the GPU efficiently"
else
    echo "✗ Failed to start MPS daemon"
    exit 1
fi

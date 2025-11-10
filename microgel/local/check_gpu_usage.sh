#!/bin/bash
# Monitor GPU usage while Ludwig runs

echo "=== Monitoring GPU usage ==="
echo "Press Ctrl+C to stop"
echo ""

# Check if Ludwig is running
if ! pgrep -f Ludwig.exe > /dev/null; then
    echo "WARNING: No Ludwig.exe process found running"
    echo "Start a simulation first, then run this script"
    echo ""
fi

# Monitor GPU every 2 seconds
while true; do
    clear
    echo "=== GPU Status at $(date) ==="
    echo ""

    nvidia-smi --query-gpu=index,name,temperature.gpu,utilization.gpu,utilization.memory,memory.used,memory.total --format=csv,noheader,nounits

    echo ""
    echo "=== Ludwig processes ==="
    ps aux | grep Ludwig.exe | grep -v grep | awk '{printf "PID: %s  CPU: %s%%  MEM: %s%%  CMD: %s\n", $2, $3, $4, $11}'

    echo ""
    echo "=== GPU Processes ==="
    nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader

    sleep 2
done

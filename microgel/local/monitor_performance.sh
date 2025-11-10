#!/bin/bash
# Monitor performance of Ludwig simulations

echo "=== System Resources ==="
echo "Date: $(date)"
echo ""

echo "=== CPU Usage per Ludwig process ==="
ps aux | grep Ludwig.exe | grep -v grep | awk '{print "PID: "$2" CPU: "$3"% MEM: "$4"%"}'
echo ""

echo "=== Total CPU cores usage ==="
mpstat -P ALL 1 1 | tail -n +4
echo ""

echo "=== Memory usage ==="
free -h
echo ""

echo "=== I/O stats ==="
iostat -x 1 2 | tail -n +4
echo ""

echo "=== WSL Memory from Windows perspective ==="
cat /proc/meminfo | grep -E "MemTotal|MemFree|MemAvailable|Cached|SwapTotal|SwapFree"

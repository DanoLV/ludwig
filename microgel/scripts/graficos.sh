#!/bin/bash
#################################################################################
# Script to collect all velocidades.png files into a single directory
#################################################################################
mkdir -p ../graficos
find . -type f -name velocidades.png | while read filepath; do
  # quita el ./ inicial y extrae la primera carpeta del path relativo
  firstdir=$(echo "$filepath" | cut -d'/' -f2)
  cp "$filepath" "../graficos/${firstdir}.png"
done
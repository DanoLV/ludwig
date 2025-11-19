#!/bin/bash
#------------------------------------------------------------------------------------
# Delete unnecesary files.
# Leave last file so that restart is possible.
#------------------------------------------------------------------------------------

find . -type d -not -path '.' | while read dir; do
  ls -1 "$dir"/dist-0* 2>/dev/null | sort | head -n -1 | xargs -r rm --
  ls -1 "$dir"/rho-0* 2>/dev/null | sort | head -n -1 | xargs -r rm --
done
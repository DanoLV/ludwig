#!/bin/bash
#------------------------------------------------------------------------------------
# Delete unnecesary files.
# Leave last file so that restart is possible.
#------------------------------------------------------------------------------------

# find . -type d -not -path '.' | while read -r dir; do
#   ls -1 "$dir"/dist-0* 2>/dev/null | sort | head -n -1 | xargs -r rm --
#   ls -1 "$dir"/rho-0* 2>/dev/null | sort | head -n -1 | xargs -r rm --
# done

find . -type d | while IFS= read -r dir; do
    # Borrar dist-0* manteniendo el último
    files=( "$dir"/dist-0* )
    if (( ${#files[@]} > 1 )); then
        printf "%s\0" "${files[@]}" \
        | sort -z \
        | head -z -n -1 \
        | xargs -0 rm -v
    fi

    # Borrar rho-0* manteniendo el último
    files=( "$dir"/rho-0* )
    if (( ${#files[@]} > 1 )); then
        printf "%s\0" "${files[@]}" \
        | sort -z \
        | head -z -n -1 \
        | xargs -0 rm -v
    fi
done

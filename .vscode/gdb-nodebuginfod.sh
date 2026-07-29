#!/usr/bin/env sh
# Wrapper de GDB para VSCode/cppdbg en WSL2.
# Anula DEBUGINFOD_URLS antes de lanzar gdb: la descarga automatica de debug
# info (debuginfod.ubuntu.com) cuelga el arranque del debugger en WSL2.
unset DEBUGINFOD_URLS
exec /usr/bin/gdb "$@"

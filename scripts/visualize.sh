#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
OUT="$SCRIPT_DIR/../tests/graph.svg"
"$1" --task "$2" --graphviz > /dev/null
neato -Tsvg -n graph.dot > "$OUT"
rm graph.dot
xdg-open "$OUT"

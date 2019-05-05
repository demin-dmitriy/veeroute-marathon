#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
EVAL_DIR="$SCRIPT_DIR/../eval"
OUT="$EVAL_DIR/graph.svg"
"$1" --task "$2" --graphviz > /dev/null
neato -Tsvg -n graph.dot > "$OUT"
mv graph.dot "$EVAL_DIR"
xdg-open "$OUT"

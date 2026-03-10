#!/bin/bash
#
# Run GROMACS + CHAP pipeline inside a Podman container.
#
# Usage:
#   ./run_chap.sh -i /path/to/input -o /path/to/output [-- extra python args]
#
# The input directory is mounted read-only as /input inside the container.
# The output directory is mounted as /output inside the container.

IMAGE="localhost/chap:latest"
NAME="chap"

# --- Defaults ---
INPUT=""
OUTPUT="$(pwd)/run"
EXTRA_ARGS=""
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Parse arguments ---
usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> [-- extra args]"
    echo ""
    echo "  -i  Input directory with CHARMM-GUI/GROMACS files (required), mounted as /input"
    echo "  -o  Output directory (default: ./run)"
    echo "  --  Everything after this is passed directly to gromacs_run.py"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) INPUT="$2"; shift 2 ;;
        -o) OUTPUT="$2"; shift 2 ;;
        --) shift; EXTRA_ARGS="$*"; break ;;
        -h|--help) usage ;;
        *)  echo "Unknown option: $1"; usage ;;
    esac
done

if [[ -z "$INPUT" ]]; then
    echo "Error: -i <input_dir> is required."
    usage
fi

INPUT="$(realpath "$INPUT")"

mkdir -p "$OUTPUT"
OUTPUT="$(realpath "$OUTPUT")"

if [[ ! -d "$INPUT" ]]; then
    echo "Error: input directory '$INPUT' does not exist."
    exit 1
fi

# --- Remove any leftover container with the same name ---
podman rm -f "$NAME" 2>/dev/null

# --- Run the pipeline inside the container ---
podman run --rm \
    --name "$NAME" \
    -e HOST_UID="$(id -u)" \
    -e HOST_GID="$(id -g)" \
    -v "${INPUT}:/input:ro,Z" \
    -v "${OUTPUT}:/output:Z" \
    -v "${SCRIPT_DIR}/gromacs_run.py:/gromacs_run.py:ro,Z" \
    --device nvidia.com/gpu=all \
    "$IMAGE" \
    python3 /gromacs_run.py

#!/bin/bash

IMAGE="localhost/chap:latest"
NAME="chap"

podman rm -f "$NAME" 2>/dev/null

podman run -it \
  --name "$NAME" \
  -v "$HOME/projects/CHAP/:/workspace:Z" \
  "$IMAGE" /bin/bash
  #--device nvidia.com/gpu=all \
  #--privileged \
  #-v /etc/OpenCL/vendors:/etc/OpenCL/vendors:ro \

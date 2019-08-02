#!/bin/bash

mkdir -p tools/mspc
wget -P tools/mspc/ https://dot.net/v1/dotnet-install.sh
bash tools/mspc/dotnet-install.sh -i tools/mspc/
wget -P tools/mspc https://github.com/Genometric/MSPC/releases/download/v3.3.1/mspc_v3.3.1.zip
unzip -d tools/mspc -o tools/mspc/mspc_v3.3.1.zip


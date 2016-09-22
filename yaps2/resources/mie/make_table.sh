#!/bin/bash

set -euo pipefail

VARIANTS=$1
MENDEL=$2

echo -e "Family\tVariants\tMIE"
join <(sort $VARIANTS) <(cut -f1,5 $MENDEL | sort )

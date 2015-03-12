#!/bin/sh

set -x

[ -z "$1" ] && echo "No argument supplied" && exit 1

mkdir -p "$1"/lib "$1"/include "$1"/scripts
cp lib/libocca.so "$1"/lib
if [ `uname` = "Darwin" ]; then
  # install_name_tool -id '@rpath/libocca.so' "$1"/lib/libocca.so
  install_name_tool -id "$1"/lib/libocca.so "$1"/lib/libocca.so
fi
cp -r include/* "$1"/include
cp  scripts/shellTools.sh "$1"/scripts

#!/bin/bash

if [ -z "$PYTHONPATH" ] ; then
  export PYTHONPATH=$PWD
else
  export PYTHONPATH=$PWD:$PYTHONPATH
fi

mkdir -p gallery
cd gallery
for fic in ../examples/*.py ; do
  echo run $fic
  ${PYTHON:-python} $fic
done

#!/bin/bash -ue

which niet >& null || pip install niet

DEBUG=false

function echoerr(){ echo $@ 1>&2; }

[ $# -eq 0 ] && IN=$(cat) || IN=$@

$DEBUG && echoerr "IN=$IN"

for i in $@ $IN
do
  $DEBUG && echoerr "i=$i"
  [[ "${i/\/}" == "$i" ]] || i="${i%/*}"
  [[ "${i/.yaml}" == "$i" ]] && i="$i.yaml"
  if [ ! -e $i ]
  then
    echoerr "ERROR: Cannot find metadata $i"
    exit 3
  fi
  for j in submetadata sources
  do
    $DEBUG && echoerr "niet $j $i"
    OUT=$(niet $j $i || echo "ERROR") 
    $DEBUG && echoerr "OUT=$OUT"
    [[ "${OUT/ERROR}" == "$OUT" ]] && $BASH_SOURCE $OUT
    echo $i
  done
done | sort -u
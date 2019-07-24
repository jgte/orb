#!/bin/bash -ue

TMP=/tmp/$BASH_SOURCE.$RANDOM

function grep-it()
{
  grep -v \# $1 | awk -F: '/:/ {print $0}' | sort
}

for i in $@
do
  if [ ! -e $TMP.1 ]
  then
    grep-it $i > $TMP.1
  else
    grep-it $i > $TMP.2
    diff -y $TMP.1 $TMP.2 | grep -v -e '[\>|\||\<]' | awk '{print $1,$2}' > $TMP.d
    mv -f $TMP.d $TMP.1
  fi
done
cat $TMP.1
rm -f $TMP*
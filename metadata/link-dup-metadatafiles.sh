#!/bin/bash -u

DIR=$(cd $(dirname $BASH_SOURCE);pwd)
SOURCE=$DIR/$1
SINK=$DIR/$2

if [ ! -d "$SOURCE" ]
then
  echo "ERROR: cannot find dir '$SOURCE'"
  exit 3
fi
if [ ! -d "$SINK" ]
then
  echo "ERROR: cannot find dir '$SINK'"
  exit 3
fi

cd $SINK

for i in $SINK/*.yaml
do
  j="$SOURCE/$(basename $i)"
  if [ -e "$j" ] && diff -q "$i" "$j" > /dev/null
  then
    ln -sfv "../$(basename $SOURCE)/$(basename $i)" .
  fi
done

cd - > /dev/null
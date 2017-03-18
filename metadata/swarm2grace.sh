#!/bin/bash

for i in *swarm*decomp*
do
  [ ! -e ${i/swarm.combined/grace.gfz} ] && cp $i ${i/swarm.combined/grace.gfz}
done

for i in *grace*decomp*
do
  sed 's/gswarm/GSWARM/g' $i  | \
  sed 's/full/rl05/g' | \
  sed 's/combined/gfz/g' | \
  sed 's/swarm/grace/g' | \
  sed 's/GSWARM/gswarm/g' > $i.tmp
  if [ -z "$(diff $i $i.tmp)" ]
  then
    rm -f $i.tmp
    continue
  fi
  echo ================
  echo $i
  echo ================
  cat $i
  echo
  echo ---
  cat $i.tmp
  echo
  echo "Proceed? [Y/n]"
  read ANSWER
  if [ "$ANSWER" == "n" ]
  then
    rm -fv $i.tmp
    exit
  fi
  mv -vf $i.tmp $i
done
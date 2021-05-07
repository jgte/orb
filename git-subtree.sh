#!/bin/bash -ue

if [ $# -lt 1 ]
then
  echo "ERROR: Need one input argument: git operation (add/push/pull)"
  exit 3
fi

[[ "${@/echo}" == "$@" ]] && ECHO= || ECHO=echo

DIR=$(cd $(dirname $BASH_SOURCE);pwd)

while IFS= read -r line || [[ -n "$line" ]]; do
  # https://www.yeahshecodes.com/git/a-simple-git-subtree-tutorial
  $ECHO git subtree $1  $line || true
done < $DIR/${BASH_SOURCE%.sh}.list

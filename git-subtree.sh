#!/bin/bash -ue

if [ $# -lt 1 ]
then
  echo "ERROR: Need one input argument: git operation (add/push/pull)"
  exit 3
fi

[[ "${@/echo}" == "$@" ]] && ECHO= || ECHO=echo

DIR=$(cd $(dirname $BASH_SOURCE);pwd)

while IFS= read -r line || [[ -n "$line" ]]; do
  SUBDIR=$(echo $line| awk -F'\\.git' '{print $1}' | awk -F/ '{print $NF}')
  # https://www.yeahshecodes.com/git/a-simple-git-subtree-tutorial
  $ECHO git subtree $1 --prefix=packages/$SUBDIR $line
done < $DIR/${BASH_SOURCE%.sh}.list

#!/bin/bash -ue

if [ $# -lt 1 ]
then
  echo "ERROR: Need one input argument: git operation (add/push/pull)"
  exit 3
fi

[[ "${@/echo}" == "$@" ]] && ECHO= || ECHO=echo

# https://gist.github.com/SKempin/b7857a6ff6bddb05717cc17a44091202
case $1 in
  add|pull) ARGS="--squash" ;;
  push)     ARGS=           ;;
  *)
    echo "ERROR: cannot handle operation '$1'"
    exit 3
  ;;
esac

DIR=$(cd $(dirname $BASH_SOURCE);pwd)
while IFS= read -r line || [[ -n "$line" ]]; do
  line=$(echo "$line" | sed 's:#.*::g')
  [ -z "$line" ] && continue
  # https://www.yeahshecodes.com/git/a-simple-git-subtree-tutorial
  $ECHO git subtree $1  $line $ARGS || true
done < $DIR/${BASH_SOURCE%.sh}.list

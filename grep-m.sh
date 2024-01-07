#!/bin/bash

DIRNOW=$( cd $(dirname $BASH_SOURCE); pwd )
project=$(
  [ -e "$DIRNOW/project.yaml" ] \
  && grep 'name:' "$DIRNOW/project.yaml" \
  || grep 'name:' "$DIRNOW/default.yaml" \
)
project=$(echo $project| awk -F\: '{print $2}')
project=${project// }

case "$1" in
  orbdir)
    grep file.orbdir\( "$DIRNOW/"*.m "$DIRNOW/metadata/$project/"*.yaml \
    | awk -F'file.orbdir' '{print $2}' \
    | awk -F"'" '{print $2}' \
    | sort -u
  ;;
  *)
    grep --color=always "$@" "$DIRNOW/"*.m "$DIRNOW/metadata/$project/"*.yaml
  ;;
esac

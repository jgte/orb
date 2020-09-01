#!/bin/bash

DIRNOW=$( cd $(dirname $BASH_SOURCE); pwd )
project=$(
  [ -e $DIRNOW/project.yaml ] \
  && grep 'name:' $DIRNOW/project.yaml \
  || grep 'name:' $DIRNOW/default.yaml \
)
project=$(echo $project| awk -F\: '{print $2}')
project=${project// }

grep --color=always "$@" $DIRNOW/*.m $DIRNOW/metadata/$project/*.yaml
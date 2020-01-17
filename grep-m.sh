#!/bin/bash

grep --color=always "$@" $( cd $(dirname $BASH_SOURCE); pwd )/*.m
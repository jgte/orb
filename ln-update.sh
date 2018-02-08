#!/bin/bash -ue

DIR=$(cd $(dirname $BASH_SOURCE);pwd)

ln -sfv $HOME/data/orb/data $DIR/data
ln -sfv $HOME/data/orb/plot $DIR/plot

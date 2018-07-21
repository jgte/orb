#!/bin/bash -ue

# Retiring data/metadata means it is no longer to be used in any way, forever (i.e. it is safe to delete)

METADATADIR=$(cd $(dirname BASH_SOURCE);pwd)
DATADIR=$(cd $METADATADIR/../data/;pwd)
RETIREDIR='retired'
OBSOLETEDIR='obsolete'
FIND_ARGS_SKIP_DIRS="-not -wholename \*$RETIREDIR\* -and -not -wholename \*$OBSOLETEDIR\*"


if [[ ! "$@" == "${@/echo/}" ]]
then
  ECHO=echo
else
  ECHO=
fi

if [[ ! "$@" == "${@/back/}" ]]
then
  BACK=true
else
  BACK=false
fi

for i in $@
do 
  case $i in 
    echo|back)
      #do nothing
    ;;
    *)
      if $BACK; then
        METADATA_LIST=$(find $METADATADIR/$RETIREDIR -name \*$i\* )
      else
        METADATA_LIST=$(find $METADATADIR $FIND_ARGS_SKIP_DIRS -name \*$i\* )
      fi
      for j in $METADATA_LIST
      do
        #gather source data dirs
        if $BACK; then
          #gather source data dirs
          DATA_LIST=$(find $DATADIR/$RETIREDIR $FIND_ARGS_SKIP_DIRS -type d -name $(basename ${j/.metadata})\*)
          #define sink dirs
              DATA_SINK_DIR=$DATADIR/
          METADATA_DINK_DIR=$METADATADIR
        else
          #make retired dirs, if needed
          [ -d     $DATADIR/$RETIREDIR ] || $ECHO mkdir -p     $DATADIR/$RETIREDIR
          [ -d $METADATADIR/$RETIREDIR ] || $ECHO mkdir -p $METADATADIR/$RETIREDIR
          #gather source data dirs
          DATA_LIST=$(find $DATADIR/ $FIND_ARGS_SKIP_DIRS -type d -name $(basename ${j/.metadata})\*)
          #define sink dirs
              DATA_SINK_DIR=$DATADIR/$RETIREDIR/ 
          METADATA_DINK_DIR=$METADATADIR/$RETIREDIR
        fi
        #move data dirs
        for k in $DATA_LIST
        do
          $ECHO mv -v $k $DATA_SINK_DIR || exit $?
        done
        #move metadata file
        $ECHO mv -v $j $METADATA_DINK_DIR
      done
      if ! $BACK; then
        #check for references to this metadata in remaining metadata files
        for j in $(find $METADATADIR -name \*$i\* -not -wholename \*$RETIREDIR\*)
        do
          OUT=$(grep -l $(basename ${j/.metadata}) *.metadata)
          [ -z "$OUT" ] || echo -e "'$(basename $j)' mentioned in:\n$OUT"
        done
      fi
    ;;
  esac
done
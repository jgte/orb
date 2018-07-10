#!/bin/bash -ue

# Retiring data/metadata means it is no longer to be used in any way, forever (i.e. it is safe to delete)

METADATADIR=$(cd $(dirname BASH_SOURCE);pwd)
DATADIR=$(cd $METADATADIR/../data/;pwd)
RETIREDIR='retired'

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
        METADATA_LIST=$(find $METADATADIR -name \*$i\* -not -wholename \*$RETIREDIR\*)
      fi
      for j in $METADATA_LIST
      do
        #gather source data dirs
        if $BACK; then
          #gather source data dirs
          DATA_LIST=$(find $DATADIR/$RETIREDIR -type d -name $(basename ${j/.metadata})\* )
          #define sink dirs
              DATA_SINK_DIR=$DATADIR/
          METADATA_DINK_DIR=$METADATADIR
        else
          #make retired dirs, if needed
          [ -d     $DATADIR/$RETIREDIR ] || $ECHO mkdir -p     $DATADIR/$RETIREDIR
          [ -d $METADATADIR/$RETIREDIR ] || $ECHO mkdir -p $METADATADIR/$RETIREDIR
          #gather source data dirs
          DATA_LIST=$(find $DATADIR/ -type d -name $(basename ${j/.metadata})\* -not -wholename \*$RETIREDIR\*)
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
        #check for references to this metadata in remaning metadata files
        for j in $(find $METADATADIR -name \*$i\* -not -wholename \*$RETIREDIR\*)
        do
          echo "== $(basename ${j/.metadata}) mentioned in:"
          grep -l $(basename ${j/.metadata}) *.metadata
        done
      fi
    ;;
  esac
done
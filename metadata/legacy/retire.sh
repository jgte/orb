#!/bin/bash -u

# Retiring data/metadata means it is no longer to be used in any way, forever (i.e. it is safe to delete)

METADATADIR=$(cd $(dirname BASH_SOURCE);pwd)
DATADIR=$(cd $METADATADIR/../data/;pwd)
RETIREDIR='retired'
OBSOLETEDIR='obsolete'
FIND_ARGS_SKIP_DIRS="-not -wholename \*$RETIREDIR\* -not -wholename \*$OBSOLETEDIR\*"

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
    enforce)
      for i in $METADATADIR/$RETIREDIR/*.metadata
      do
        $ECHO mv $i $METADATADIR/
	$ECHO $0 $(basename ${i/.metadata})
	$ECHO mv $METADATADIR/$(basename $i) $METADATADIR/$RETIREDIR/
      done
    ;;
    *)
      if $BACK; then
        METADATA_LIST=$(find $METADATADIR/$RETIREDIR $FIND_ARGS_SKIP_DIRS  -name \*$i\*)
      else
        METADATA_LIST=$(find $METADATADIR $FIND_ARGS_SKIP_DIRS -name \*$i\*  | grep -v "/$RETIREDIR/" | grep -v "/$OBSOLETEDIR/")
      fi
      [ -z "$ECHO" ] || echo METADATA_LIST=$METADATA_LIST
      for j in $METADATA_LIST
      do
        #gather source data dirs
        if $BACK; then
          #gather source data dirs
          DATA_LIST=$(find $DATADIR/$RETIREDIR -type d -name $(basename ${j/.metadata})\* | grep -v "/$RETIREDIR/" | grep -v "/$OBSOLETEDIR/")
          #define sink dirs
              DATA_SINK_DIR=$DATADIR/
          METADATA_DINK_DIR=$METADATADIR
        else
          #make retired dirs, if needed
          [ -d     $DATADIR/$RETIREDIR ] || $ECHO mkdir -p     $DATADIR/$RETIREDIR
          [ -d $METADATADIR/$RETIREDIR ] || $ECHO mkdir -p $METADATADIR/$RETIREDIR
          #gather source data dirs
          DATA_LIST=$(find $DATADIR/ -type d -name $(basename ${j/.metadata})\*  | grep -v "/$RETIREDIR/" | grep -v "/$OBSOLETEDIR/")
          #define sink dirs
              DATA_SINK_DIR=$DATADIR/$RETIREDIR/ 
          METADATA_DINK_DIR=$METADATADIR/$RETIREDIR
        fi
        [ -z "$ECHO" ] || echo DATA_LIST=$DATA_LIST
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
        for j in $(find $METADATADIR -name \*$i\* | grep -v "/$RETIREDIR/" | grep -v "/$OBSOLETEDIR/")
        do
          OUT=$(grep -l $(basename ${j/.metadata}) *.metadata)
          [ -z "$OUT" ] || echo -e "'$(basename $j)' mentioned in:\n$OUT"
        done
      fi
    ;;
  esac
done

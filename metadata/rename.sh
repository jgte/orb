#!/bin/bash -ue

# Retiring data/metadata means it is no longer to be used in any way, forever (i.e. it is safe to delete)

METADATADIR=$(cd $(dirname BASH_SOURCE);pwd)
DATADIR=$(cd $METADATADIR/../data/;pwd)
ECHO=

case $# in
  2)
    #do nothing
  ;;
  3)
    if [[ "$3" == "echo" ]] &&
    then
      ECHO=echo
    else
      echo "ERROR: if there are 3 input arguments, the third must be 'echo'"
      exit 3
    fi
  ;;
  *)
    echo "ERROR: need 2 or 3 input arguments, not $#"
    exit 3
  ;;
esac

METADATA=$(find $METADATADIR -name \*$1\* -not -wholename \*$RETIREDIR\*)
if [ $(echo "$METADATA" | wc -l) -eq 0 ]
then  
  echo "ERROR: cannot find any metadata called *$i*"
  exit 3
elif [ $(echo "$METADATA" | wc -l) -gt 1 ]; then
  echo "ERROR: there are several metadata called *$i*:"
  echo "$METADATA"
  echo "Can only handle one metadata at a time"
  exit 3
fi

#loop over source data dir(s)
for i in $(find $DATADIR/ -type d -name $(basename ${METADATA/.metadata})\*)
do
  #define sink dir
  DATA_SINK_DIR=$(basename $i)
  DATA_SINK_DIR=${DATA_SINK_DIR/$1/$2}
  #rename data dir
  $ECHO mv -v $i $DATADIR/$i
done

#

continuar aqui

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
    $ECHO mv -vf $k $DATA_SINK_DIR || exit $?
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

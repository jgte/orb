#!/bin/bash -ue

# making data obsolete means that it's definition/processing has changed and the old version of the data should not be used any more.

METADATADIR=$(cd $(dirname BASH_SOURCE);pwd)
DATADIR=$(cd $METADATADIR/../data/;pwd)
OBSOLETEDIR='obsolete'

GSTAT=$(which gstat || which stat)

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
      METADATA_LIST=$(find $METADATADIR -name \*$i\* -not -wholename \*$OBSOLETEDIR\*)
      for j in $METADATA_LIST
      do
        #get list of files
        if $BACK; then
          DATA_LIST=$(find $DATADIR/$OBSOLETEDIR/ -type d -name $(basename ${j/.metadata})\*)
        else
          DATA_LIST=$(find $DATADIR/ -not -wholename \*$OBSOLETEDIR\* -type d -name $(basename ${j/.metadata})\*)
        fi
        #move data dirs
        for k in $DATA_LIST
        do
          if [ -z "$(ls $k)" ]
          then
            #if the directory is empty, then just get rid of it
            echo "Directory $k is empty, removing it..."
            $ECHO rmdir $k || exit $?
          else
            if $BACK; then
              #move back to the data dir
              SINK=$DATADIR/$(basename ${k/.metadata})/
              #update k
              k=$(find $k -maxdepth 1 -type d | sort | tail -n1)
            else
              #determine latest date of the files in this dir
              DATA_VERSION=$($GSTAT -c '%y' $(ls -t $k/*.mat | head -n1) | awk '{print $1}')
              #build sink dir
              SINK=$DATADIR/$OBSOLETEDIR/$(basename ${k/.metadata})/$DATA_VERSION
            fi
            #make destination
            [ -d $SINK ] || $ECHO mkdir -p $SINK
            if $BACK && [ -z "$(ls $k/*.mat 2> /dev/null)" ]
            then
              #if the directory is empty, then just get rid of it
              echo "Directory $k is empty, removing it..."
              $ECHO rmdir $k || exit $?
            else
              $ECHO mv -v $k/*.mat $SINK/ || exit $?
            fi
          fi
        done
      done 
      #check for references to this metadata in remaning metadata files
      for j in $METADATA_LIST
      do
        OUT=$(grep -l $(basename ${j/.metadata}) *.metadata)
        [ -z "$OUT" ] ||echo "$(basename $j) mentioned in:"
      done
    ;;
  esac
done
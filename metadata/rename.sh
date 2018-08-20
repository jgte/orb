#!/bin/bash -ue

# renames data and metadata

METADATADIR=$(cd $(dirname BASH_SOURCE);pwd)
DATADIR=$(cd $METADATADIR/../data/;pwd)
RETIREDIR='retired'
OBSOLETEDIR='obsolete'

ECHO=

case $# in
  2)
    #do nothing
  ;;
  3)
    if [[ "$3" == "echo" ]]
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

METADATA_OLD=$(find $METADATADIR -name \*$1\* | grep -v "/$RETIREDIR/" | grep -v "/$OBSOLETEDIR/")
METADATA_NEW=$METADATADIR/$(basename ${METADATA_OLD/$1/$2})

if [ $(echo "$METADATA_OLD" | wc -l) -eq 0 ]
then  
  echo "ERROR: cannot find any metadata called *$METADATA_OLD*"
  exit 3
elif [ $(echo "$METADATA_OLD" | wc -l) -gt 1 ]; then
  echo "ERROR: there are several metadata called *$METADATA_OLD*:"
  echo "$METADATA_OLD"
  echo "Can only handle one metadata at a time"
  exit 3
fi

# if [ -z "$ECHO" ]
# then
#   read -p "Renaming $(basename $METADATA_OLD) to $(basename $METADATA_NEW)? [Y/n]" ANSWER
#   if [ "$ANSWER" == "N" ] || [ "$ANSWER" == "n" ]
#   then
#     echo "Nothing done..."
#     exit
#   fi
# fi


#define the old and new metadata roots
MTDR_OLD=$(basename ${METADATA_OLD/.metadata})
MTDR_NEW=$(basename ${METADATA_NEW/.metadata})

if [ -e $DATADIR/$MTDR_NEW ]
then
  echo "ERROR: directory $DATADIR/$MTDR_NEW exists, cannot renamed"
  exit 3
fi

#make new data dir
$ECHO mkdir -p $DATADIR/$MTDR_NEW

#rename all files in this new dir
for j in $(find $DATADIR/$MTDR_OLD -name $MTDR_OLD\* -not -type d)
do
  #rename data 
  $ECHO mv -v $j $DATADIR/$MTDR_NEW/$(basename ${j//$MTDR_OLD/$MTDR_NEW})
done

#rename metadata
$ECHO mv -v $METADATA_OLD $METADATA_NEW

#remove old dir
rmdir $DATADIR/$MTDR_OLD 

#check for references to this metadata in remaining metadata files
OUT=$(grep -l $(basename ${MTDR_OLD/.metadata}) *.metadata)
[ -z "$OUT" ] || echo -e "'$MTDR_OLD' mentioned in:\n$OUT"

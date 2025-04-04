#!/bin/bash -ue

# renames data and metadata

DIRPROJECT=$(cd $(dirname BASH_SOURCE)/..;pwd)
DATADIR=$(cd $DIRPROJECT/data/;pwd)
PLOTDIR=$(cd $DIRPROJECT/plot/;pwd)
[ -e $DIRPROJECT/project.yaml ] \
&& PROJECT_METADATA=$DIRPROJECT/project.yaml \
|| PROJECT_METADATA=$DIRPROJECT/default.yaml
METADATADIR=$DIRPROJECT
METADATADIR+="/$(awk '/metadata_dir:/ {print $2}' $PROJECT_METADATA)"
METADATADIR+="/$(awk '/name:/         {print $2}' $PROJECT_METADATA)"
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

if [ -z "$ECHO" ]
then
  read -p "Renaming $(basename $METADATA_OLD) to $(basename $METADATA_NEW)? [Y/n]" ANSWER
  if [ "$ANSWER" == "N" ] || [ "$ANSWER" == "n" ]
  then
    echo "Nothing done..."
    exit
  fi
fi

#define the old and new metadata roots
MTDR_OLD=$(basename ${METADATA_OLD/.yaml})
MTDR_NEW=$(basename ${METADATA_NEW/.yaml})

if [ -e $DATADIR/$MTDR_NEW ]
then
  echo "ERROR: directory $DATADIR/$MTDR_NEW exists, cannot renamed"
  exit 3
fi

for dir in $DATADIR $PLOTDIR
do
  #make new data dir
  $ECHO mkdir -p $dir/$MTDR_NEW
  #rename all files in this new dir
  for j in $(find $dir/$MTDR_OLD -name $MTDR_OLD\* -not -type d || echo "NOTICE: no data found for $MTDR_OLD, ignoring data move" 1>&2)
  do
    #rename data
    $ECHO mv -v $j $dir/$MTDR_NEW/$(basename ${j//$MTDR_OLD/$MTDR_NEW})
  done
  #remove old dir
  $ECHO rmdir $dir/$MTDR_OLD || echo "NOTICE: no data found for $MTDR_OLD, ignoring dir delete" 1>&2
done

#rename metadata
$ECHO mv -v $METADATA_OLD $METADATA_NEW

#check for references to this metadata in remaining metadata files
OUT=$(grep -l $(basename ${MTDR_OLD/.yaml}) $METADATADIR/*.yaml)
[ -z "$OUT" ] || $ECHO file-find-replace.sh -from=$MTDR_OLD -to=$MTDR_NEW $OUT

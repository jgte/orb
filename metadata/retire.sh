#!/bin/bash -ue

# Retiring data/metadata means it is no longer to be used in any way, forever (i.e. it is safe to delete)

METADATA_EXT=yaml

DIRPROJECT=$(cd $(dirname BASH_SOURCE)/..;pwd)
DATADIR=$(cd $DIRPROJECT/data/;pwd)
PLOTDIR=$(cd $DIRPROJECT/plot/;pwd)
[ -e $DIRPROJECT/project.yaml ] \
&& PROJECT_METADATA=$DIRPROJECT/project.yaml \
|| PROJECT_METADATA=$DIRPROJECT/default.yaml
SOURCE_METADIR=$DIRPROJECT
SOURCE_METADIR+="/$(awk '/metadata_dir:/ {print $2}' $PROJECT_METADATA)"
SOURCE_METADIR+="/$(awk '/name:/         {print $2}' $PROJECT_METADATA)"
SOURCE_DATADIR=$DIRPROJECT
SOURCE_DATADIR+="/$(awk '/^data_dir:/ {print $2}' $PROJECT_METADATA)"
SOURCE_PLOTDIR=$DIRPROJECT
SOURCE_PLOTDIR+="/$(awk '/plot_dir:/ {print $2}' $PROJECT_METADATA)"

ECHO=
VERBOSE=false
BACK=false
DELETE=true
RETIREDIR=$PWD/retired
IGNOREDIRS=(obsolete retired to-delete)
ARGS=()
INCLUDE_DATA=true
INCLUDE_PLOT=true
INCLUDE_META=true
for i in $@
do 
  case $i in 
    -echo|echo)             ECHO=echo          ;;
    -verbose|verbose)       VERBOSE=true       ;;
    -back|back)             BACK=true          ;;
    -no-data)               INCLUDE_DATA=false ;;
    -no-plot|-no-plots)     INCLUDE_PLOT=false ;;
    -no-metadata)           INCLUDE_META=false ;;
    -no-delete|-keep|-copy) DELETE=false       ;;
    -retired-dir=*)  RETIREDIR=${i/-retired-dir=} ;;
    -ignore-dirs=*) IGNOREDIRS=${i/-ignore-dirs=}; IGNOREDIRS=(${IGNOREDIRS//,/ }) ;;
    *) ARGS+=($i) ;;
  esac
done

#define copy/move operation
$DELETE \
&& OP="mv -v" \
|| OP="cp -av"

$VERBOSE && echo "\
ECHO       : $ECHO
BACK       : $BACK
DELETE     : $DELETE
RETIREDIR  : $RETIREDIR
IGNOREDIRS : $IGNOREDIRS
OP         : $OP
"

#sanity
if [ ${#ARGS[@]} -eq 0 ]
then
  echo "ERROR: need a valid metadata filename part"
  exit 3
fi
#more sanity: back and enforce cannot go together
if $BACK && [[ ! "${ARGS[@]/enforce}" == "${ARGS[@]}" ]]
then
  echo "ERROR: input options 'force' and 'back' are incomptible."
  exit 3
fi

#building file arguments on the ignore dirs
if [ ${#IGNOREDIRS[@]} -gt 0 ]
then
  FIND_ARGS_SKIP_DIRS=()
  for i in "${IGNOREDIRS[@]}"
  do
    FIND_ARGS_SKIP_DIRS+=("-not -wholename *$i*")
  done
  FIND_ARGS_SKIP_DIRS="${FIND_ARGS_SKIP_DIRS[@]}"
else
  FIND_ARGS_SKIP_DIRS=
fi
$VERBOSE && echo "\
FIND_ARGS_SKIP_DIRS : ${FIND_ARGS_SKIP_DIRS[@]}
"

#build and create sink dirs
SINK_METADIR=$RETIREDIR/metadata/$(echo "$SOURCE_METADIR" | awk -F/ '{print $NF}')
SINK_DATADIR=$RETIREDIR/data
SINK_PLOTDIR=$RETIREDIR/plot

#honour the back flag
if $BACK
then
     TMP_METADIR=$SINK_METADIR
     TMP_DATADIR=$SINK_DATADIR
     TMP_PLOTDIR=$SINK_PLOTDIR
    SINK_METADIR=$SOURCE_METADIR
    SINK_DATADIR=$SOURCE_DATADIR
    SINK_PLOTDIR=$SOURCE_PLOTDIR
  SOURCE_METADIR=$TMP_METADIR
  SOURCE_DATADIR=$TMP_DATADIR
  SOURCE_PLOTDIR=$TMP_PLOTDIR
fi
$VERBOSE && echo "\
SINK_METADIR   : $SINK_METADIR
SINK_DATADIR   : $SINK_DATADIR
SINK_PLOTDIR   : $SINK_PLOTDIR
SOURCE_METADIR : $SOURCE_METADIR
SOURCE_DATADIR : $SOURCE_DATADIR
SOURCE_PLOTDIR : $SOURCE_PLOTDIR
"

#make sure sources exist
for i in $SOURCE_METADIR $SOURCE_DATADIR $SOURCE_PLOTDIR
do
  if [ ! -d "$i" ]
  then
    echo "ERROR: cannot find source directory $i"
    exit 3
  fi
done
#create sinks if they don't exist
$INCLUDE_META && [ ! -d "$SINK_METADIR" ] && $ECHO mkdir -p $SINK_METADIR
$INCLUDE_DATA && [ ! -d "$SINK_DATADIR" ] && $ECHO mkdir -p $SINK_DATADIR
$INCLUDE_PLOT && [ ! -d "$SINK_PLOTDIR" ] && $ECHO mkdir -p $SINK_PLOTDIR

for i in "${ARGS[@]}"
do 
  case $i in 
    enforce)
      for i in $SINK_METADIR/*.$METADATA_EXT
      do
        $ECHO mv $i $SOURCE_METADIR/
	      $ECHO $BASH_SOURCE $(basename ${i/.$METADATA_EXT})
      done
      if [ ! "${ARGS[@]}" == "enforce" ]
      then
        echo -n "WARNING: when option 'enforce' is given, no other metadata operations are allowed; "
        echo "ignored the following arguments: ${ARGS[@]/enforce}"
      fi
      exit
    ;;
    *)
      $VERBOSE && echo "----- handling metadata filename part: $i"
      METADATA_LIST=$(find $SOURCE_METADIR/ ${FIND_ARGS_SKIP_DIRS[@]} -name \*$i\*)
      for j in $METADATA_LIST
      do
        $VERBOSE && echo "----- moving metadata file: $(basename $j)"
        #move metadata file
        $INCLUDE_META && $ECHO $OP $j $SINK_METADIR
        #move data dirs
        if $INCLUDE_DATA 
        then
          #gather source data dirs
          DATA_LIST=$(find $SOURCE_DATADIR/ ${FIND_ARGS_SKIP_DIRS[@]} -type d -name $(basename ${j/.$METADATA_EXT}))
          for k in $DATA_LIST
          do
            $VERBOSE && echo "----- moving data dir: $(basename $k)"
            $ECHO $OP $k $SINK_DATADIR || exit $?
          done
        fi
        #move PLOT dirs
        if $INCLUDE_PLOT
        then
          #gather source plot dirs
          PLOT_LIST=$(find $SOURCE_PLOTDIR/ ${FIND_ARGS_SKIP_DIRS[@]} -type d -name $(basename ${j/.$METADATA_EXT}))
          for k in $PLOT_LIST
          do
            $VERBOSE && echo "----- moving plot dir: $(basename $k)"
            $ECHO $OP $k $SINK_PLOTDIR || exit $?
          done
        fi
      done
      if ! $BACK; then
        #check for references to this metadata in remaining metadata files
        for j in $METADATA_LIST
        do
          OUT=$(grep -l $(basename ${j/.$METADATA_EXT}) $SOURCE_METADIR/*.$METADATA_EXT)
          [ -z "$OUT" ] || echo -e "'$(basename $j)' mentioned in:\n$OUT"
        done
      fi
    ;;
  esac
done

#!/bin/bash -u

DIR=$(cd $(dirname $BASH_SOURCE);pwd)

#inits
SOURCE=
SINK=
ECHO=
for i in "$@"
do
  case "$i" in
    help)
      echo "Create links between metadata subdirs. Options are:"
      grep ') #' $BASH_SOURCE \
        | grep -v grep \
        | grep -v sed \
        | sed 's:) #:|:g' \
        | column -t -s \|
      exit
    ;;
    -x) #set bash's -x
      set -x
    ;;
    echo) #show what would have been done but don't do anything
      ECHO=echo
    ;;
    source=*) #mandatory: define the directory to which the links will point to
      SOURCE=${i/source=}
      if [ ! -d "$SOURCE" ]
      then
        echo "ERROR: cannot find dir '$SOURCE'"
        exit 3
      fi
    ;;
    sink=*) #mandatory: define the directory inside which the links will reside
      SINK=${i/sink=}
      if [ ! -d "$SINK" ]
      then
        echo "ERROR: cannot find dir '$SINK'"
        exit 3
      fi
    ;;
    *)
      echo "ERROR: cannot handle input '$i'"
    ;;
  esac
done

if [ -z "$SOURCE" ]
then
  echo "ERROR: need input 'source=...; call '$BASH_SOURCE help' for more info"
  exit 3
fi
if [ -z "$SINK" ]
then
  echo "ERROR: need input 'sink=...; call '$BASH_SOURCE help' for more info"
  exit 3
fi


cd $SINK

for i in *.yaml
do
  j="../$SOURCE/$(basename $i)"
  if [ -e "$j" ] && diff -q "$i" "$j" > /dev/null
  then
    $ECHO ln -sfv "../$(basename $SOURCE)/$(basename $i)" .
  fi
done

cd - > /dev/null
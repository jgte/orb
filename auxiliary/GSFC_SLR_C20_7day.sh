#!/bin/bash -ue

DIR=$(cd $(dirname $BASH_SOURCE);pwd)
DATA_FILE=gsfc_slr_5x5c61s61.txt
URL="https://earth.gsfc.nasa.gov/sites/default/files/geo/slr-weekly/$DATA_FILE"
DATE=$(which gdate 2> /dev/null || which date)
TOUCH=$(which gtouch 2> /dev/null || which touch)

#download data file if needed
DATE_FILE=/tmp/$$.30daysago
$TOUCH -d "$($DATE -d '30 days ago')" $DATE_FILE
[ "$DIR/$DATA_FILE" -ot $DATE_FILE ] && wget $URL -O $DIR/$DATA_FILE
rm -f $DATE_FILE

awk '
/^[0-9]/ {printf("%s ",$2)}
/^  2  0 /{printf("%s\n",$3)}
' $DIR/$DATA_FILE \
> "$DIR/GSFC_SLR_C20_7day.txt"
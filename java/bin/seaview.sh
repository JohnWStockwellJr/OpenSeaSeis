#!/bin/bash
# Launch Seismic Viewer (Seaview)

# Increase if needed (when Out-of-Memory problems occur):
MEMSIZE=3824

DIR=$(dirname $0)
# Capture the case when script is executed from current directory
LETTER=$(echo $DIR | cut -c 1)
if [ $LETTER == '.' ]; then
   DIR=$(pwd)/$DIR
fi

java -cp .:${DIR}/../lib/CSeisLib.jar:${DIR}/../lib/SeaView.jar -Xmx${MEMSIZE}m -Djava.library.path=${DIR}/../lib cseis/seaview/SeaView  $*

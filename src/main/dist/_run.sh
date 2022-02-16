#!/usr/bin/env bash
# shell script to run object-mapper-pipeline
. /etc/profile

APPNAME="object-mapper-pipeline"
APPDIR=/home/rgddata/pipelines/$APPNAME

cd $APPDIR
java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$APPDIR/properties/log4j2.xml \
    -jar lib/$APPNAME.jar "$@"

#!/usr/bin/env bash
# shell script to run ObjectMapper pipeline
. /etc/profile

APPNAME=ObjectMapper
APPDIR=/home/rgddata/pipelines/$APPNAME

cd $APPDIR
java -Dspring.config=$APPDIR/../properties/default_db.xml \
    -Dlog4j.configuration=file://$APPDIR/properties/log4j.properties \
    -jar lib/$APPNAME.jar "$@"

#!/usr/bin/env bash
# shell script to run ObjectMapper pipeline
. /etc/profile

APPNAME=ObjectMapper
APPDIR=/home/rgddata/pipelines/$APPNAME

cd $APPDIR
echo ""
DB_OPTS="-Dspring.config=$APPDIR/../properties/default_db.xml"
LOG4J_OPTS="-Dlog4j.configuration=file://$APPDIR/properties/log4j.properties"
export OBJECT_MAPPER_OPTS="$DB_OPTS $LOG4J_OPTS"

bin/$APPNAME "$@"

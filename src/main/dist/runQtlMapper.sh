# analyze all strain markers and compute qtl positions for rat, mouse and human
APPHOME=/home/rgddata/pipelines/ObjectMapper
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAIL_LIST=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST=mtutaj@mcw.edu,sjwang@mcw.edu,jrsmith@mcw.edu
fi

$APPHOME/_run.sh -qtls
mailx -s "[$SERVER] QTL Mapper report" $EMAIL_LIST < $APPHOME/logs/qtl_mapper.log
#mailx -s "[$SERVER] QTL Mapper position changes" $EMAIL_LIST < $APPHOME/logs/position_changes.log

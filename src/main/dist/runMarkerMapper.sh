# analyze rat DB_SNP markers and compute/validate their genomic positions on assemblies 3.4, 5.0, 6.0 and 7.2
APPHOME=/home/rgddata/pipelines/object-mapper-pipeline
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAIL_LIST=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST=mtutaj@mcw.edu,sjwang@mcw.edu,jrsmith@mcw.edu
fi

$APPHOME/_run.sh -markers -rat
mailx -s "[$SERVER] Marker Mapper report" $EMAIL_LIST < $APPHOME/logs/marker_mapper.log

# analyze rat DB_SNP markers and compute/validate their genomic positions on assemblies 3.4, 5.0 and 6.0
APPHOME=/home/rgddata/pipelines/ObjectMapper
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAIL_LIST=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST=mtutaj@mcw.edu,sjwang@mcw.edu,jrsmith@mcw.edu
fi

$APPHOME/run.sh -sslps -rat
mailx -s "[$SERVER] SSLP Mapper report" $EMAIL_LIST < $APPHOME/logs/sslp_mapper.log

# generate a report file for qtls with all out-of-region genes
#
APPHOME=/home/rgddata/pipelines/ObjectMapper
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAIL_LIST=mtutaj@mcw.edu
if [ "$SERVER" = "REED" ]; then
  EMAIL_LIST=mtutaj@mcw.edu,sjwang@mcw.edu,jrsmith@mcw.edu
fi

TODAY=`date +"%Y-%m-%d"`
$APPHOME/run.sh -qtls -reportOutOfRegionGenes
mailx -s "[$SERVER] QTL Mapper out-of-region genes" $EMAIL_LIST < $APPHOME/data/qtlsWithOutOfRegionGenes_${TODAY}.txt

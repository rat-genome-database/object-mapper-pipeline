# compute strain positions on reference assembly based on their marker positions
APPHOME=/home/rgddata/pipelines/object-mapper-pipeline
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAIL_LIST=mtutaj@mcw.edu

$APPHOME/_run.sh -strains -rat
mailx -s "[$SERVER] Strain Mapper report" $EMAIL_LIST < $APPHOME/logs/strain_mapper_summary.log

package edu.mcw.rgd.dataload.ObjectMapper;

import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.datamodel.SSLP;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;


/**
 * analyze rat DB_SNP sslps and recompute their positions for RGSC 3.4, Rnor5.0 and Rnor6.0 assemblies
 */
public class SslpMapper extends BaseMapper {
    Logger log = Logger.getLogger("sslpMapper");

    public void run(int speciesType) throws Exception {
        long startTime = System.currentTimeMillis();

        log.info(getDao().getConnectionInfo());
        log.info(getVersion());

        // only rat is supported
        if( speciesType!= SpeciesType.RAT ) {
            String msg = "SslpMapper can only run for rat! You tried to run it for "+SpeciesType.getCommonName(speciesType);
            log.error(msg);
            return;
        }

        // get list of all strains for given species
        List<SSLP> sslps = dao.getActiveSSLPsForDbSnp(speciesType);
        log.info("number of active DB_SNP sslps: "+sslps.size());

        // run processing for every map key specified in property file
        for( Integer mapKey: getMapKeysForSpecies(speciesType) ) {
            run(sslps, mapKey);
        }
        log.info("OK -- elapsed "+ Utils.formatElapsedTime(startTime, System.currentTimeMillis()));
    }

    public void run(List<SSLP> sslps, int mapKey) throws Exception {

        // counters
        int matchingMapPositions = 0;
        int deletedMapPositions = 0;
        int insertedMapPositions = 0;

        // process all sslps
        for( SSLP sslp: sslps ) {

            // compare the map positions with existing map positions in rgd
            List<MapData> sslpPositionsInRgd = dao.getMapPositions(sslp.getRgdId(), mapKey);
            List<MapData> dbSnpPositionsInRgd = dao.getDbSnpPositions(sslp.getName(), mapKey, sslp.getRgdId());

            // read positions from database
            List<MapData> mdsToBeInserted = new ArrayList<>();
            List<MapData> mdsToBeDeleted = new ArrayList<>();

            int matchingPositions = qcMapData(dbSnpPositionsInRgd, sslpPositionsInRgd, mdsToBeInserted, mdsToBeDeleted);

            updateMapData(mdsToBeInserted, mdsToBeDeleted);

            matchingMapPositions += matchingPositions;
            insertedMapPositions += mdsToBeInserted.size();
            deletedMapPositions += mdsToBeDeleted.size();
        }

        // display statistics
        if( matchingMapPositions>0 ) {
            log.info("map_key="+mapKey+": matching db_snp sslp map positions: "+matchingMapPositions);
        }
        if( insertedMapPositions>0 ) {
            log.info("map_key="+mapKey+": inserted db_snp sslp map positions: "+insertedMapPositions);
        }
        if( deletedMapPositions>0 ) {
            log.info("map_key="+mapKey+": deleted db_snp sslp map positions: "+deletedMapPositions);
        }
    }

}
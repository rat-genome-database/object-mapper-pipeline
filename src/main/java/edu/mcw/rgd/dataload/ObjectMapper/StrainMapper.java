package edu.mcw.rgd.dataload.ObjectMapper;

import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.Strain;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * @author mtutaj
 * @since Mar 11, 2011
 * Computes strain positions on given assemblies based on position of flanking markers
 */
public class StrainMapper extends BaseMapper {

    Logger log = Logger.getLogger("strainMapper");

    public void run(int speciesType) throws Exception {
        log.info(getDao().getConnectionInfo());

        long time0 = System.currentTimeMillis();

        // only rat is supported
        if( speciesType!= SpeciesType.RAT ) {
            log.warn("StrainMapper can only be run for rat");
            System.exit(1);
        }
        log.info(getVersion());

        // get list of all strains for given species
        List<Strain> strains = dao.getActiveStrains(speciesType);
        log.info("number of active strains: "+strains.size());

        // run processing for every map key specified in property file
        for( Integer mapKey: getMapKeysForSpecies(speciesType) ) {
            run(strains, mapKey);
        }
        log.info("OK! - elapsed "+Utils.formatElapsedTime(time0, System.currentTimeMillis()));
    }

    public void run(List<Strain> strains, int mapKey) throws Exception {

        // counters
        int matchingMapPositions = 0;
        int deletedMapPositions = 0;
        int insertedMapPositions = 0;

        // process all strains
        for( Strain strain: strains ) {

            // get validated rgd ids for flanking markers
            Map<String, List<Integer>> markers;
            String strainType = Utils.defaultString(strain.getStrainTypeName());
            if( strainType.equals("mutant")
                    || strainType.equals("mutant_knockin")
                    || strainType.equals("mutant_knockout") ) {
                markers = dao.getRgdIdsForMutantStrains(strain.getRgdId());
            } else {
                markers = dao.getRgdIdsForFlankingMarkers(strain.getRgdId());
            }

            // compare the map positions with existing map positions in rgd
            List<MapData> strainMapPositions = dao.getMapPositions(strain.getRgdId(), mapKey);

            // no data to process?
            if( markers.isEmpty() && strainMapPositions.isEmpty() ) {
                continue;
            }

            // combine positions within regions
            List<MapData> incomingMapPositions = new ArrayList<>();
            for( String region: markers.keySet() ) {
                // process all map positions
                List<MapData> sslpMapPositions = new ArrayList<>(2); // combined map positions
                for( Integer markerRgdId: markers.get(region) ) {
                    List<MapData> mapPosList = dao.getMapPositions(markerRgdId, mapKey);
                    for( MapData md: mapPosList ) {
                        combineMapPositions(sslpMapPositions, md, strain.getRgdId(), region);
                    }
                }
                incomingMapPositions.addAll(sslpMapPositions);
            }

            // read positions from database
            List<MapData> mdsToBeInserted = new ArrayList<>();
            List<MapData> mdsToBeDeleted = new ArrayList<>();

            int matchingPositions = qcMapData(incomingMapPositions, strainMapPositions, mdsToBeInserted, mdsToBeDeleted);

            updateMapData(mdsToBeInserted, mdsToBeDeleted);

            matchingMapPositions += matchingPositions;
            insertedMapPositions += mdsToBeInserted.size();
            deletedMapPositions += mdsToBeDeleted.size();
        }

        // display statistics
        if( matchingMapPositions>0 ) {
            log.info("map_key="+mapKey+": matching strain map positions: "+matchingMapPositions);
        }
        if( insertedMapPositions>0 ) {
            log.info("map_key="+mapKey+": inserted strain map positions: "+insertedMapPositions);
        }
        if( deletedMapPositions>0 ) {
            log.info("map_key="+mapKey+": deleted  strain map positions: "+deletedMapPositions);
        }
    }

    void combineMapPositions(List<MapData> mds, MapData mdSslp, int strainRgdId, String region) {

        // the map position to be combined must have a chromosome and start & stop positions
        if( mdSslp.getChromosome()==null || mdSslp.getChromosome().length()==0 ||
            mdSslp.getStartPos()==null || mdSslp.getStartPos()<0 ||
            mdSslp.getStopPos()==null || mdSslp.getStopPos()<0 ) {
            log.error("invalid position for sslp_rgd_id="+mdSslp.getRgdId());
            return;
        }

        // check if chromosomes match
        for( MapData mdStrain: mds ) {
            // combine new position with existing position for same chromosome
            if( mdStrain.getChromosome().equals(mdSslp.getChromosome()) ) {
                // combine positions
                if( mdSslp.getStartPos() < mdStrain.getStartPos() )
                    mdStrain.setStartPos(mdSslp.getStartPos());
                if( mdSslp.getStopPos() > mdStrain.getStopPos() )
                    mdStrain.setStopPos(mdSslp.getStopPos());
                mdStrain.setNotes(region);
                return;
            }
        }

        // if positions have not been combined yet, add another map entry
        MapData mdStrain = new MapData();
        mdStrain.setMapKey(mdSslp.getMapKey());
        mdStrain.setRgdId(strainRgdId);
        mdStrain.setChromosome(mdSslp.getChromosome());
        mdStrain.setStartPos(mdSslp.getStartPos());
        mdStrain.setStopPos(mdSslp.getStopPos());
        mdStrain.setMapsDataPositionMethodId(1); // position by flanking markers
        mdStrain.setSrcPipeline(getSrcPipeline());
        mdStrain.setNotes(region);
        mds.add(mdStrain);
    }
}

package edu.mcw.rgd.dataload.ObjectMapper;

import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.Strain;
import edu.mcw.rgd.process.Utils;

import java.util.*;

/**
 * @author mtutaj
 * @since Mar 11, 2011
 * Computes strain positions on given assemblies based on position of flanking markers
 */
public class StrainMapper extends BaseMapper {

    public void run(int speciesType) throws Exception {

        // only rat is supported
        if( speciesType!= SpeciesType.RAT ) {
            log.warn("StrainMapper can only be run for rat");
            return;
        }

        // get list of all strains for given species
        List<Strain> strains = dao.getActiveStrains(speciesType);
        Collections.shuffle(strains);
        log.info("number of active strains: "+strains.size());

        // run processing for every map key specified in property file
        for( Integer mapKey: getMapKeysForSpecies(speciesType) ) {
            run(strains, mapKey);
        }
    }

    public void run(List<Strain> strains, int mapKey) throws Exception {

        // counters
        int matchingMapPositions = 0;
        int deletedMapPositions = 0;
        int insertedMapPositions = 0;
        int methodIdUpdated = 0;

        // process all strains
        for( Strain strain: strains ) {

            String strainType = Utils.defaultString(strain.getStrainTypeName());
            boolean isMutantStrain = strainType.equals("mutant")
                    || strainType.equals("mutant_knockin")
                    || strainType.equals("mutant_knockout");

            // get validated rgd ids for flanking markers
            Map<String, List<Integer>> markers;
            if( isMutantStrain ) {
                markers = dao.getRgdIdsForMutantStrains(strain.getRgdId());
            } else {
                markers = dao.getRgdIdsForFlankingMarkers(strain.getRgdId());
            }

            // compare the map positions with existing map positions in rgd
            //
            int positionMethodId = 1; // by flanking markers
            if( isMutantStrain ) {
                positionMethodId = 7; // position from mutated gene
            }

            List<MapData> strainMapPositions = dao.getMapPositions(strain.getRgdId(), mapKey);

            List<MapData> allellicVariantsPositions = Collections.emptyList();
            if( isMutantStrain ) {
                allellicVariantsPositions = dao.getAllelicVariantPositionsForStrain(strain.getRgdId(), mapKey);
                if( allellicVariantsPositions!=null && allellicVariantsPositions.size()>0 ) {
                    positionMethodId = 8; // Position from associated allelic variant
                }
            }

            // no data to process?
            if( markers.isEmpty() && strainMapPositions.isEmpty() && allellicVariantsPositions.isEmpty() ) {
                continue;
            }

            // combine positions within regions
            List<MapData> incomingMapPositions = new ArrayList<>();
            if( positionMethodId == 8 ) {
                // for allelic variants: just use the position as-is
                // and ensure the rgd id is correct
                for( MapData md: allellicVariantsPositions ) {

                    MapData md2 = md.clone();
                    md2.setRgdId( strain.getRgdId() );
                    md2.setKey(0);
                    md2.setSrcPipeline( getSrcPipeline() );
                    md2.setStrand(null);
                    md2.setNotes("allelic variant RGD ID: "+md.getRgdId());
                    md2.setMapsDataPositionMethodId(positionMethodId);

                    incomingMapPositions.add(md2);
                }
            } else {
                for (String region : markers.keySet()) {
                    // process all map positions
                    List<MapData> sslpMapPositions = new ArrayList<>(2); // combined map positions
                    for (Integer markerRgdId : markers.get(region)) {
                        List<MapData> mapPosList = dao.getMapPositions(markerRgdId, mapKey);
                        for (MapData md : mapPosList) {
                            combineMapPositions(sslpMapPositions, md, strain.getRgdId(), region, positionMethodId);
                        }
                    }
                    incomingMapPositions.addAll(sslpMapPositions);
                }
            }

            // read positions from database
            List<MapData> mdsToBeInserted = new ArrayList<>();
            List<MapData> mdsToBeDeleted = new ArrayList<>();

            // ensure that incoming positions have the correct mapping method
            for( MapData md: incomingMapPositions ) {
                if( !Utils.intsAreEqual(md.getMapsDataPositionMethodId(), positionMethodId) ) {
                    md.setMapsDataPositionMethodId(positionMethodId);
                    methodIdUpdated++;
                }
            }

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
        if( methodIdUpdated>0 ) {
            log.warn("WARN (unexpected):   map_key="+mapKey+": method id updated: "+methodIdUpdated);
        }
    }

    void combineMapPositions(List<MapData> mds, MapData mdSslp, int strainRgdId, String region, int positionMethodId) {

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
        mdStrain.setMapsDataPositionMethodId(positionMethodId);
        mdStrain.setSrcPipeline(getSrcPipeline());
        mdStrain.setNotes(region);
        mds.add(mdStrain);
    }
}

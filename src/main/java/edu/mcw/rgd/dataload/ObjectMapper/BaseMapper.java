package edu.mcw.rgd.dataload.ObjectMapper;

import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.text.SimpleDateFormat;
import java.util.*;

/**
 * @author mtutaj
 * @since Mar 11, 2011
 * base class for all mappers
 */
public abstract class BaseMapper {

    DAO dao;
    Logger log;
    private Map<String,String> mapKeys;
    private String version;
    private String logName;
    private String srcPipeline;
    private Map<String,String> params;

    public void start(int speciesType) throws Exception {
        long time0 = System.currentTimeMillis();
        log.info("   "+getVersion());
        log.info("   "+getDao().getConnectionInfo());

        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        log.info("   started at "+sdt.format(new Date()));

        run(speciesType);

        log.info("OK! - elapsed "+Utils.formatElapsedTime(time0, System.currentTimeMillis()));
        log.info("===");
    }

    abstract public void run(int speciesType) throws Exception;

    /**
     * return set of map keys for given species
     * @param speciesTypeKey species type key
     * @return set of map keys for given species or null if there are no map keys defined
     */
    public Set<Integer> getMapKeysForSpecies(int speciesTypeKey) {

        String speciesName = SpeciesType.getCommonName(speciesTypeKey);

        Set<Integer> mapKeySet = new HashSet<>();
        for( Map.Entry<String,String> entry: mapKeys.entrySet() ) {
            if( Utils.stringsAreEqualIgnoreCase(speciesName, entry.getValue()) ) {
                mapKeySet.add(Integer.parseInt(entry.getKey()));
            }
        }

        return mapKeySet.isEmpty() ? null : mapKeySet;
    }

    public int qcMapData(List<MapData> mdsIncoming, List<MapData> mdsInRgd, List<MapData> mdsToBeInserted, List<MapData> mdsToBeDeleted) {

        int matchCount = 0;
        for( MapData mdIncoming: mdsIncoming ) {

            boolean wasMatch = false;
            Iterator<MapData> it = mdsInRgd.listIterator();
            while( it.hasNext() ) {
                MapData mdInRgd = it.next();
                if( mdInRgd.equalsByGenomicCoords(mdIncoming)
                 && Utils.intsAreEqual(mdInRgd.getMapsDataPositionMethodId(), mdIncoming.getMapsDataPositionMethodId())) {

                    matchCount++;
                    it.remove();
                    wasMatch = true;
                    break;
                }
            }

            if( !wasMatch ) {
                // incoming did not match rgd -- must be added to RGD
                mdsToBeInserted.add(mdIncoming);
            }
        }

        // any non-matching mds in RGD must be removed
        if( mdsInRgd.size()>0 )
            mdsToBeDeleted.addAll(mdsInRgd);
        return matchCount;
    }

    /**
     * update map data
     * @param insMapData list of MapData objects to be inserted into RGD
     * @param delMapData list of MapData objects to be removed from RGD
     * @return number of updated/inserted/deleted rows
     * @throws Exception thrown when something is wrong in DAO framework
     */
    public int updateMapData(List<MapData> insMapData, List<MapData> delMapData) throws Exception {

        // total number of records updated or inserted
        int changeCount = 0;

        // insert new maps data
        if( insMapData!=null && insMapData.size()>0 ) {

            dao.insertMapPositions(insMapData);
            changeCount += insMapData.size();
        }

        // delete map data
        if( delMapData!=null && delMapData.size()>0 ) {

            dao.deleteMapPositions(delMapData);
            changeCount += delMapData.size();
        }

        return changeCount;
    }

    public DAO getDao() {
        return dao;
    }

    public void setDao(DAO dao) {
        this.dao = dao;
    }

    public void setMapKeys(Map<String,String> mapKeys) {
        this.mapKeys = mapKeys;
    }

    public Map<String,String> getMapKeys() {
        return mapKeys;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public String getLogName() {
        return logName;
    }

    public void setLogName(String logName) {
        this.logName = logName;

        log = LogManager.getLogger(logName);
    }

    public void setSrcPipeline(String srcPipeline) {
        this.srcPipeline = srcPipeline;
    }

    public String getSrcPipeline() {
        return srcPipeline;
    }

    public Map<String, String> getParams() {
        return params;
    }

    public void setParams(Map<String, String> params) {
        this.params = params;
    }
}

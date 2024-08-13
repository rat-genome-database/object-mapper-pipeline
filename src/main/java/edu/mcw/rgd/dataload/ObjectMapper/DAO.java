package edu.mcw.rgd.dataload.ObjectMapper;

import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.dao.impl.variants.VariantDAO;
import edu.mcw.rgd.dao.spring.IntStringMapQuery;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.variants.VariantMapData;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.util.*;
import java.util.Map;

/**
 * @author mtutaj
 * @since Mar 11, 2011
 * all database related code lands here -- to be centrally managed
 */
public class DAO  {

    Logger logPos = LogManager.getLogger("position_changes");

    AssociationDAO assocDAO = new AssociationDAO();
    MapDAO mapDAO = new MapDAO();
    QTLDAO qtlDAO = assocDAO.getQtlDAO();
    SSLPDAO sslpdao = assocDAO.getSslpDAO();
    StrainDAO strainDAO = assocDAO.getStrainDAO();
    RgdVariantDAO variantDAO = new RgdVariantDAO();

    VariantDAO vdao = new VariantDAO();
    GWASCatalogDAO gdao = new GWASCatalogDAO();

    public String getConnectionInfo() {
        return assocDAO.getConnectionInfo();
    }

    /**
     * get all active strain for given species
     * @param speciesTypeKey species type key
     * @return list of strains
     * @throws Exception when something really bad happens in spring framework
     */
    public List<Strain> getActiveStrains(int speciesTypeKey) throws Exception {

        // get all strains
        List<Strain> strains = strainDAO.getActiveStrains();

        // remove all strains not belonging to our species of interest
        strains.removeIf(strain -> strain.getSpeciesTypeKey() != speciesTypeKey);

        return strains;
    }

    List<Gene> getGeneAssociationsByQTL(int qtlRgdId) throws Exception {
        return assocDAO.getGeneAssociationsByQTL(qtlRgdId);
    }

    /**
     * get rgd ids for flanking markers (sslps, genes or other genomic elements) associated with given strain
     * @param strainRgdId strain rgd id
     * @return map of region name to list of rgd ids for region markers
     * @throws Exception when something really bad happens in spring framework
     */
    public java.util.Map<String, List<Integer>> getRgdIdsForFlankingMarkers(int strainRgdId) throws Exception {

        // load strain2object associations from db
        List<Strain2MarkerAssociation> strainAssocs = assocDAO.getStrain2SslpAssociations(strainRgdId);
        strainAssocs.addAll(assocDAO.getStrain2GeneAssociations(strainRgdId));
        strainAssocs.addAll(assocDAO.getStrain2StrainAssociations(strainRgdId));

        // remove alleles from the list (many genes are alleles)
        Iterator<Strain2MarkerAssociation> its2m = strainAssocs.iterator();
        while( its2m.hasNext() ) {
            if( Utils.stringsAreEqual(Utils.NVL(its2m.next().getMarkerType(), "allele"), "allele") ) {
                its2m.remove();
            }
        }

        // group markers by regions
        java.util.Map<String, List<Integer>> regionMap = new HashMap<>();
        for( Strain2MarkerAssociation assoc: strainAssocs ) {
            String region = "REGION"+assoc.getRegionName();
            List<Integer> rgdIdsOfMarkersInRegion = regionMap.get(region);
            if( rgdIdsOfMarkersInRegion==null ) {
                rgdIdsOfMarkersInRegion = new ArrayList<>();
                regionMap.put(region, rgdIdsOfMarkersInRegion);
            }
            rgdIdsOfMarkersInRegion.add(assoc.getMarkerRgdId());
        }
        return regionMap;
    }

    /**
     * get rgd ids for mutant strains (genes) associated with given strain
     * @param strainRgdId strain rgd id
     * @return map of region name to list of rgd ids for region markers
     * @throws Exception when something really bad happens in spring framework
     */
    public java.util.Map<String, List<Integer>> getRgdIdsForMutantStrains(int strainRgdId) throws Exception {

        // load strain2object associations from db
        List<Strain2MarkerAssociation> strainAssocs = assocDAO.getStrain2GeneAssociations(strainRgdId);

        // remove non-alleles from the list (f.e. flanking markers)
        Iterator<Strain2MarkerAssociation> its2m = strainAssocs.iterator();
        while( its2m.hasNext() ) {
            Strain2MarkerAssociation assoc = its2m.next();
            if( !Utils.stringsAreEqual(Utils.NVL(assoc.getMarkerType(), "allele"), "allele") ) {
                its2m.remove();
            } else {
                // a gene allele: replace rgd id of variation with the rgd id of the parent gene
                for( Gene geneAllele: assocDAO.getGeneDAO().getGeneFromVariant(assoc.getMarkerRgdId()) ) {
                    assoc.setMarkerRgdId(geneAllele.getRgdId());
                    assoc.setMarkerSymbol(geneAllele.getSymbol());
                }
            }
        }

        // group markers by regions; a region is a gene
        java.util.Map<String, List<Integer>> regionMap = new HashMap<>();
        for( Strain2MarkerAssociation assoc: strainAssocs ) {
            String region = "REGION_"+assoc.getMarkerSymbol().toUpperCase();
            List<Integer> rgdIdsOfMarkersInRegion = regionMap.get(region);
            if( rgdIdsOfMarkersInRegion==null ) {
                rgdIdsOfMarkersInRegion = new ArrayList<>();
                regionMap.put(region, rgdIdsOfMarkersInRegion);
            }
            rgdIdsOfMarkersInRegion.add(assoc.getMarkerRgdId());
        }
        return regionMap;
    }

    public List<MapData> getAllelicVariantPositionsForStrain( int strainRgdId, int mapKey ) throws Exception {

        List<Strain2MarkerAssociation> geneAlleles = assocDAO.getStrain2GeneAssociations(strainRgdId);
        geneAlleles.removeIf(i -> !Utils.stringsAreEqual(Utils.NVL(i.getMarkerType(), "allele"), "allele"));
        if( geneAlleles.isEmpty() ) {
            return Collections.emptyList();
        }

        // are any of these gene alleles have phenotypic allele variants?
        List<MapData> mds = new ArrayList<>();
        for( Strain2MarkerAssociation a: geneAlleles ) {

            List<RgdVariant> variants = variantDAO.getVariantsFromGeneRgdId(a.getMarkerRgdId());
            for( RgdVariant v: variants ) {
                mds.addAll( mapDAO.getMapData(v.getRgdId(), mapKey));
            }
        }
        return mds;
    }

    /**
     * get map data for given object and map key
     * @param rgdId object rgd id
     * @param mapKey map key
     * @return list of MapData objects
     * @throws Exception when something really bad happens in spring framework
     */
    public List<MapData> getMapPositions(int rgdId, int mapKey) throws Exception {
        return mapDAO.getMapData(rgdId, mapKey);
    }

    public List<MapData> getMapPositions(String rsId, int mapKey) throws Exception {
        return mapDAO.getMapData(rsId, mapKey);
    }

    /**
     * delete a number of map positions from database
     * @param mds list of MapData objects
     * @return nr of rows affected
     * @throws Exception when something really bad happens in spring framework
     */
    public int deleteMapPositions(List<MapData> mds ) throws Exception {
        for( MapData md: mds ) {
            logPos.info("DELETE|"+md.dump("|"));
        }
        return mapDAO.deleteMapData(mds);
    }

    /**
     * insert a number of map positions to database
     * @param mds list of MapData objects
     * @return nr of rows affected
     * @throws Exception when something really bad happens in spring framework
     */
    public int insertMapPositions(List<MapData> mds ) throws Exception {
        for( MapData md: mds ) {
            logPos.info("INSERT|"+md.dump("|"));
        }
        return mapDAO.insertMapData(mds);
    }

    /**
     * update a single MapData object
     * @param md MapData object
     * @return number of rows affected
     * @throws Exception when unexpected error in spring framework occurs
     */
    public int updateMapData(MapData md) throws Exception{
        return mapDAO.updateMapData(md);
    }

    /**
     * return map of chromosome sizes for given map
     * @param mapKey map key
     * @return map of chromosome sizes
     * @throws Exception
     */
    java.util.Map<String,Integer> getChromosomeSizes(int mapKey) throws Exception {
        return mapDAO.getChromosomeSizes(mapKey);
    }

    /**
     * get positions of given DB_SNP for given assembly
     * @param rsId DB_SNP rsId, like 'rs8143345'
     * @param mapKey map key
     * @param rgdId rgd id of an object to receive DB_SNP positions (marker)
     * @return list of map positions for this DB_SNP on the specified assembly
     * @throws Exception when something really bad happens
     */
    public List<MapData> getDbSnpPositions(String rsId, int mapKey, int rgdId) throws Exception {
        List<IntStringMapQuery.MapPair> posList = mapDAO.getDbSnpPositions(rsId, mapKey);
        List<MapData> results = new ArrayList<>(posList.size());
        for( IntStringMapQuery.MapPair pair: posList ) {
            MapData md = new MapData();
            md.setStartPos(pair.keyValue);
            md.setStopPos(pair.keyValue);
            md.setChromosome(pair.stringValue);
            md.setRgdId(rgdId);
            md.setMapKey(mapKey);
            results.add(md);
        }
        return results;
    }

    /**
     * get active qtls for given species
     * @param speciesType species type
     * @return list of active qtls for given species
     * @throws Exception if something wrong happens in spring framework
     */
    public List<QTL> getActiveQtls(int speciesType) throws Exception {
        return qtlDAO.getActiveQTLs(speciesType);
    }

    /**
     * get all active SSLP objects for DB_SNP (SSLP symbol must be a valid rsId symbol)
     * @param speciesTypeKey species type key
     * @return list of SSLP objects
     * @throws Exception when unexpected error in spring framework occurs
     */
    public List<SSLP> getActiveSSLPsForDbSnp(int speciesTypeKey) throws Exception {
        return sslpdao.getActiveSSLPsForDbSnp(speciesTypeKey);
    }

    public List<MapData> createMapDataWithDbSNP(QTL qtl, String dbSnpSource, int mapKey) throws Exception{
        String sql = "select snp_name,chromosome, position, allele, ref_allele from db_snp where map_key=? and source=? and snp_name=?";
        ArrayList<MapData> maps = new ArrayList<>();

        int rgdId = 0;
        try (Connection c = mapDAO.getConnection()){

            GWASCatalog g = gdao.getGwasCatalogByQTLRgdId(qtl.getRgdId());
            if (g==null || Utils.isStringEmpty(g.getStrongSnpRiskallele()))
                return null;

            PreparedStatement ps = c.prepareStatement(sql);
            ps.setInt(1,mapKey);
            ps.setString(2,dbSnpSource);
            ps.setString(3,qtl.getPeakRsId());
            ResultSet rs = ps.executeQuery();
            List<VariantMapData> vars = vdao.getAllVariantUnionByRsId(qtl.getPeakRsId());
            if (vars.isEmpty())
                return null;
            else
            {
                for (VariantMapData vmd : vars){
                    if (vmd.getVariantNucleotide().equals(g.getStrongSnpRiskallele()))
                        rgdId = (int) vmd.getId();
                }
            }
            if (rgdId==0)
                return null;
            while (rs.next()){
                MapData md = new MapData();
                md.setRgdId(rgdId);
                md.setChromosome(rs.getString("CHROMOSOME"));
                md.setStartPos(rs.getInt("POSITION"));
                if (rs.getString("ALLELE")==null)
                    continue;
                String allele = rs.getString("ALLELE").trim();
                md.setStopPos(md.getStartPos()+allele.length());
                md.setNotes(rs.getString("SNP_NAME"));
                md.setMapKey(mapKey);
                md.setSrcPipeline("RGD_MAPPER_PIPELINE");
                maps.add(md);
            }
            Set<MapData> set = new HashSet<>(maps);
            maps.clear();
            maps.addAll(set);
        }

        if (maps.isEmpty())
            return null;
        else
            insertMapPositions(maps);
        return maps;
    }
}

package edu.mcw.rgd.dataload.ObjectMapper;

import edu.mcw.rgd.dao.impl.GWASCatalogDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.Map;

/**
 * @author mtutaj
 * @since Oct 15, 2010
 * analyze all qtls and recalculate their positions using estimated qtl size, when possible
 */
public class QtlMapper extends BaseMapper {
    Logger logPosChanges = LogManager.getLogger("position_changes");
    BufferedWriter outOfRegionGenesFile = null;

    Map<String, Integer> chromosomeSizes;
    private Map<String,String> qtlSizeEstimate;

    // summary variables
    private int mouseQtlPipelinePositions;
    private int nonPositionableQtls;
    private int obsoletePositions;
    private int insertedPositions;
    private int matchedPositions;
    private int updatedPositions;
    private int consensusChrPositions;
    private int outOfRegionGenes;
    private String outOfRegionGenesFileName;
    private Map<Integer,String> dbSnpSource;

    public void run(int speciesType) throws Exception {

        // check if we should run out-of-region-genes reporter
        if( getParams().containsKey("reportOutOfRegionGenes") ) {
            outOfRegionGenesFile = createFileForOutOfRegionGenes();
            outOfRegionGenesFile.write("=====\nGENERATION STARTED AT "+new Date()+"\n=====\n");
        }

        String qtlSizeEstimate = getQtlSizeEstimate().get(SpeciesType.getCommonName(speciesType).toLowerCase());
        int estimatedQtlSize = Integer.parseInt(qtlSizeEstimate);

        Set<Integer> mapKeys = getMapKeysForSpecies(speciesType);
        if( mapKeys != null ) {
            // run processing for every map key specified in property file
            for( Integer mapKey: mapKeys ) {
                generateHeaderForOutOfRegionGenesFile(mapKey);
                run(speciesType, mapKey, estimatedQtlSize);
            }
        }

        if( outOfRegionGenesFile!=null ) {
            outOfRegionGenesFile.write("=====\nGENERATION FINISHED AT "+new Date()+"\n=====\n");
            outOfRegionGenesFile.close();
        }
    }

    public BufferedWriter createFileForOutOfRegionGenes() throws IOException {
        // create file name for out-of-region-genes
        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd");
        String fileName = getOutOfRegionGenesFileName().replace("{DATE}", sdt.format(new Date()));
        return new BufferedWriter(new FileWriter(fileName, true));
    }

    public void generateHeaderForOutOfRegionGenesFile(int mapKey) throws Exception {
        if( outOfRegionGenesFile!=null ) {
            outOfRegionGenesFile.write("\n");
            outOfRegionGenesFile.write("======\n");

            edu.mcw.rgd.datamodel.Map map = MapManager.getInstance().getMap(mapKey);
            int speciesTypeKey = map.getSpeciesTypeKey();
            outOfRegionGenesFile.write("SPECIES: "+SpeciesType.getCommonName(speciesTypeKey));
            outOfRegionGenesFile.write(" ["+SpeciesType.getTaxonomicName(speciesTypeKey)+"]");
            outOfRegionGenesFile.write("  "+map.getName()+" [MAPKEY="+mapKey+"]\n");
            outOfRegionGenesFile.write("======\n");
        }
    }
    /**
     * examine all qtls and calculate their sizes
     * @param speciesType species type key
     * @param mapKey map key
     * @throws Exception
     */
    public void run(int speciesType, int mapKey, int qtlSizeEstimate) throws Exception {

        mouseQtlPipelinePositions = 0;
        nonPositionableQtls = 0;
        obsoletePositions = 0;
        insertedPositions = 0;
        matchedPositions = 0;
        updatedPositions = 0;
        consensusChrPositions = 0;
        outOfRegionGenes = 0;

        log.info("Starting "+getVersion()+" for species_type="+speciesType+" and map_key="+mapKey);

        // run qtl mapper for old assembly
        this.chromosomeSizes = dao.getChromosomeSizes(mapKey);
        List<QTL> qtls = dao.getActiveQtls(speciesType);
        // process active qtls for given species
        for( QTL qtl:  qtls) {
//            boolean isGwas = false;
            int estimateSize = qtlSizeEstimate;
            // create a new qtl record object
            QtlRecord rec = new QtlRecord();
            rec.qtl = qtl;
            if (qtl.getSymbol().startsWith("GWAS")){
                estimateSize = 1;
                rec.isGwas = true;
            }
            // load genomic positions for qtl and markers
            loadPositions(rec, mapKey);
            if( rec.consensusChr!=null )
                consensusChrPositions++;

            // can the qtl be positioned at all ?
            if( !rec.isQtlSomehowPositionable() ) {
                nonPositionableQtls++;

                if( rec.qtlHasPosition() ) {
                    // qtl has position in RGD, but it is not positionable
                    //
                    // warn about this problem if the position in RGD did not come through the MOUSEQTL pipeline
                    if( !Utils.stringsAreEqual(rec.posQtl.getSrcPipeline(), "MOUSEQTL") ) {
                        log.warn("QTL "+qtl.getSymbol()+", rgd_id="+qtl.getRgdId()+" cannot be positioned but it has a position in RGD.");
                        obsoletePositions++;
                    } else {
                        mouseQtlPipelinePositions++;
                    }
                }
                else {
                    log.debug("QTL "+qtl.getSymbol()+", rgd_id="+qtl.getRgdId()+" cannot be positioned.");
                }
                continue;
            }

            if( rec.qtlHasNoPosition() ) {
                insertQtl(rec, estimateSize, mapKey);
                insertedPositions++;
                continue;
            }

            if( rec.isQtlPositionUpToDate(estimateSize, mapKey) ) {
                log.debug("QTL "+qtl.getSymbol()+", rgd_id="+qtl.getRgdId()+" is up to date. Position information not changed");
                matchedPositions++;
            } else {
                updateQtl(rec, qtlSizeEstimate, mapKey);
                updatedPositions++;
            }

            if( rec.hasOutOfRegionGenes(mapKey) ) {
                dumpOutOfRegionGenes(mapKey, rec);
                outOfRegionGenes++;
            }
        }

        printSummary(mapKey);
    }

    void printSummary(int mapKey) {
        String summary = "SUMMARY FOR MAP_KEY="+mapKey+":\n";
        if( nonPositionableQtls>0 )
            summary += "Non-positionable qtls: "+nonPositionableQtls+"\n";
        if( mouseQtlPipelinePositions>0 )
            summary += " --with positions provided by MouseQtl pipeline: "+mouseQtlPipelinePositions+"\n";
        if( obsoletePositions>0 )
            summary += " --with obsolete positions:"+obsoletePositions+"\n";
        if( insertedPositions>0 )
            summary += "Qtls with inserted positions:"+insertedPositions+"\n";
        if( updatedPositions>0 )
            summary += "Qtls with updated positions:"+updatedPositions+"\n";
        if( matchedPositions>0 )
            summary+= "Qtls with matching (unchanged) positions:"+matchedPositions+"\n";
        if( consensusChrPositions>0 )
            summary+= "Qtls with consensus chromosome positions:"+consensusChrPositions+"\n";
        if( outOfRegionGenes>0 )
            summary += "Qtls with out-of-region genes:"+outOfRegionGenes+"\n";

        log.info(summary);
    }

    void updateQtl(QtlRecord rec, int qtlSizeEstimate, int mapKey) throws Exception {

         MapData md = rec.calculatePosition(qtlSizeEstimate, mapKey);
         md = rec.boundPosition(md); // make sure qtl fits within chromosome bounds

         if( !md.equalsByGenomicCoords(rec.posQtl) ) {
             md.setKey(rec.posQtl.getKey());
             if( rec.consensusChr!=null ) {
                 md.setNotes("qtl positioned by markers consensus chromosome");
             }
             md.setSrcPipeline(getSrcPipeline());

             log.debug("UPDATE_OLD_POS "+rec.posQtl.dump("|"));
             log.debug("UPDATE_NEW_POS "+md.dump("|"));

             // separate log file with human readable position changes
             String msg = "Position update for QTL "+rec.qtl.getSymbol()+", RGD_ID="+rec.qtl.getRgdId()+", MAP_KEY="+mapKey+"\n";
             msg += "  OLD POS: "+getHumanReadablePosition(rec.posQtl)+"\n";
             msg += "  NEW POS: "+getHumanReadablePosition(md)+"\n";
             logPosChanges.info(msg);

             dao.updateMapData(md);
         }
     }

    String getHumanReadablePosition(MapData md) {
        String msg = "chr"+md.getChromosome()+":"+ Utils.formatThousands(md.getStartPos())+".."+Utils.formatThousands(md.getStopPos());
        if( md.getSrcPipeline()!=null )
            msg += "  source="+md.getSrcPipeline();
        if( md.getMapsDataPositionMethodId()!=null )
            msg += "  map_pos_method="+md.getMapPositionMethod();
        if( md.getNotes()!=null )
            msg += "  notes=("+md.getNotes()+")";
        return msg;
    }

    void insertQtl(QtlRecord rec, int estimatedQtlSize, int mapKey) throws Exception {
        if( rec.qtlHasPosition() ) {
            throw new Exception("Attempting to insert a position for an already mapped QTL, rgd_id = "+rec.qtl.getRgdId()+", map_key = "+mapKey);
        }

        MapData md = rec.calculatePosition(estimatedQtlSize, mapKey);
        md = rec.boundPosition(md); // make sure qtl fits within chromosome bounds
        if( rec.consensusChr!=null ) {
            md.setNotes("qtl positioned by markers consensus chromosome");
        }
        md.setSrcPipeline(getSrcPipeline());

        List<MapData> mds = new ArrayList<>(1);
        mds.add(md);

        log.debug("INSERT "+md.dump("|"));

        // separate log file with human readable position changes
        String msg = "Position insert for QTL "+rec.qtl.getSymbol()+", RGD_ID="+rec.qtl.getRgdId()+", MAP_KEY="+mapKey+"\n";
        msg += "  NEW POS: "+getHumanReadablePosition(md)+"\n";
        logPosChanges.info(msg);

        dao.insertMapPositions(mds);
    }

     // load positions for qtl and all markers
     public void loadPositions(QtlRecord rec, int mapKey) throws Exception {

         // qtl position
         List<MapData> mds = dao.getMapPositions(rec.qtl.getRgdId(), mapKey);
         List<MapData> mdCombined = combinePositions(mds);
         rec.posQtl = mdCombined==null ? null : mdCombined.get(0);
         if( mds.size()>1 )
             log.warn("Multiple positions for qtl "+rec.qtl.getRgdId());

         // marker positions
         List<MapData> mdsFlank1 = null;
         List<MapData> mdsFlank2 = null;
         List<MapData> mdsPeak = null;
         List<MapData> mdsMarkers = new ArrayList<MapData>(); // positions of all markers
         boolean hasRGDIds = false;
        if (rec.qtl.getFlank1RgdId()!=null || rec.qtl.getFlank2RgdId()!=null || rec.qtl.getPeakRgdId()!=null)
            hasRGDIds = true;

         // flanking region 1
         if( rec.qtl.getFlank1RgdId()!=null ) {
             mds = dao.getMapPositions(rec.qtl.getFlank1RgdId(), mapKey);
             mdsFlank1 =  combinePositions(mds);
             if( mds.size()>1 )
                 log.debug("Multiple positions for flanking region 1 "+rec.qtl.getFlank1RgdId());
             if( mdsFlank1!=null ) {
                 if( mdsFlank1.size()==1 )
                    rec.posFlank1 = mdsFlank1.get(0);
                 mdsMarkers.addAll(mdsFlank1);
             }
         }

         // flanking region 2
         if( rec.qtl.getFlank2RgdId()!=null ) {
             mds = dao.getMapPositions(rec.qtl.getFlank2RgdId(), mapKey);
             mdsFlank2 = combinePositions(mds);
             if( mds.size()>1 )
                 log.debug("Multiple positions for flanking region 2 "+rec.qtl.getFlank2RgdId());
             if( mdsFlank2!=null ) {
                 if( mdsFlank2.size()==1 )
                    rec.posFlank2 = mdsFlank2.get(0);
                 mdsMarkers.addAll(mdsFlank2);
             }
         }

         // peak region
         if( rec.qtl.getPeakRgdId()!=null ) {
             mds = dao.getMapPositions(rec.qtl.getPeakRgdId(), mapKey);
             mdsPeak = combinePositions(mds);
             if( mds.size()>1 )
                 log.debug("Multiple positions for peak marker "+rec.qtl.getPeakRgdId());
             if( mdsPeak!=null ) {
                 if( mdsPeak.size()==1 )
                    rec.posPeak = mdsPeak.get(0);
                 mdsMarkers.addAll(mdsPeak);
             }
         }

         // find/check pos in db_snp table
        if (rec.qtl.getPeakRsId()!=null && rec.isGwas && !hasRGDIds){
            mds = dao.getMapPositions(rec.qtl.getPeakRsId(),mapKey);
            mdsPeak = combinePositions(mds);
            if( mdsPeak!=null ) {
                if( mdsPeak.size()==1 )
                    rec.posPeak = mdsPeak.get(0);
                mdsMarkers.addAll(mdsPeak);
            }
            else {
                // create marker map data for qtl map data
                // grab from dbSnp
                List<MapData> peakRsMaps = null;
                if (mapKey == 38) {
                    peakRsMaps = dao.createMapDataWithDbSNP(rec.qtl, dbSnpSource.get(mapKey), mapKey);
                }
                else { // or grab from variant tables
                    peakRsMaps = dao.createMapDataWOdbSnp(rec.qtl, mapKey);
                }

                if (peakRsMaps != null) {
                    if (peakRsMaps.size() == 1)
                        rec.posPeak = peakRsMaps.get(0);
                    mdsMarkers.addAll(peakRsMaps);
                }
            }
        }

         // check if marker positions are same chromosomes
         if( positionsOnOneChromosome(mdsMarkers) )
             return; // all markers are positioned on one chromosome

         // determine consensus chromosome:
         // one chromosome must be most frequent; if there are ties, there is no consensus chromosome
         rec.consensusChr = getConsensusChromosome(mdsMarkers);
         if( rec.consensusChr==null ) {
             // there is no consensus chromosome
             // clear all flanking markers positions -- we cannot position this qtl
             log.debug("No consensus chromosome for qtl "+rec.qtl.getSymbol()+" "+rec.qtl.getRgdId());
             rec.posFlank1 = rec.posFlank2 = rec.posPeak = null;
             return;
         }

         // there is a consensus chromosome
         // process positions for markers on consensus chromosome
         rec.posFlank1 = getConsensusPosition(rec.consensusChr, mdsFlank1);
         rec.posFlank2 = getConsensusPosition(rec.consensusChr, mdsFlank2);
         rec.posPeak = getConsensusPosition(rec.consensusChr, mdsPeak);
     }

    /**
     *  if there multiple positions on same chromosome, pick the lowest and highest positions as marker position;
     *  return list of combined positions per chromosome
     */
    List<MapData> combinePositions( List<MapData> mds ) throws CloneNotSupportedException {

        // handle trivial cases
        if( mds==null || mds.isEmpty() )
            return null;

        // positions per chromosome
        List<MapData> mdCombined = new ArrayList<MapData>();
        for( MapData md: mds ) {
            if( md.getStartPos()>0 && md.getStopPos()>0 )
                combinePositions(md, mdCombined);
        }
        return mdCombined.isEmpty() ? null : mdCombined;
    }

    void combinePositions(MapData md, List<MapData> mdCombined) {

        // check if 'md' chromosome is on 'mdCombined' list, if not, add it
        if( mdCombined.isEmpty() ) {
            mdCombined.add(md);
            return;
        }
        MapData mdExisting = null; // existing combined position on same chr as 'md'
        for( MapData md2: mdCombined ) {
            if( md2.getChromosome().equals(md.getChromosome()) ) {
                mdExisting = md2;
                break;
            }
        }
        if( mdExisting==null ) {
            mdCombined.add(md);
            return;
        }

        // 'md' chr is same as existing combined positions on the same chromosome:
        // merge position data
        if( md.getStartPos() < mdExisting.getStartPos() )
            mdExisting.setStartPos(md.getStartPos()); // lower combined start pos
        if( md.getStopPos() > mdExisting.getStopPos() )
            mdExisting.setStopPos(md.getStopPos());  // higher combined end pos
    }

    boolean positionsOnOneChromosome(List<MapData> mds) {

        if( mds.size()<=1 )
            return true;
        String chr = mds.get(0).getChromosome();
        for( int i=1; i<mds.size(); i++ ) {
            if( !chr.equals(mds.get(i).getChromosome()) )
                return false;
        }
        return true;
    }

    String getConsensusChromosome(List<MapData> mds) {

        // assertion
        if( mds==null || mds.size()<=1 )
            return null;

        // build map of chromosome frequencies
        Map<String, Integer> chrFreq = new HashMap<String, Integer>();
        for( MapData md: mds ) {
            String chr = md.getChromosome();
            Integer freq = chrFreq.get(chr);
            if( freq==null )
                freq = 1;
            else
                freq++;
            chrFreq.put(chr, freq);
        }

        // extract most frequent chromosome
        String chr1 = "";
        Integer freq1 = 0;
        for( Map.Entry<String, Integer> entry: chrFreq.entrySet() ) {
            if( entry.getValue()>freq1 ) {
                // new most frequent chromosome found
                chr1 = entry.getKey();
                freq1 = entry.getValue();
            }
        }
        // remove most frequent chromosome from freq map
        chrFreq.remove(chr1);

        // extract second most frequent chromosome
        //String chr2 = "";
        Integer freq2 = 0;
        for( Map.Entry<String, Integer> entry: chrFreq.entrySet() ) {
            if( entry.getValue()>freq2 ) {
                // new most frequent chromosome found
                //chr2 = entry.getKey();
                freq2 = entry.getValue();
            }
        }

        // if both first and second most frequent chromosomes have the same frequency,
        // consensus chromosome cannot be determined
        if( freq1.equals(freq2) )
            return null;

        // we found consensus chromosome!
        return chr1;
    }

    MapData getConsensusPosition(String consensusChr, List<MapData>mds) {

        if( mds==null )
            return null;

        for( MapData md: mds ) {
            if( md.getChromosome().equals(consensusChr) ) {
                return md;
            }
        }
        return null;
    }

    public void setQtlSizeEstimate(Map<String,String> qtlSizeEstimate) {
        this.qtlSizeEstimate = qtlSizeEstimate;
    }

    public Map<String,String> getQtlSizeEstimate() {
        return qtlSizeEstimate;
    }

    void dumpOutOfRegionGenes(int mapKey, QtlRecord rec ) throws Exception {

        if( outOfRegionGenesFile==null ) {
            return;
        }

        String geneInfo = "";
        for( Gene gene: rec.outOfRegionGenes ) {
            for( MapData posGene: dao.getMapPositions(gene.getRgdId(), mapKey) ) {
                // skip notes and source for genes
                posGene.setNotes(null);
                posGene.setSrcPipeline(null);
                geneInfo +="   GENE "+gene.getSymbol()+" RGDID:"+gene.getRgdId()+" POS "+getHumanReadablePosition(posGene)+"\n";
            }
        }

        if( !geneInfo.isEmpty() ) {
            outOfRegionGenesFile.write("QTL " + rec.qtl.getSymbol() + " RGDID:" + rec.qtl.getRgdId());
            outOfRegionGenesFile.write(" POS " + getHumanReadablePosition(rec.posQtl) + "\n");
            outOfRegionGenesFile.write(geneInfo);
        }
    }

    public void setOutOfRegionGenesFileName(String outOfRegionGenesFileName) {
        this.outOfRegionGenesFileName = outOfRegionGenesFileName;
    }

    public String getOutOfRegionGenesFileName() {
        return outOfRegionGenesFileName;
    }

    public void setDbSnpSource(Map<Integer,String> dbSnpSource) {
        this.dbSnpSource = dbSnpSource;
    }

    public Map<Integer,String> getDbSnpSource() {
        return dbSnpSource;
    }


    class QtlRecord {
        public QTL qtl;
        public MapData posQtl;    // genomic position for qtl
        public MapData posFlank1; // genomic position for flank1 marker
        public MapData posFlank2; // genomic position for flank2 marker
        public MapData posPeak;   // genomic position for peak marker
        public String consensusChr;
        public List<Gene> outOfRegionGenes;
        public boolean isGwas = false;

        // qtl has position on a map
        boolean qtlHasPosition() {
            return posQtl!=null;
        }

        // qtl does not have position on a map
        boolean qtlHasNoPosition() {
            return posQtl==null;
        }

        boolean hasPositionableMarker() {
            return posFlank1!=null || posFlank2!=null || posPeak!=null;
        }

        boolean hasSameMarkerChromosomes() {
            Set<String> chrSet = new HashSet<>();
            if( posFlank1!=null ) {
                chrSet.add(posFlank1.getChromosome());
            }
            if( posFlank2!=null ) {
                chrSet.add(posFlank2.getChromosome());
            }
            if( posPeak!=null ) {
                chrSet.add(posPeak.getChromosome());
            }
            return chrSet.size()==1;
        }

        boolean qtlHasFlankingMarkerPosition() {
            return posFlank1!=null && posFlank2!=null;
        }

        MapData calculatePosition(int estimatedQtlSize, int mapKey) throws Exception {

             MapData md = new MapData();
             md.setRgdId(qtl.getRgdId());
             md.setMapKey(mapKey);

             // *NEW* check if the QTL is a GWAS QTL because it has a peak rs Id
            if (isGwas)
            {
                md.setMapsDataPositionMethodId(3); // POS_BY_PEAK_ONLY
//                if (calculateFromPeakRsId(md, qtl.getRgdId()))
                md.setStartPos(posPeak.getStartPos());
                md.setStopPos(posPeak.getStopPos());
                md.setChromosome(posPeak.getChromosome());
                return md;
            }

             // first try qtl with both flanking markers
             if( qtlHasFlankingMarkerPosition() ) {
                 md.setMapsDataPositionMethodId(1); // POS_BY_FLANKS
                 md.setChromosome(posFlank1.getChromosome());
                 md.setStartPos(posFlank1.getStartPos());
                 md.setStopPos(posFlank2.getStopPos());
                 return md;
             }

             // Do we have a peak position ?
             if( posPeak!=null ) {

                 // do we have one peak positionable?
                 MapData posFlank = posFlank1!=null ? posFlank1 : posFlank2;

                 // it is possible for the flank and peak marker to be the same
                 // in this case we probably shouldn't the flank as it makes the QTL
                 // appear very small
                 if( posFlank!=null && posFlank.getRgdId()!=posPeak.getRgdId() ) {
                     // calculate on flank and peak
                     md.setMapsDataPositionMethodId(2); // POS_BY_SINGLE_FLANK_AND_PEAK
                     calculateFromFlankAndPeak(md, posFlank, posPeak);
                     return md;
                 }
                 else {
                     // else peak alone with adjusted size
                     md.setMapsDataPositionMethodId(5); // POS_BY_PEAK_WITH_ADJ_SIZE
                     calculateFromPeakOnly(md, estimatedQtlSize);
                     return md;
                 }
             }

             // No peak position
             // then try 1 flank alone
             if( posFlank1!=null || posFlank2!=null ) {
                 md.setMapsDataPositionMethodId(4); // POS_BY_SINGLE_FLANK_ONLY
                 calculateFromFlankOnly(md, estimatedQtlSize);
                 return md;
             }

             // if we get here we have a logic error
             // couldn't use any fo the positioning methods
             throw new Exception("Couldn't calculate position for QTL rgd_id="+qtl.getRgdId());
        }

        boolean calculateFromPeakRsId(MapData md, int qtlRgdId) throws Exception {
            GWASCatalogDAO gdao = new GWASCatalogDAO();
            GWASCatalog g = gdao.getGwasCatalogByQTLRgdId(qtlRgdId);
            if (g == null)
                return false;
            try {
                int start = Integer.parseInt(g.getPos());
                md.setStartPos(start);
                md.setStopPos(start + 1);
                md.setChromosome(g.getChr());
                return true;
            }
            catch (Exception e){ // gwas position is null
                return false;
            }
        }

        void calculateFromFlankAndPeak(MapData md, MapData mdFlank, MapData mdPeak) {

             md.setChromosome(mdPeak.getChromosome());

             int mid_point = mdPeak.getStartPos() + ( (mdPeak.getStopPos() - mdPeak.getStartPos()) / 2 );
             // Don't know if flank is upstream or downstream of peak
             // so just find the size of the flank to peak distance
             // and make our final size symmetrical about the peak
             int downstream = Math.min(mdFlank.getStartPos(), mid_point);
             int upstream = Math.max(mdFlank.getStopPos(), mid_point);
             int flank2peak = upstream - downstream;

             md.setStartPos(mid_point - flank2peak);
             md.setStopPos(mid_point + flank2peak);
        }

        void calculateFromPeakOnly(MapData md, int estimatedQtlSize) {

             md.setChromosome(posPeak.getChromosome());
             int mid_point = posPeak.getStartPos() + ( (posPeak.getStopPos() - posPeak.getStartPos()) / 2 );
             int half_size = estimatedQtlSize / 2;
             md.setStartPos(mid_point - half_size);
             md.setStopPos(mid_point + half_size);
        }

        void calculateFromFlankOnly(MapData md, int estimatedQtlSize) {

             if( posFlank1 != null ) {
                 md.setChromosome(posFlank1.getChromosome());
                 md.setStartPos(posFlank1.getStartPos());
                 md.setStopPos(md.getStartPos() + estimatedQtlSize);
                 return;
             }
             else if( posFlank2!=null ) {
                 md.setChromosome(posFlank2.getChromosome());
                 md.setStopPos(posFlank2.getStopPos());
                 md.setStartPos(md.getStopPos() - estimatedQtlSize);
                 return;
             }

             // if we get here logic error
             log.warn("Logic error calculating position for QTL rgd_id = "+qtl.getRgdId());
        }

        MapData boundPosition(MapData md) {

            // ensure start < stop
            if( md.getStartPos() > md.getStopPos() ) {
                log.debug("swapping start-stop positions for qtl=" + qtl.getRgdId());
                int val = md.getStartPos();
                md.setStartPos(md.getStopPos());
                md.setStopPos(val);
            }

            int startBound = 1;
            int stopBound = chromosomeSizes.get(md.getChromosome());
            if( md.getStartPos() < startBound ) {
                md.setStartPos(startBound);
                log.debug("adjusting start position for qtl="+qtl.getRgdId());
            }
            if( md.getStopPos() > stopBound ) {
                md.setStopPos(stopBound);
                log.debug("adjusting stop position for qtl="+qtl.getRgdId());
            }
            return md;
        }

        boolean isQtlSomehowPositionable() {
            // not Positionable if no positioned markers
            // or if markers are on different chromosomes
            if( !hasPositionableMarker() ) {
                log.debug("No positionable markers for "+qtl.getSymbol());
                return false;
            }

            // check if all markers are on same chromosome
            if( !hasSameMarkerChromosomes() ) {
                log.warn("Different chromosomes for markers for "+qtl.getSymbol());
                return false;
            }

            // if qtl mapping position method is 5 or 6 (set by mouseQtl pipeline only)
            // we should not try to reposition the qtl
            if( qtlHasPosition() ) {
                int mapPosMethod = posQtl.getMapsDataPositionMethodId()!=null ? posQtl.getMapsDataPositionMethodId() : 0;
                if( mapPosMethod==5 || mapPosMethod==6 ) {
                    // MouseQtl:
                    // 5: the qtl size is properly adjusted by mouse qtl pipeline
                    // 6: the qtl position is taken directly from incoming data
                    log.debug("Qtl "+qtl.getSymbol()+" won't be repositioned, because its mapPosMethod is="+mapPosMethod);
                    return false;
                }
            }
            // passed all tests so we're good
            return true;
        }

        boolean isQtlPositionUpToDate(int qtlSizeEstimate, int mapKey) throws Exception {

            if( posQtl==null )
                return false;

            int current_pos_method = posQtl.getMapsDataPositionMethodId();
            MapData md = boundPosition(calculatePosition(qtlSizeEstimate, mapKey));
            int calc_pos_method = md.getMapsDataPositionMethodId();

            if( current_pos_method==calc_pos_method && posQtl.equalsByGenomicCoords(md) )
                return true;

            if( current_pos_method > calc_pos_method ) {
                log.warn("Position information for QTL rgd_id = "+qtl.getRgdId()+
                        " appears to be downgraded. Current pos_method="+current_pos_method+
                        ", new pos_method="+calc_pos_method);
            }
            return false;
        }

        boolean hasOutOfRegionGenes(int mapKey) throws Exception {
            if( posQtl==null )
                return false;

            for( Gene gene: dao.getGeneAssociationsByQTL(qtl.getRgdId()) ) {
                boolean geneInRegion = false;
                for( MapData posGene: dao.getMapPositions(gene.getRgdId(), mapKey) ) {
                    if( posGene.getChromosome().equals(posQtl.getChromosome()) &&
                        posGene.getStopPos() >= posQtl.getStartPos() &&
                        posGene.getStartPos() <= posQtl.getStopPos() ) {
                        geneInRegion = true;
                        break;
                    }
                }
                if( !geneInRegion ) {
                    if( this.outOfRegionGenes==null )
                        this.outOfRegionGenes = new ArrayList<Gene>();
                    this.outOfRegionGenes.add(gene);
                }
            }
            return this.outOfRegionGenes!=null;
        }
    }

}
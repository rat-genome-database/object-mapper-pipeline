package edu.mcw.rgd.dataload.ObjectMapper;

import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.Utils;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;

import java.util.HashMap;
import java.util.Map;

/**
 * @author mtutaj
 * @since Mar 11, 2011
 */
public class Manager {

    DAO dao;
    private String version;

    static public void main(String[] args) throws Exception {

        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        Manager manager = (Manager) (bf.getBean("manager"));

        String usageMsg = "Please specify a module (or modules) and optionally species: [-qtls | -strains | -markers] [-rat | -mouse | -human]";
        if( args.length<1 ) {
            System.out.println(usageMsg);
            return;
        }
        String beanId = null;
        int speciesType = SpeciesType.ALL;
        Map<String,String> params = new HashMap<>();

        // parse cmdline params
        for( String arg: args ) {
            switch (arg) {
                case "-strains":
                    beanId = "strainMapper";
                    break;
                case "-qtls":
                    beanId = "qtlMapper";
                    break;
                case "-markers":
                    beanId = "markerMapper";
                    break;
                case "-rat":
                    speciesType = SpeciesType.RAT;
                    break;
                case "-mouse":
                    speciesType = SpeciesType.MOUSE;
                    break;
                case "-human":
                    speciesType = SpeciesType.HUMAN;
                    break;
                case "-reportOutOfRegionGenes":
                    params.put("reportOutOfRegionGenes", "1");
                    break;
                case "-test_email":
                    String mailServer = "localhost";
                    String[] recipients = new String[]{"mtutaj@mcw.edu"};
                    String mailFrom = "rgddata@travis.mcw.edu";
                    Utils.sendMail(mailServer, mailFrom, recipients, "TEST", "A message from travis");
                    return;
            }
        }

        BaseMapper mapper = (BaseMapper) bf.getBean(beanId);
        mapper.setDao(manager.getDao());
        mapper.setParams(params);
        try {
            mapper.log.info(manager.getVersion());
            if (speciesType == SpeciesType.ALL) {
                mapper.start(SpeciesType.RAT);
                mapper.start(SpeciesType.MOUSE);
                mapper.start(SpeciesType.HUMAN);
            } else {
                mapper.start(speciesType);
            }
        } catch(Exception e) {
            Utils.printStackTrace(e, mapper.log);
            throw new Exception(e);
        }
    }

    public DAO getDao() {
        return dao;
    }

    public void setDao(DAO dao) {
        this.dao = dao;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }
}

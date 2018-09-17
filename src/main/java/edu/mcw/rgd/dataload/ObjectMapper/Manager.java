package edu.mcw.rgd.dataload.ObjectMapper;

import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;

import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: Mar 11, 2011
 * Time: 3:20:44 PM
 */
public class Manager {

    Logger log = Logger.getLogger("detail");
    DAO dao;
    private String version;

    static public void main(String[] args) throws Exception {

        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        Manager manager = (Manager) (bf.getBean("manager"));
        manager.log.info(manager.getVersion());

        String usageMsg = "Please specify a module (or modules) and optionally species: [-qtls | -strains] [-rat | -mouse | -human]";
        if( args.length<1 ) {
            manager.log.info(usageMsg);
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
                case "-sslps":
                    beanId = "sslpMapper";
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
                    String mailFrom = "rgddata@kirwan.mcw.edu";
                    Utils.sendMail(mailServer, mailFrom, recipients, "TEST", "A message from kirwan");
                    return;
            }
        }

        BaseMapper mapper = (BaseMapper) bf.getBean(beanId);
        mapper.setDao(manager.getDao());
        mapper.setParams(params);
        if( speciesType==SpeciesType.ALL ) {
            mapper.run(SpeciesType.RAT);
            mapper.run(SpeciesType.MOUSE);
            mapper.run(SpeciesType.HUMAN);
        } else {
            mapper.run(speciesType);
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

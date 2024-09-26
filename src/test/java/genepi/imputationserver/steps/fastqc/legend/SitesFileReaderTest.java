package genepi.imputationserver.steps.fastqc.legend;

import junit.framework.TestCase;
import org.junit.Test;

import java.util.List;

public class SitesFileReaderTest extends TestCase {

    @Test
    public void testReadFile() throws Exception {
        SitesFileReader fileReader = new SitesFileReader("test-data/configs/hapmap-chr1/ref-panels/hapmap_r22.chr1.CEU.hg19_impute.legend.gz", "eur");
        //List<SitesEntry> sites = fileReader.findByPosition("1", 1060174);
        //assertEquals("rs7548798", sites.get(sites.size() - 1).getRsId());
        fileReader.close();
    }

}
package genepi.imputationserver.steps.fastqc.legend;

import junit.framework.TestCase;
import org.junit.Test;

public class SitesFileReaderTest extends TestCase {

    @Test
    public void testReadFile() throws Exception {
        SitesFileReader fileReader = new SitesFileReader("test-data/configs/hapmap-chr1/ref-panels/hapmap_r22.chr1.CEU.hg19_impute.legend.gz", "eur");
        assertEquals("rs7548798", fileReader.findByPosition("1", 1060174).getRsId());
        fileReader.close();
    }

}
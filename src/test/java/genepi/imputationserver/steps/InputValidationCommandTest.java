package genepi.imputationserver.steps;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.junit.Test;

import genepi.imputationserver.util.AbstractTestcase;
import genepi.imputationserver.util.RefPanelUtil;
import genepi.imputationserver.util.report.CloudgeneReport;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class InputValidationCommandTest extends AbstractTestcase {

	private static final String TABIX_HOME = "files/bin/tabix";

	private static final String CLOUDGENE_LOG = "cloudgene.report.json";

	public static final boolean VERBOSE = true;

	@Test
	public void testHg19DataWithBuild38() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/three";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setBuild("hg38");
		command.setupTabix(TABIX_HOME);

		assertEquals(-1, (int) command.call());

		CloudgeneReport CloudgeneLog = new CloudgeneReport(CLOUDGENE_LOG);

		assertTrue(CloudgeneLog.hasInMemory("This is not a valid hg38 encoding."));
		assertTrue(CloudgeneLog.hasInMemory("[ERROR]"));

	}

	@Test(expected = IOException.class)
	public void testWithWrongReferencePanel() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
//		TODO: test how command handles wrong filename --> nor RefPanelUtil.
		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "missing-reference-panel");

	}

	@Test
	public void testHg38DataWithBuild19() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr20";
		String inputFolder = "test-data/data/chr20-unphased-hg38";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setBuild("hg19");
		command.setPopulation("eur");
		command.setupTabix(TABIX_HOME);

		assertEquals(-1, (int) command.call());

		CloudgeneReport CloudgeneLog = new CloudgeneReport(CLOUDGENE_LOG);
		assertTrue(CloudgeneLog.hasInMemory("This is not a valid hg19 encoding."));
		assertTrue(CloudgeneLog.hasInMemory("[ERROR]"));

	}

	@Test
	public void testWrongVcfFile() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/wrong_vcf";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setPopulation("eur");
		command.setupTabix(TABIX_HOME);

		assertEquals(-1, (int) command.call());

		// check error message
		CloudgeneReport CloudgeneLog = new CloudgeneReport(CLOUDGENE_LOG);
		assertTrue(CloudgeneLog.hasInMemory("[ERROR] Unable to parse header with error"));

	}

	@Test
	public void testMixedPopulation() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setPopulation("mixed");
		command.setupTabix(TABIX_HOME);

		assertEquals(0, (int) command.call());

	}

	@Test
	public void testCorrectHrcPopulation() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hrc-fake");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setPopulation("mixed");
		command.setupTabix(TABIX_HOME);

		assertEquals(0, (int) command.call());

	}

	@Test
	public void testWrongHrcPopulation() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hrc-fake");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setPopulation("aas");
		command.setupTabix(TABIX_HOME);

		assertEquals(-1, (int) command.call());

		// check error message
		CloudgeneReport CloudgeneLog = new CloudgeneReport(CLOUDGENE_LOG);
		assertTrue(CloudgeneLog.hasInMemory("[ERROR] Population 'aas' is not supported by reference panel 'hrc-fake'"));

	}

	@Test
	public void testWrong1KP3Population() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "phase3-fake");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setPopulation("asn");
		command.setupTabix(TABIX_HOME);

		assertEquals(-1, (int) command.call());

		// check error message
		CloudgeneReport CloudgeneLog = new CloudgeneReport(CLOUDGENE_LOG);
		assertTrue(
				CloudgeneLog.hasInMemory("[ERROR] Population 'asn' is not supported by reference panel 'phase3-fake'"));

	}

	@Test
	public void testWrongTopmedPopulation() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "TOPMedfreeze6-fake");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setPopulation("asn");
		command.setupTabix(TABIX_HOME);

		assertEquals(-1, (int) command.call());

		CloudgeneReport CloudgeneLog = new CloudgeneReport(CLOUDGENE_LOG);
		assertTrue(CloudgeneLog
				.hasInMemory("[ERROR] Population 'asn' is not supported by reference panel 'TOPMedfreeze6-fake'"));

	}

	@Test
	public void testUnorderedVcfFile() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		// input folder contains no vcf or vcf.gz files
		String inputFolder = "test-data/data/unorderd";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setPopulation("eur");
		command.setupTabix(TABIX_HOME);

		assertEquals(-1, (int) command.call());

		CloudgeneReport CloudgeneLog = new CloudgeneReport(CLOUDGENE_LOG);
		assertTrue(CloudgeneLog.hasInMemory("[ERROR] The provided VCF file is malformed"));
		assertTrue(CloudgeneLog.hasInMemory("Error during index creation"));

	}

	@Test
	public void testWrongChromosomes() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/wrong_chrs";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setPopulation("eur");
		command.setupTabix(TABIX_HOME);

		assertEquals(-1, (int) command.call());

		CloudgeneReport CloudgeneLog = new CloudgeneReport(CLOUDGENE_LOG);
		assertTrue(CloudgeneLog.hasInMemory("[ERROR] The provided VCF file contains more than one chromosome."));

	}

	@Test
	public void testSingleUnphasedVcfWithEagle() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setPopulation("eur");
		command.setupTabix(TABIX_HOME);

		assertEquals(0, (int) command.call());

		CloudgeneReport CloudgeneLog = new CloudgeneReport(CLOUDGENE_LOG);
		assertTrue(CloudgeneLog.hasInMemory("[OK] 1 valid VCF file(s) found."));
		assertTrue(CloudgeneLog.hasInMemory("Samples: 41"));
		assertTrue(CloudgeneLog.hasInMemory("Chromosomes: 1"));
		assertTrue(CloudgeneLog.hasInMemory("SNPs: 905"));
		assertTrue(CloudgeneLog.hasInMemory("Chunks: 1"));
		assertTrue(CloudgeneLog.hasInMemory("Datatype: unphased"));
		assertTrue(CloudgeneLog.hasInMemory("Reference Panel: hapmap2"));
		assertTrue(CloudgeneLog.hasInMemory("Phasing: eagle"));

	}

	@Test
	public void testThreeUnphasedVcfWithEagle() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/three";
		// create workflow context
		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setPopulation("eur");
		command.setupTabix(TABIX_HOME);

		assertEquals(0, (int) command.call());

		CloudgeneReport CloudgeneLog = new CloudgeneReport(CLOUDGENE_LOG);
		assertTrue(CloudgeneLog.hasInMemory("[OK] 3 valid VCF file(s) found."));
		assertTrue(CloudgeneLog.hasInMemory("Samples: 41"));
		assertTrue(CloudgeneLog.hasInMemory("Chromosomes: 2 3 4"));
		assertTrue(CloudgeneLog.hasInMemory("SNPs: 2715"));
		assertTrue(CloudgeneLog.hasInMemory("Chunks: 3"));
		assertTrue(CloudgeneLog.hasInMemory("Datatype: unphased"));
		assertTrue(CloudgeneLog.hasInMemory("Reference Panel: hapmap2"));
		assertTrue(CloudgeneLog.hasInMemory("Phasing: eagle"));

	}

	@Test
	public void testTabixIndexCreationChr20() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		// input folder contains no vcf or vcf.gz files
		String inputFolder = "test-data/data/chr20-phased";

		// create workflow context
		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setPopulation("eur");
		command.setupTabix(TABIX_HOME);

		assertEquals(0, (int) command.call());

		CloudgeneReport CloudgeneLog = new CloudgeneReport(CLOUDGENE_LOG);
		assertTrue(CloudgeneLog.hasInMemory("[OK] 1 valid VCF file(s) found."));

		// test tabix index and count snps
		String vcfFilename = inputFolder + "/chr20.R50.merged.1.330k.recode.small.vcf.gz";
		VCFFileReader vcfReader = new VCFFileReader(new File(vcfFilename),
				new File(vcfFilename + TabixUtils.STANDARD_INDEX_EXTENSION), true);
		CloseableIterator<VariantContext> snps = vcfReader.query("20", 1, 1000000000);
		int count = 0;
		while (snps.hasNext()) {
			snps.next();
			count++;
		}
		snps.close();
		vcfReader.close();

		// check snps
		assertEquals(7824, count);

	}

	@Test
	public void testTabixIndexCreationChr1() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		// input folder contains no vcf or vcf.gz files
		String inputFolder = "test-data/data/single";

		// create workflow context
		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(getFiles(inputFolder));
		command.setReference(panel);
		command.setPopulation("eur");
		command.setupTabix(TABIX_HOME);

		assertEquals(0, (int) command.call());

		CloudgeneReport CloudgeneLog = new CloudgeneReport(CLOUDGENE_LOG);
		assertTrue(CloudgeneLog.hasInMemory("[OK] 1 valid VCF file(s) found."));
		// test tabix index and count snps
		String vcfFilename = inputFolder + "/minimac_test.50.vcf.gz";
		VCFFileReader vcfReader = new VCFFileReader(new File(vcfFilename),
				new File(vcfFilename + TabixUtils.STANDARD_INDEX_EXTENSION), true);
		CloseableIterator<VariantContext> snps = vcfReader.query("1", 1, 1000000000);
		int count = 0;
		while (snps.hasNext()) {
			snps.next();
			count++;
		}
		snps.close();
		vcfReader.close();

		// check snps
		assertEquals(905, count);

	}

}

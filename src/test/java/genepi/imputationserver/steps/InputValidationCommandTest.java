package genepi.imputationserver.steps;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import genepi.imputationserver.util.AbstractTestcase;
import genepi.imputationserver.util.RefPanel;
import genepi.imputationserver.util.OutputReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class InputValidationCommandTest extends AbstractTestcase {

	private static final String CLOUDGENE_LOG = "cloudgene.report";

	@Test
	public void testHg19DataWithBuild38() throws Exception {

		String inputFolder = "test-data/data/three";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr1/hapmap2.json");
		command.setBuild("hg38");

		assertEquals(1, (int) command.call());

		OutputReader log = new OutputReader(CLOUDGENE_LOG);
		assertTrue(log.hasInMemory("This is not a valid hg38 encoding."));
		assertTrue(log.hasInMemory("::error::"));

	}

	@Test
	public void testHg38DataWithBuild19() throws Exception {

		String inputFolder = "test-data/data/chr20-unphased-hg38";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr20/hapmap2.json");
		command.setBuild("hg19");
		command.setPopulation("eur");

		assertEquals(1, (int) command.call());

		OutputReader log = new OutputReader(CLOUDGENE_LOG);
		log.view();
		assertTrue(log.hasInMemory("This is not a valid hg19 encoding."));
		assertTrue(log.hasInMemory("::error::"));

	}

	@Test
	public void testWrongVcfFile() throws Exception {

		String inputFolder = "test-data/data/wrong_vcf";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr1/hapmap2.json");
		command.setPopulation("eur");

		assertEquals(1, (int) command.call());

		// check error message
		OutputReader log = new OutputReader(CLOUDGENE_LOG);
		assertTrue(log.hasInMemory("::error:: Unable to parse header with error"));

	}

	@Test
	public void testMixedPopulation() throws Exception {

		String inputFolder = "test-data/data/single";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr1/hapmap2.json");
		command.setPopulation("mixed");

		assertEquals(0, (int) command.call());

	}

	@Test
	public void testCorrectHrcPopulation() throws Exception {

		String inputFolder = "test-data/data/single";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr1/hrc-fake.json");
		command.setPopulation("mixed");

		assertEquals(0, (int) command.call());

	}

	@Test
	public void testWrongHrcPopulation() throws Exception {

		String inputFolder = "test-data/data/single";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr1/hrc-fake.json");
		command.setPopulation("aas");

		assertEquals(1, (int) command.call());

		// check error message
		OutputReader log = new OutputReader(CLOUDGENE_LOG);
		assertTrue(log.hasInMemory("::group type=error::"));
		assertTrue(log.hasInMemory("Population 'aas' is not supported by reference panel 'hrc-fake'"));

	}

	@Test
	public void testWrong1KP3Population() throws Exception {

		String inputFolder = "test-data/data/single";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr1/phase3-fake.json");
		command.setPopulation("asn");

		assertEquals(1, (int) command.call());

		// check error message
		OutputReader log = new OutputReader(CLOUDGENE_LOG);
		assertTrue(log.hasInMemory("::group type=error::"));
		assertTrue(log.hasInMemory("Population 'asn' is not supported by reference panel 'phase3-fake'"));

	}

	@Test
	public void testWrongTopmedPopulation() throws Exception {

		String inputFolder = "test-data/data/single";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr1/TOPMedfreeze6-fake.json");
		command.setPopulation("asn");

		assertEquals(1, (int) command.call());

		OutputReader log = new OutputReader(CLOUDGENE_LOG);
		assertTrue(log.hasInMemory("::group type=error::"));
		assertTrue(log.hasInMemory("Population 'asn' is not supported by reference panel 'TOPMedfreeze6-fake'"));
	}

	@Test
	public void testUnorderedVcfFile() throws Exception {

		String inputFolder = "test-data/data/unorderd";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr1/hapmap2.json");
		command.setPopulation("eur");

		assertEquals(1, (int) command.call());

		OutputReader log = new OutputReader(CLOUDGENE_LOG);
		log.view();
		assertTrue(log.hasInMemory("::error:: The provided VCF file is malformed"));
		assertTrue(log.hasInMemory("Error during index creation"));

	}

	@Test
	public void testWrongChromosomes() throws Exception {

		String inputFolder = "test-data/data/wrong_chrs";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr1/hapmap2.json");
		command.setPopulation("eur");

		assertEquals(1, (int) command.call());

		OutputReader log = new OutputReader(CLOUDGENE_LOG);
		assertTrue(log.hasInMemory("::error:: The provided VCF file contains more than one chromosome."));

	}

	@Test
	public void testSingleUnphasedVcfWithEagle() throws Exception {

		String inputFolder = "test-data/data/single";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr1/hapmap2.json");
		command.setPopulation("eur");

		assertEquals(0, (int) command.call());

		OutputReader log = new OutputReader(CLOUDGENE_LOG);
		assertTrue(log.hasInMemory("1 valid VCF file(s) found."));
		assertTrue(log.hasInMemory("Samples: 41"));
		assertTrue(log.hasInMemory("Chromosomes: 1"));
		assertTrue(log.hasInMemory("SNPs: 905"));
		assertTrue(log.hasInMemory("Chunks: 1"));
		assertTrue(log.hasInMemory("Datatype: unphased"));
		assertTrue(log.hasInMemory("Reference Panel: hapmap2"));
		assertTrue(log.hasInMemory("Phasing: eagle"));

	}

	@Test
	public void testThreeUnphasedVcfWithEagle() throws Exception {

		String inputFolder = "test-data/data/three";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr1/hapmap2.json");
		command.setPopulation("eur");

		assertEquals(0, (int) command.call());

		OutputReader log = new OutputReader(CLOUDGENE_LOG);
		assertTrue(log.hasInMemory("3 valid VCF file(s) found."));
		assertTrue(log.hasInMemory("Samples: 41"));
		assertTrue(log.hasInMemory("Chromosomes: 2 3 4"));
		assertTrue(log.hasInMemory("SNPs: 2715"));
		assertTrue(log.hasInMemory("Chunks: 3"));
		assertTrue(log.hasInMemory("Datatype: unphased"));
		assertTrue(log.hasInMemory("Reference Panel: hapmap2"));
		assertTrue(log.hasInMemory("Phasing: eagle"));

	}

	@Test
	public void testTabixIndexCreationChr20() throws Exception {

		String inputFolder = "test-data/data/chr20-phased";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr1/hapmap2.json");
		command.setPopulation("eur");

		assertEquals(0, (int) command.call());

		OutputReader log = new OutputReader(CLOUDGENE_LOG);
		assertTrue(log.hasInMemory("1 valid VCF file(s) found."));
		assertTrue(log.hasInMemory("Samples: 51"));

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

		String inputFolder = "test-data/data/single";

		InputValidationCommand command = buildCommand(inputFolder);
		command.setReference("test-data/configs/hapmap-chr1/hapmap2.json");
		command.setPopulation("eur");

		assertEquals(0, (int) command.call());

		OutputReader log = new OutputReader(CLOUDGENE_LOG);
		assertTrue(log.hasInMemory("1 valid VCF file(s) found."));
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

	private InputValidationCommand buildCommand(String inputFolder) {
		InputValidationCommand command = new InputValidationCommand();
		command.setMinSamples(1);
		command.setFiles(getFiles(inputFolder));
		command.setReport(CLOUDGENE_LOG);
		return command;
	}

}

package genepi.imputationserver.steps;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.junit.Test;

import genepi.imputationserver.util.AbstractTestcase;
import genepi.imputationserver.util.WorkflowTestContext;
import genepi.io.FileUtil;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class InputValidationTest extends AbstractTestcase {

	public static final boolean VERBOSE = true;
	
	@Test
	public void testInputValidation() throws Exception {
		InputValidationCommand command = new InputValidationCommand();
		command.setFiles(Arrays.asList("/home/seb/Desktop/minimac_test2.50.vcf.gz"));
		command.setReference("/home/seb/Desktop/ref.json");
		command.setOutput("cloudgene.log");
		command.setChunksize(20000000);
		command.setPhasing("eagle");
		command.setupTabix("files/bin/tabix");
		assertEquals(0, (int) command.call());
	}
	
	/*
	@Test(expected = IOException.class)
	public void testWithWrongReferencePanel() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		// create workflow context
		buildContext(inputFolder, configFolder, "missing-reference-panel");

	}

	@Test
	public void testHg19DataWithBuild38() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/three";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "hapmap2");
		context.setInput("build", "hg38");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(false, result);

		// check analyze task
		assertTrue(context.hasInMemory("[RUN] Analyze file minimac_test2.50.vcf.gz"));

		// check error message
		assertTrue(context.hasInMemory("[ERROR]"));
		assertTrue(context.hasInMemory("This is not a valid hg38 encoding."));

	}

	public void testHg38DataWithBuild19() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr20";
		String inputFolder = "test-data/data/chr20-unphased-hg38";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "hapmap2");
		context.setInput("build", "hg19");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(false, result);

		// check analyze task
		assertTrue(context.hasInMemory("[RUN] Analyze file chr20.R50.merged.1.330k.recode.unphased.small.hg38.vcf.gz"));

		// check error message
		// check error message
		assertTrue(context.hasInMemory("[ERROR]"));
		assertTrue(context.hasInMemory("This is not a valid hg19 encoding."));

	}

	@Test
	public void testWrongFiles() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		// input folder contains no vcf or vcf.gz files
		String inputFolder = "test-data/data/wrong_files";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "hapmap2");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(false, result);

		// check error message
		assertTrue(context.hasInMemory("[ERROR] The provided files are not VCF files"));

	}

	@Test
	public void testWrongVcfFile() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		// input folder contains no vcf or vcf.gz files
		String inputFolder = "test-data/data/wrong_vcf";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "hapmap2");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(false, result);

		// check error message
		assertTrue(context.hasInMemory("[ERROR] Unable to parse header with error"));

	}

	@Test
	public void testMixedPopulation() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "hrc-fake");

		context.setInput("phasing", "eagle");
		context.setInput("population", "mixed");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(true, result);

	}

	@Test
	public void testCorrectHrcPopulation() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "hrc-fake");

		context.setInput("phasing", "eagle");
		context.setInput("population", "eur");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(true, result);

	}

	@Test
	public void testWrongHrcPopulation() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "hrc-fake");

		context.setInput("phasing", "eagle");
		context.setInput("population", "aas");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(false, result);

		assertTrue(context.hasInMemory("[ERROR] Population 'aas' is not supported by reference panel 'hrc-fake'"));

	}

	@Test
	public void testWrong1KP3Population() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "phase3-fake");

		context.setInput("phasing", "eagle");
		context.setInput("population", "asn");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(false, result);

		assertTrue(context.hasInMemory("[ERROR] Population 'asn' is not supported by reference panel 'phase3-fake'"));

	}

	@Test
	public void testWrongTopmedPopulation() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "TOPMedfreeze6-fake");

		context.setInput("phasing", "eagle");
		context.setInput("population", "asn");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(false, result);

		assertTrue(context
				.hasInMemory("[ERROR] Population 'asn' is not supported by reference panel 'TOPMedfreeze6-fake'"));

	}

	@Test
	public void testUnorderedVcfFile() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		// input folder contains no vcf or vcf.gz files
		String inputFolder = "test-data/data/unorderd";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "hapmap2");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(false, result);

		// check error message
		assertTrue(context.hasInMemory(" [ERROR] The provided VCF file is malformed"));
		assertTrue(context.hasInMemory(" Error during index creation"));

	}

	@Test
	public void testWrongChromosomes() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/wrong_chrs";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "hapmap2");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(false, result);

		// check analyze task
		assertTrue(context.hasInMemory("[RUN] Analyze file minimac_test.50.vcf.gz"));

		// check error message
		assertTrue(context.hasInMemory("[ERROR] The provided VCF file contains more than one chromosome."));

	}

	@Test
	public void testSingleUnphasedVcfWithEagle() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "hapmap2");

		context.setInput("phasing", "eagle");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(true, result);

		// check analyze task and results
		assertTrue(context.hasInMemory("[RUN] Analyze file minimac_test.50.vcf.gz"));
		assertTrue(context.hasInMemory("[OK] 1 valid VCF file(s) found."));

		// check statistics
		assertTrue(context.hasInMemory("Samples: 41"));
		assertTrue(context.hasInMemory("Chromosomes: 1"));
		assertTrue(context.hasInMemory("SNPs: 905"));
		assertTrue(context.hasInMemory("Chunks: 1"));
		assertTrue(context.hasInMemory("Datatype: unphased"));
		assertTrue(context.hasInMemory("Reference Panel: hapmap2"));
		assertTrue(context.hasInMemory("Phasing: eagle"));

	}

	@Test
	public void testThreeUnphasedVcfWithEagle() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/three";
		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "hapmap2");
		context.setInput("phasing", "eagle");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(true, result);

		// check analyze task and results
		assertTrue(context.hasInMemory("[RUN] Analyze file minimac_test2.50.vcf.gz"));
		assertTrue(context.hasInMemory("[RUN] Analyze file minimac_test3.50.vcf.gz"));
		assertTrue(context.hasInMemory("[RUN] Analyze file minimac_test4.50.vcf.gz"));
		assertTrue(context.hasInMemory("[OK] 3 valid VCF file(s) found."));

		// check statistics
		assertTrue(context.hasInMemory("Samples: 41"));
		assertTrue(context.hasInMemory("Chromosomes: 2 3 4"));
		assertTrue(context.hasInMemory("SNPs: 2715"));
		assertTrue(context.hasInMemory("Chunks: 3"));
		assertTrue(context.hasInMemory("Datatype: unphased"));
		assertTrue(context.hasInMemory("Reference Panel: hapmap2"));
		assertTrue(context.hasInMemory("Phasing: eagle"));

	}

	@Test
	public void testTabixIndexCreationChr20() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		// input folder contains no vcf or vcf.gz files
		String inputFolder = "test-data/data/chr20-phased";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "hapmap2");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(true, result);
		assertTrue(context.hasInMemory("[OK] 1 valid VCF file(s) found."));

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
	public void testTabixIndexCreationChr1() throws IOException {

		String configFolder = "test-data/configs/hapmap-chr1";
		// input folder contains no vcf or vcf.gz files
		String inputFolder = "test-data/data/single";

		// create workflow context
		WorkflowTestContext context = buildContext(inputFolder, configFolder, "hapmap2");
		context.setInput("phasing", "eagle");

		// create step instance
		InputValidation inputValidation = new InputValidationMock(configFolder);

		// run and test
		boolean result = run(context, inputValidation);

		// check if step is failed
		assertEquals(true, result);
		assertTrue(context.hasInMemory("[OK] 1 valid VCF file(s) found."));

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
	
	*/

	

}

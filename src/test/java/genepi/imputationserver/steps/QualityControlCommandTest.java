package genepi.imputationserver.steps;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.Arrays;
import java.util.Map;

import org.junit.Test;

import genepi.imputationserver.util.AbstractTestcase;
import genepi.imputationserver.util.RefPanelUtil;
import genepi.imputationserver.util.report.CloudgeneReport;
import genepi.io.FileUtil;
import genepi.io.text.LineReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import net.lingala.zip4j.exception.ZipException;

public class QualityControlCommandTest extends AbstractTestcase {

	public static final boolean VERBOSE = true;

	private static final String TABIX_HOME = "files/bin/tabix";

	private static final String CLOUDGENE_LOG = "cloudgene.report.json";

	private static final String TEST_DATA_TMP = "test-data/tmp";

	private QualityControlCommand buildCommand(String inputFolder) {

		File tmp = new File(TEST_DATA_TMP);
		if (tmp.exists()) {
			FileUtil.deleteDirectory(tmp);
		}
		tmp.mkdirs();

		QualityControlCommand command = new QualityControlCommand();
		command.setFiles(getFiles(inputFolder));
		command.setupTabix(TABIX_HOME);
		command.setMafOutput(FileUtil.path(TEST_DATA_TMP, "maf.txt"));
		command.setMetafilesOutput(TEST_DATA_TMP);
		command.setChunksOutput(TEST_DATA_TMP);
		command.setStatisticsOutput(TEST_DATA_TMP);
		return command;
	}

	@Test
	public void testQcStatistics() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		assertEquals(-1, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		assertTrue(context.hasInMemory("Alternative allele frequency > 0.5 sites: 185"));
		assertTrue(context.hasInMemory("Excluded sites in total: 336"));
		assertTrue(context.hasInMemory("Remaining sites in total: 96"));
		assertTrue(context.hasInMemory("Monomorphic sites: 331"));

	}

	@Test
	public void testQcStatisticAllChunksExcluded() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");
		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);

		assertEquals(-1, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);
		// check statistics
		assertTrue(context.hasInMemory("No chunks passed the QC step"));

	}

	@Test
	public void testQcStatisticsAllChunksFailed() throws Exception {

		String configFolder = "test-data/configs/hapmap-3chr";
		String inputFolder = "test-data/data/simulated-chip-3chr-qc";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);

		assertEquals(-1, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		// check statistics
		assertTrue(context.hasInMemory("Alternative allele frequency > 0.5 sites: 37,503"));
		assertTrue(context.hasInMemory("Duplicated sites: 618"));
		assertTrue(context.hasInMemory("36 Chunk(s) excluded"));
		assertTrue(context.hasInMemory("No chunks passed the QC step"));

	}

	@Test
	public void testCountLinesInFailedChunkFile() throws Exception {

		String configFolder = "test-data/configs/hapmap-3chr";
		String inputFolder = "test-data/data/simulated-chip-3chr-qc";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);

		assertEquals(-1, (int) command.call());

		LineReader reader = new LineReader(FileUtil.path(TEST_DATA_TMP, "chunks-excluded.txt"));

		int count = 0;
		String testLine = null;

		while (reader.next()) {
			count++;
			if (count == 8) {
				testLine = reader.get();
			}
		}

		assertEquals(37, count);

		assertEquals("chunk_1_0120000001_0140000000" + "\t" + "108" + "\t" + "0.9391304347826087" + "\t" + "19",
				testLine);

	}

	@Test
	public void testQcStatisticsAllChunksPassed() throws Exception {

		String configFolder = "test-data/configs/hapmap-3chr";
		String inputFolder = "test-data/data/simulated-chip-3chr-imputation";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");
		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		assertEquals(0, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		// check statistics
		assertTrue(context.hasInMemory("Excluded sites in total: 3,058"));
		assertTrue(context.hasInMemory("Remaining sites in total: 117,498"));

	}

	@Test
	public void testCountSitesForOneChunkedContig() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/simulated-chip-1chr-imputation";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);

		assertEquals(0, (int) command.call());

		File[] files = new File(TEST_DATA_TMP).listFiles();
		Arrays.sort(files);
		// baseline from a earlier job execution
		int[] array = { 4750, 5174, 5106, 5832, 5318, 4588, 968, 3002, 5781, 5116, 5699, 6334, 3188 };
		int pos = 0;

		for (File file : files) {
			int count = 0;
			if (file.getName().endsWith(".gz")) {
				VCFFileReader vcfReader = new VCFFileReader(file, false);
				CloseableIterator<VariantContext> it = vcfReader.iterator();
				while (it.hasNext()) {
					it.next();
					count++;
				}
				assertEquals(array[pos], count);
				vcfReader.close();
				pos++;
			}

		}

	}

	@Test
	public void testCountAmountSplitsForSeveralContigs() throws Exception {

		String configFolder = "test-data/configs/hapmap-3chr";
		String inputFolder = "test-data/data/simulated-chip-3chr-imputation";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);

		assertEquals(0, (int) command.call());

		File[] files = new File(TEST_DATA_TMP).listFiles();
		Arrays.sort(files);

		// baseline from a earlier job execution

		int count = 0;
		for (File file : files) {
			if (file.getName().endsWith(".gz")) {
				count++;
			}

		}
		// https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
		assertEquals(13 + 13 + 10, count);

	}

	@Test
	public void testCountLinesInChunkMetaFile() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/simulated-chip-1chr-imputation";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		assertEquals(0, (int) command.call());

		LineReader reader = new LineReader(FileUtil.path(TEST_DATA_TMP, "1"));

		int count = 0;
		while (reader.next()) {
			count++;
		}

		assertEquals(13, count);

	}

	@Test
	public void testCountSamplesInCreatedChunk() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/simulated-chip-1chr-imputation";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		assertEquals(0, (int) command.call());

		for (File file : new File(TEST_DATA_TMP).listFiles()) {
			if (file.getName().endsWith("chunk_1_80000001_100000000.vcf.gz")) {
				VCFFileReader vcfReader = new VCFFileReader(file, false);
				CloseableIterator<VariantContext> it = vcfReader.iterator();
				if (it.hasNext()) {
					VariantContext a = it.next();
					assertEquals(255, a.getNSamples());
				}
				vcfReader.close();
			}

		}

	}

	@Test
	public void testMonomorphicSnps() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr20";
		String inputFolder1 = "test-data/data/chr20-phased-1sample";
		String inputFolder50 = "test-data/data/chr20-phased";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder1);
		command.setReference(panel);
		assertEquals(0, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		assertTrue(context.hasInMemory("Monomorphic sites: 0"));

		command.setFiles(getFiles(inputFolder50));

		assertEquals(0, (int) command.call());

		context = new CloudgeneReport(CLOUDGENE_LOG);

		assertTrue(context.hasInMemory("Monomorphic sites: 11"));

	}

	@Test
	public void testChrXSplits() throws Exception, ZipException {

		String configFolder = "test-data/configs/hapmap-chrX";
		String inputFolder = "test-data/data/chrX-unphased";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "phase1");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);

		assertEquals(0, (int) command.call());

		int count = 0;
		for (File file2 : new File(TEST_DATA_TMP).listFiles()) {
			if (file2.getName().endsWith("vcf.gz")) {
				count++;
			}
		}

		assertEquals(13, count);

	}

	@Test
	public void testChrXInvalidAlleles() throws Exception, ZipException {

		String configFolder = "test-data/configs/hapmap-chrX";
		String inputFolder = "test-data/data/chrX-phased-invalid";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "phase1");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);

		assertEquals(0, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		assertTrue(context.hasInMemory("Invalid alleles: 190"));

	}

	@Test
	public void testChrXMixedGenotypes() throws Exception, ZipException {

		String configFolder = "test-data/configs/hapmap-chrX";
		String inputFolder = "test-data/data/chrX-unphased-mixed";
		File tmp = new File(TEST_DATA_TMP);

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "phase1");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		assertEquals(-1, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);
		assertTrue(context.hasInMemory("Chromosome X nonPAR region includes > 10 % mixed genotypes."));

	}

	@Test
	public void testChrXPloidyError() throws Exception, ZipException {

		String configFolder = "test-data/configs/hapmap-chrX";
		String inputFolder = "test-data/data/chrX-unphased-ploidy";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "phase1");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		assertEquals(-1, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		assertTrue(context.hasInMemory("ChrX nonPAR region includes ambiguous samples"));

	}

	@Test
	public void testAlleleFrequencyCheckWithWrongPopulation() throws Exception {
		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/single";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		command.setPopulation("afr");
		assertEquals(-1, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		// check statistics
		assertTrue(context.hasInMemory("Population 'afr' is not supported by reference panel 'hapmap2'."));
	}

	@Test
	public void testAlleleFrequencyCheckWithNoSamplesForPopulation() throws Exception {

		String configFolder = "test-data/configs/hapmap-3chr";
		String inputFolder = "test-data/data/simulated-chip-3chr-imputation";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		assertEquals(0, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		// check statistics
		assertTrue(context.hasInMemory("[WARN] Skip allele frequency check."));
	}

	@Test
	public void testQcStatisticsAllowStrandFlips() throws Exception {

		String configFolder = "test-data/configs/hapmap-3chr";
		String inputFolder = "test-data/data/simulated-chip-3chr-imputation";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		assertEquals(0, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		// check statistics
		assertTrue(context.hasInMemory("Excluded sites in total: 3,058"));
		assertTrue(context.hasInMemory("Remaining sites in total: 117,498"));

	}

	@Test
	public void testQcStatisticsDontAllowStrandFlips() throws Exception {

		String configFolder = "test-data/configs/hapmap-3chr";
		String inputFolder = "test-data/data/simulated-chip-3chr-imputation";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt",
				"hapmap2-qcfilter-strandflips");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		assertEquals(-1, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		// check statistics
		assertTrue(context.hasInMemory("Excluded sites in total: 3,058"));
		assertTrue(context.hasInMemory("Remaining sites in total: 117,498"));
		assertTrue(context.hasInMemory(
				"<b>Error:</b> More than -1 obvious strand flips have been detected. Please check strand. Imputation cannot be started!"));

	}

	@Test
	public void testQcStatisticsFilterOverlap() throws Exception {

		String configFolder = "test-data/configs/hapmap-3chr";
		String inputFolder = "test-data/data/simulated-chip-3chr-imputation";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt",
				"hapmap2-qcfilter-ref-overlap");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		command.call();

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		// check statistics
		assertTrue(context.hasInMemory("<b>Warning:</b> 36 Chunk(s) excluded: reference overlap < 99.0%"));

	}

	@Test
	public void testQcStatisticsFilterMinSnps() throws Exception {

		String configFolder = "test-data/configs/hapmap-3chr";
		String inputFolder = "test-data/data/simulated-chip-3chr-imputation";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt",
				"hapmap2-qcfilter-min-snps");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		assertEquals(0, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);
		// check statistics
		assertTrue(context.hasInMemory("<b>Warning:</b> 2 Chunk(s) excluded: < 1000 SNPs"));

	}

	@Test
	public void testQcStatisticsFilterSampleCallrate() throws Exception {

		String configFolder = "test-data/configs/hapmap-3chr";
		String inputFolder = "test-data/data/simulated-chip-3chr-imputation";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt",
				"hapmap2-qcfilter-low-callrate");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		command.call();

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		// check statistics
		assertTrue(context
				.hasInMemory("<b>Warning:</b> 36 Chunk(s) excluded: at least one sample has a call rate < 101.0%"));

	}

	@Test
	public void testChr23PipelineLifting() throws Exception, ZipException {

		String configFolder = "test-data/configs/hapmap-chrX-hg38";
		String inputFolder = "test-data/data/chr23-unphased";

		// maybe git large files?
		if (!new File(
				"test-data/configs/hapmap-chrX-hg38/ref-panels/ALL.X.nonPAR.phase1_v3.snps_indels_svs.genotypes.all.noSingleton.recode.hg38.bcf")
				.exists()) {
			System.out.println("chrX bcf nonPAR file not available");
			return;
		}

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "phase1");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		assertEquals(0, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		assertTrue(context.hasInMemory("Alternative allele frequency > 0.5 sites: 8,973"));
		assertTrue(context.hasInMemory("[MESSAGE] [WARN] Excluded sites in total: 18,076"));

	}

	@Test
	public void testChrXPipelineLifting() throws Exception, ZipException {

		String configFolder = "test-data/configs/hapmap-chrX-hg38";
		String inputFolder = "test-data/data/chrX-unphased";

		// maybe git large files?
		if (!new File(
				"test-data/configs/hapmap-chrX-hg38/ref-panels/ALL.X.nonPAR.phase1_v3.snps_indels_svs.genotypes.all.noSingleton.recode.hg38.bcf")
				.exists()) {
			System.out.println("chrX bcf nonPAR file not available");
			return;
		}

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "phase1");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		assertEquals(0, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);
		assertTrue(context.hasInMemory("Alternative allele frequency > 0.5 sites: 8,973"));
		assertTrue(context.hasInMemory("[MESSAGE] [WARN] Excluded sites in total: 18,076"));

	}

	@Test
	public void testRegionImputationSimple() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/simulated-chip-1chr-imputation";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2-region-simple");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		assertEquals(-1, (int) command.call());

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);

		// check statistics
		assertTrue(context.hasInMemory("Remaining sites in total: 1"));

	}

	@Test
	public void testRegionImputationComplex() throws Exception {

		String configFolder = "test-data/configs/hapmap-chr1";
		String inputFolder = "test-data/data/simulated-chip-1chr-imputation";

		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", "hapmap2-region-complex");

		resolveEnvVariable(panel, configFolder);

		QualityControlCommand command = buildCommand(inputFolder);
		command.setReference(panel);
		command.call();

		CloudgeneReport context = new CloudgeneReport(CLOUDGENE_LOG);
		// check statistics
		assertTrue(context.hasInMemory("Remaining sites in total: 2"));

	}

}

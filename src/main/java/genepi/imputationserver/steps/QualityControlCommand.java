package genepi.imputationserver.steps;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.Callable;

import org.apache.commons.lang.StringEscapeUtils;

import genepi.imputationserver.steps.fastqc.ITask;
import genepi.imputationserver.steps.fastqc.ITaskProgressListener;
import genepi.imputationserver.steps.fastqc.LiftOverTask;
import genepi.imputationserver.steps.fastqc.RangeEntry;
import genepi.imputationserver.steps.fastqc.StatisticsTask;
import genepi.imputationserver.steps.fastqc.TaskResults;
import genepi.imputationserver.steps.vcf.VcfFileUtil;
import genepi.imputationserver.util.RefPanel;
import genepi.imputationserver.util.RefPanelPopulation;
import genepi.imputationserver.util.report.CloudgeneReport;
import genepi.io.FileUtil;
import genepi.io.text.LineWriter;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

@Command
public class QualityControlCommand implements Callable<Integer> {

	@Parameters(description = "VCF files")
	private List<String> files;

	@Option(names = "--reference", description = "Reference Panel", required = true)
	private String reference;

	@Option(names = "--population", description = "Reference Population", required = true)
	private String population;

	@Option(names = "--maf-output", description = "MAF file", required = true)
	private String mafOutput = "";

	@Option(names = "--chunks-out", description = "VCF Chunks Folder", required = true)
	private String chunksOutput = "";

	@Option(names = "--statistics-out", description = "Statistics Folder", required = true)
	private String statisticsOutput = "";

	@Option(names = "--metafiles-out", description = "Metafiles Folder", required = true)
	private String metafilesOutput = "";

	@Option(names = "--build", description = "Build", required = false)
	private String build = "hg19";

	@Option(names = "--chain", description = "chainFile", required = false)
	private String chainFile = "";

	@Option(names = "--chunksize", description = "VCF chunksize", required = false)
	private int chunksize = 20000000;

	@Option(names = "--phasing-window", description = "Phasing window", required = false)
	private int phasingWindow = 5000000;

	@Option(names = "--report", description = "Cloudgene Report Output", required = false)
	private String report = "cloudgene.report.json";

	private CloudgeneReport context = new CloudgeneReport();

	private RefPanel panel = null;

	public QualityControlCommand() {
		VcfFileUtil.setTabixBinary("tabix");
	}

	public void setFiles(List<String> files) {
		this.files = files;
	}

	public void setRefPanel(RefPanel panel) {
		this.panel = panel;
	}

	public void setChainFile(String chainFile) {
		this.chainFile = chainFile;
	}

	public void setPopulation(String population) {
		this.population = population;
	}

	public void setMafOutput(String mafOutput) {
		this.mafOutput = mafOutput;
	}

	public void setChunksOutput(String chunksOutput) {
		this.chunksOutput = chunksOutput;
	}

	public void setStatisticsOutput(String statisticsOutput) {
		this.statisticsOutput = statisticsOutput;
	}

	public void setMetafilesOutput(String metafilesOutput) {
		this.metafilesOutput = metafilesOutput;
	}

	public void setReport(String report) {
		this.report = report;
	}

	public String getReport() {
		return report;
	}

	protected void setupTabix(String path) {
		VcfFileUtil.setTabixBinary(path);
	}

	@Override
	public Integer call() throws Exception {

		context.setFilename(report);

		if (panel == null) {
			try {
				panel = RefPanel.loadFromJson(reference);
				if (panel == null) {
					context.error("Reference not found.");
					return -1;
				}
			} catch (Exception e) {
				context.error("Unable to parse reference panel: " + StringEscapeUtils.escapeHtml(e.getMessage()));
				return -1;
			}
		}

		if (!analyzeFiles()) {
			return -1;
		} else {
			return 0;
		}

	}

	private boolean analyzeFiles() {
		Collections.sort(files);

		LineWriter excludedSnpsWriter = null;

		String excludedSnpsFile = FileUtil.path(statisticsOutput, "snps-excluded.txt");

		try {
			excludedSnpsWriter = new LineWriter(excludedSnpsFile);
			excludedSnpsWriter.write("#Position" + "\t" + "FilterType" + "\t" + " Info", false);
		} catch (Exception e) {
			context.error("Error creating file writer");
			return false;
		}

		String[] vcfFilenames = files.toArray(new String[files.size()]);
		// check if liftover is needed
		if (!build.equals(panel.getBuild())) {
			context.warning("Uploaded data is " + build + " and reference is " + panel.getBuild() + ".");

			if (chainFile == null) {
				context.error("Currently we do not support liftOver from " + build + " to " + panel.getBuild());
				return false;
			}

			String fullPathChainFile = FileUtil.path(chainFile);
			if (!new File(fullPathChainFile).exists()) {
				context.error("Chain file " + fullPathChainFile + " not found.");
				return false;
			}

			LiftOverTask task = new LiftOverTask();
			task.setVcfFilenames(vcfFilenames);
			task.setChainFile(fullPathChainFile);
			task.setChunksDir(chunksOutput);
			task.setExcludedSnpsWriter(excludedSnpsWriter);

			TaskResults results = runTask(context, task);

			if (results.isSuccess()) {
				vcfFilenames = task.getNewVcfFilenames();
			} else {
				return false;
			}

		}

		// calculate statistics
		StatisticsTask task = new StatisticsTask();
		task.setVcfFilenames(vcfFilenames);
		task.setExcludedSnpsWriter(excludedSnpsWriter);
		task.setChunkSize(chunksize);
		task.setPhasingWindow(phasingWindow);
		task.setPopulation(population);
		// support relative path
		String legend = panel.getLegend();
		if (!legend.startsWith("/") && !legend.startsWith("./")) {
			legend = FileUtil.path(legend);
		}

		// check chromosomes
		if (!panel.supportsPopulation(population)) {
			StringBuilder report = new StringBuilder();
			report.append(
					"Population '" + population + "' is not supported by reference panel '" + panel.getId() + "'.\n");
			if (panel.getPopulations() != null) {
				report.append("Available populations:");
				for (RefPanelPopulation pop : panel.getPopulations()) {
					report.append("\n - " + pop.getId());
				}
			}
			context.error(report.toString());
			return false;
		}

		int refSamples = panel.getSamplesByPopulation(population);
		if (refSamples <= 0) {
			context.warning("Skip allele frequency check.");
			task.setAlleleFrequencyCheck(false);
		}

		task.setLegendFile(legend);
		task.setRefSamples(refSamples);
		task.setMafFile(mafOutput);
		task.setChunkFileDir(metafilesOutput);
		task.setChunksDir(chunksOutput);
		task.setStatDir(statisticsOutput);
		task.setBuild(panel.getBuild());

		double referenceOverlap = panel.getQcFilterByKey("overlap");
		int minSnps = (int) panel.getQcFilterByKey("minSnps");
		double sampleCallrate = panel.getQcFilterByKey("sampleCallrate");
		double mixedGenotypesChrX = panel.getQcFilterByKey("mixedGenotypeschrX");
		int strandFlips = (int) (panel.getQcFilterByKey("strandFlips"));
		int alleleSwitches = (int) (panel.getQcFilterByKey("alleleSwitches"));
		String ranges = panel.getRange();
		
		if (ranges != null) {
			HashSet<RangeEntry> rangeEntries = new HashSet<RangeEntry>();

			for (String range : panel.getRange().split(",")) {
				String chromosome = range.split(":")[0].trim();
				String region = range.split(":")[1].trim();
				int start = Integer.valueOf(region.split("-")[0].trim());
				int end = Integer.valueOf(region.split("-")[1].trim());
				RangeEntry entry = new RangeEntry();
				entry.setChromosome(chromosome);
				entry.setStart(start);
				entry.setEnd(end);
				rangeEntries.add(entry);
			}

			task.setRanges(rangeEntries);
			context.log("Reference Panel Ranges: " + rangeEntries);
		} else {
			context.log("Reference Panel Ranges: genome-wide");
		}

		task.setReferenceOverlap(referenceOverlap);
		task.setMinSnps(minSnps);
		task.setSampleCallrate(sampleCallrate);
		task.setMixedGenotypeschrX(mixedGenotypesChrX);

		TaskResults results = runTask(context, task);

		if (!results.isSuccess()) {
			return false;
		}

		try {

			excludedSnpsWriter.close();

			if (!excludedSnpsWriter.hasData()) {
				FileUtil.deleteFile(excludedSnpsFile);
			}

		} catch (IOException e) {
			e.printStackTrace();
		}

		DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.US);
		DecimalFormat df = new DecimalFormat("#.00", symbols);
		DecimalFormat formatter = new DecimalFormat("###,###.###", symbols);

		StringBuffer text = new StringBuffer();

		text.append("<b>Statistics:</b> \n");
		if (ranges != null) {
			text.append("Ref. Panel Range: " + ranges + "\n");
		}
		text.append(
				"Alternative allele frequency > 0.5 sites: " + formatter.format(task.getAlternativeAlleles()) + "\n");
		text.append("Reference Overlap: "
				+ df.format(
						task.getFoundInLegend() / (double) (task.getFoundInLegend() + task.getNotFoundInLegend()) * 100)
				+ " %" + "\n");

		text.append("Match: " + formatter.format(task.getMatch()) + "\n");
		text.append("Allele switch: " + formatter.format(task.getAlleleSwitch()) + "\n");
		text.append("Strand flip: " + formatter.format(task.getStrandFlipSimple()) + "\n");
		text.append("Strand flip and allele switch: " + formatter.format(task.getStrandFlipAndAlleleSwitch()) + "\n");
		text.append("A/T, C/G genotypes: " + formatter.format(task.getComplicatedGenotypes()) + "\n");

		text.append("<b>Filtered sites:</b> \n");
		text.append("Filter flag set: " + formatter.format(task.getFilterFlag()) + "\n");
		text.append("Invalid alleles: " + formatter.format(task.getInvalidAlleles()) + "\n");
		text.append("Multiallelic sites: " + formatter.format(task.getMultiallelicSites()) + "\n");
		text.append("Duplicated sites: " + formatter.format(task.getDuplicates()) + "\n");
		text.append("NonSNP sites: " + formatter.format(task.getNoSnps()) + "\n");
		text.append("Monomorphic sites: " + formatter.format(task.getMonomorphic()) + "\n");
		text.append("Allele mismatch: " + formatter.format(task.getAlleleMismatch()) + "\n");
		text.append("SNPs call rate < 90%: " + formatter.format(task.getLowCallRate()));

		context.ok(text.toString());

		text = new StringBuffer();

		text.append("Excluded sites in total: " + formatter.format(task.getFiltered()) + "\n");
		text.append("Remaining sites in total: " + formatter.format(task.getOverallSnps()) + "\n");

		if (task.getFiltered() > 0) {
			text.append("See snps-excluded.txt for details" + "\n");
		}

		if (task.getNotFoundInLegend() > 0) {
			text.append("Typed only sites: " + formatter.format(task.getNotFoundInLegend()) + "\n");
			text.append("See typed-only.txt for details" + "\n");
		}

		if (task.getRemovedChunksSnps() > 0) {

			text.append("\n<b>Warning:</b> " + formatter.format(task.getRemovedChunksSnps())

					+ " Chunk(s) excluded: < " + minSnps + " SNPs (see chunks-excluded.txt for details).");
		}

		if (task.getRemovedChunksCallRate() > 0) {

			text.append("\n<b>Warning:</b> " + formatter.format(task.getRemovedChunksCallRate())

					+ " Chunk(s) excluded: at least one sample has a call rate < " + (sampleCallrate * 100) + "% (see "
					+ "chunks-excluded.txt for details).");
		}

		if (task.getRemovedChunksOverlap() > 0) {

			text.append("\n<b>Warning:</b> " + formatter.format(task.getRemovedChunksOverlap())

					+ " Chunk(s) excluded: reference overlap < " + (referenceOverlap * 100) + "% (see "
					+ " chunks-excluded.txt for details).");
		}

		long excludedChunks = task.getRemovedChunksSnps() + task.getRemovedChunksCallRate()
				+ task.getRemovedChunksOverlap();

		long overallChunks = task.getOverallChunks();

		if (excludedChunks > 0) {
			text.append("\nRemaining chunk(s): " + formatter.format(overallChunks - excludedChunks));

		}

		if (excludedChunks == overallChunks) {

			text.append("\n<b>Error:</b> No chunks passed the QC step. Imputation cannot be started!");
			context.error(text.toString());

			return false;

		}
		// strand flips (normal flip & allele switch + strand flip)
		else if (task.getStrandFlipSimple() + task.getStrandFlipAndAlleleSwitch() > strandFlips) {
			text.append("\n<b>Error:</b> More than " + strandFlips
					+ " obvious strand flips have been detected. Please check strand. Imputation cannot be started!");
			context.error(text.toString());

			return false;
		}

		// Check if too many allele switches are detected
		else if (task.getAlleleSwitch() + task.getStrandFlipAndAlleleSwitch() > alleleSwitches) {
			text.append("<br><b>Error:</b> More than " + alleleSwitches
					+ " allele switches have been detected. Imputation cannot be started!");
			context.error(text.toString());

			return false;
		}

		else if (task.isChrXMissingRate()) {
			text.append(
					"\n<b>Error:</b> Chromosome X nonPAR region includes > 10 % mixed genotypes. Imputation cannot be started!");
			context.error(text.toString());

			return false;
		}

		else if (task.isChrXPloidyError()) {
			text.append(
					"\n<b>Error:</b> ChrX nonPAR region includes ambiguous samples (haploid and diploid positions). Imputation cannot be started! See "
							+ "chrX-info.txt");
			context.error(text.toString());

			return false;
		}

		else {

			text.append(results.getMessage());
			context.warning(text.toString());
			return true;

		}
	}

	protected TaskResults runTask(final CloudgeneReport context2, ITask task) {
		context2.beginTask("Running " + task.getName() + "...");
		TaskResults results;
		try {
			results = task.run(new ITaskProgressListener() {

				@Override
				public void progress(String message) {
					context2.updateTask(message, CloudgeneReport.RUNNING);

				}
			});

			if (results.isSuccess()) {
				context2.endTask(task.getName(), CloudgeneReport.OK);
			} else {
				context2.endTask(task.getName() + "\n" + results.getMessage(), CloudgeneReport.ERROR);
			}
			return results;
		} catch (Exception | Error e) {
			e.printStackTrace();
			TaskResults result = new TaskResults();
			result.setSuccess(false);
			result.setMessage(e.getMessage());
			StringWriter s = new StringWriter();
			e.printStackTrace(new PrintWriter(s));
			context2.println("Task '" + task.getName() + "' failed.\nException:" + s.toString());
			context2.endTask(e.getMessage(), CloudgeneReport.ERROR);
			return result;
		}

	}

}

package genepi.imputationserver.steps;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Vector;
import java.util.concurrent.Callable;
import genepi.imputationserver.steps.fastqc.ITask;
import genepi.imputationserver.steps.fastqc.LiftOverTask;
import genepi.imputationserver.steps.fastqc.RangeEntry;
import genepi.imputationserver.steps.fastqc.StatisticsTask;
import genepi.imputationserver.steps.fastqc.TaskResults;
import genepi.imputationserver.steps.vcf.VcfFileUtil;
import genepi.imputationserver.util.OutputWriter;
import genepi.imputationserver.util.RefPanel;
import genepi.imputationserver.util.RefPanelPopulation;
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
	private String report = null;

	private OutputWriter output = null;

	private RefPanel panel = null;

	public QualityControlCommand() {
		VcfFileUtil.setTabixBinary("tabix");
	}

	public void setFiles(List<String> files) {
		this.files = files;
	}

	public void setChainFile(String chainFile) {
		this.chainFile = chainFile;
	}

	public void setReference(String reference) {
		this.reference = reference;
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

		if (report != null) {
			output = new OutputWriter(report);
		} else {
			output = new OutputWriter();
		}

		try {
			panel = RefPanel.loadFromJson(reference);
		} catch (Exception e) {
			output.error("Unable to parse reference panel:", e);
			return -1;
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
			output.error("Error creating file writer:", e);
			return false;
		}

		String[] vcfFilenames = files.toArray(new String[files.size()]);
		// check if liftover is needed
		if (!build.equals(panel.getBuild())) {
			output.warning("Uploaded data is " + build + " and reference is " + panel.getBuild() + ".");

			if (chainFile == null) {
				output.error("Currently we do not support liftOver from " + build + " to " + panel.getBuild());
				return false;
			}

			String fullPathChainFile = FileUtil.path(chainFile);
			if (!new File(fullPathChainFile).exists()) {
				output.error("Chain file " + fullPathChainFile + " not found.");
				return false;
			}

			LiftOverTask task = new LiftOverTask();
			task.setVcfFilenames(vcfFilenames);
			task.setChainFile(fullPathChainFile);
			task.setChunksDir(chunksOutput);
			task.setExcludedSnpsWriter(excludedSnpsWriter);

			TaskResults results = runTask(output, task);

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
		String sites = panel.getSites();

		// check chromosomes
		if (!panel.supportsPopulation(population)) {
			List<String> report = new Vector<String>();
			report.add("Population '" + population + "' is not supported by reference panel '" + panel.getId() + "'.");
			if (panel.getPopulations() != null) {
				report.add("Available populations:");
				for (RefPanelPopulation pop : panel.getPopulations()) {
					report.add(" - " + pop.getId());
				}
			}
			output.error(report);
			return false;
		}

		int refSamples = panel.getSamplesByPopulation(population);
		if (refSamples <= 0) {
			output.warning("Skip allele frequency check.");
			task.setAlleleFrequencyCheck(false);
		}

		task.setSitesFile(sites);
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
			output.log("Reference Panel Ranges: " + rangeEntries);
		} else {
			output.log("Reference Panel Ranges: genome-wide");
		}

		task.setReferenceOverlap(referenceOverlap);
		task.setMinSnps(minSnps);
		task.setSampleCallrate(sampleCallrate);
		task.setMixedGenotypeschrX(mixedGenotypesChrX);

		TaskResults results = runTask(output, task);

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

		List<String> text = new Vector<String>();

		text.add("<b>Statistics:</b>");
		if (ranges != null) {
			text.add("Ref. Panel Range: " + ranges);
		}
		text.add("Alternative allele frequency > 0.5 sites: " + formatter.format(task.getAlternativeAlleles()));
		text.add("Reference Overlap: "
				+ df.format(
						task.getFoundInLegend() / (double) (task.getFoundInLegend() + task.getNotFoundInLegend()) * 100)
				+ " %");

		text.add("Match: " + formatter.format(task.getMatch()));
		text.add("Allele switch: " + formatter.format(task.getAlleleSwitch()));
		text.add("Strand flip: " + formatter.format(task.getStrandFlipSimple()));
		text.add("Strand flip and allele switch: " + formatter.format(task.getStrandFlipAndAlleleSwitch()));
		text.add("A/T, C/G genotypes: " + formatter.format(task.getComplicatedGenotypes()));

		text.add("<b>Filtered sites:</b>");
		text.add("Filter flag set: " + formatter.format(task.getFilterFlag()));
		text.add("Invalid alleles: " + formatter.format(task.getInvalidAlleles()));
		text.add("Multiallelic sites: " + formatter.format(task.getMultiallelicSites()));
		text.add("Duplicated sites: " + formatter.format(task.getDuplicates()));
		text.add("NonSNP sites: " + formatter.format(task.getNoSnps()));
		text.add("Monomorphic sites: " + formatter.format(task.getMonomorphic()));
		text.add("Allele mismatch: " + formatter.format(task.getAlleleMismatch()));
		text.add("SNPs call rate < 90%: " + formatter.format(task.getLowCallRate()));

		output.message(text);

		text = new Vector<String>();

		text.add("Excluded sites in total: " + formatter.format(task.getFiltered()));
		text.add("Remaining sites in total: " + formatter.format(task.getOverallSnps()));

		if (task.getFiltered() > 0) {
			text.add("See snps-excluded.txt for details");
		}

		if (task.getNotFoundInLegend() > 0) {
			text.add("Typed only sites: " + formatter.format(task.getNotFoundInLegend()));
			text.add("See typed-only.txt for details");
		}

		if (task.getRemovedChunksSnps() > 0) {

			text.add("\n<b>Warning:</b> " + formatter.format(task.getRemovedChunksSnps())

					+ " Chunk(s) excluded: < " + minSnps + " SNPs (see chunks-excluded.txt for details).");
		}

		if (task.getRemovedChunksCallRate() > 0) {

			text.add("\n<b>Warning:</b> " + formatter.format(task.getRemovedChunksCallRate())

					+ " Chunk(s) excluded: at least one sample has a call rate < " + (sampleCallrate * 100) + "% (see "
					+ "chunks-excluded.txt for details).");
		}

		if (task.getRemovedChunksOverlap() > 0) {

			text.add("\n<b>Warning:</b> " + formatter.format(task.getRemovedChunksOverlap())

					+ " Chunk(s) excluded: reference overlap < " + (referenceOverlap * 100) + "% (see "
					+ " chunks-excluded.txt for details).");
		}

		long excludedChunks = task.getRemovedChunksSnps() + task.getRemovedChunksCallRate()
				+ task.getRemovedChunksOverlap();

		long overallChunks = task.getOverallChunks();

		if (excludedChunks > 0) {
			text.add("\nRemaining chunk(s): " + formatter.format(overallChunks - excludedChunks));

		}

		if (excludedChunks == overallChunks) {

			text.add("\n<b>Error:</b> No chunks passed the QC step. Imputation cannot be started!");
			output.error(text);

			return false;

		}
		// strand flips (normal flip & allele switch + strand flip)
		else if (task.getStrandFlipSimple() + task.getStrandFlipAndAlleleSwitch() > strandFlips) {
			text.add("\n<b>Error:</b> More than " + strandFlips
					+ " obvious strand flips have been detected. Please check strand. Imputation cannot be started!");
			output.error(text);

			return false;
		}

		// Check if too many allele switches are detected
		else if (task.getAlleleSwitch() + task.getStrandFlipAndAlleleSwitch() > alleleSwitches) {
			text.add("<br><b>Error:</b> More than " + alleleSwitches
					+ " allele switches have been detected. Imputation cannot be started!");
			output.error(text.toString());

			return false;
		}

		else if (task.isChrXMissingRate()) {
			text.add(
					"\n<b>Error:</b> Chromosome X nonPAR region includes > 10 % mixed genotypes. Imputation cannot be started!");
			output.error(text);

			return false;
		}

		else if (task.isChrXPloidyError()) {
			text.add(
					"\n<b>Error:</b> ChrX nonPAR region includes ambiguous samples (haploid and diploid positions). Imputation cannot be started! See "
							+ "chrX-info.txt");
			output.error(text);

			return false;
		}

		else {
			text.add(results.getMessage());
			
			output.warning(text);
			return true;

		}
	}

	protected TaskResults runTask(final OutputWriter output, ITask task) {

		TaskResults results;
		try {
			results = task.run();

			if (results.isSuccess()) {
				output.message(task.getName());
			} else {
				output.error(task.getName() + "\n" + results.getMessage());
			}
			return results;
		} catch (Exception e) {
			System.out.println("dfdfd " +e.getMessage());
			e.printStackTrace();
			TaskResults result = new TaskResults();
			result.setSuccess(false);
			result.setMessage(e.getMessage());
			output.error("Task '" + task.getName(), e);
			return result;
		}

	}

}

package genepi.imputationserver.steps;

import java.io.File;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Callable;
import genepi.imputationserver.steps.fastqc.ITask;
import genepi.imputationserver.steps.fastqc.LiftOverTask;
import genepi.imputationserver.steps.fastqc.RangeEntry;
import genepi.imputationserver.steps.fastqc.StatisticsTask;
import genepi.imputationserver.steps.fastqc.TaskResults;
import genepi.imputationserver.util.OutputWriter;
import genepi.imputationserver.util.RefPanel;
import genepi.imputationserver.util.RefPanelPopulation;
import genepi.imputationserver.util.StringUtils;
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
			return 1;
		}

		if (files.isEmpty()) {
			output.error("No vcf files provided.");
			return 1;
		}

		Collections.sort(files);
		String[] vcfFilenames = files.toArray(new String[files.size()]);

		if (!build.equals(panel.getBuild())) {
			vcfFilenames = liftOver(vcfFilenames);
			if (vcfFilenames == null) {
				return 1;
			}
		}

		if (analyzeFiles(vcfFilenames)) {
			return 0;
		} else {
			return 1;
		}

	}

	private String[] liftOver(String[] vcfFilenames) {

		output.warning("Uploaded data is " + build + " and reference is " + panel.getBuild() + ".");

		if (chainFile == null) {
			output.error("Currently we do not support liftOver from " + build + " to " + panel.getBuild());
			return null;
		}

		if (!new File(chainFile).exists()) {
			output.error("Chain file " + chainFile + " not found.");
			return null;
		}

		LiftOverTask task = new LiftOverTask();
		task.setVcfFilenames(vcfFilenames);
		task.setChainFile(chainFile);
		task.setChunksDir(chunksOutput);
		task.setStatDir(statisticsOutput);
		TaskResults results = runTask(output, task);

		if (results.isSuccess()) {
			return task.getNewVcfFilenames();
		} else {
			return null;
		}

	}

	private boolean analyzeFiles(String[] vcfFilenames) {

		StatisticsTask task = new StatisticsTask();
		task.setVcfFilenames(vcfFilenames);
		task.setChunkSize(chunksize);
		task.setPhasingWindow(phasingWindow);
		task.setPopulation(population);

		// check populations
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

		task.setSitesFile(panel.getSites());
		task.setRefSamples(refSamples);
		task.setMafFile(mafOutput);
		task.setChunkFileDir(metafilesOutput);
		task.setChunksDir(chunksOutput);
		task.setStatDir(statisticsOutput);
		task.setBuild(panel.getBuild());

		double referenceOverlap = panel.getQcFilterByKey("overlap");
		task.setMinReferenceOverlap(referenceOverlap);

		int minSnps = (int) panel.getQcFilterByKey("minSnps");
		task.setMinSnps(minSnps);

		double sampleCallrate = panel.getQcFilterByKey("sampleCallrate");
		task.setSampleCallrate(sampleCallrate);

		double mixedGenotypesChrX = panel.getQcFilterByKey("mixedGenotypeschrX");
		task.setMixedGenotypeschrX(mixedGenotypesChrX);

		int strandFlips = (int) (panel.getQcFilterByKey("strandFlips"));
		int alleleSwitches = (int) (panel.getQcFilterByKey("alleleSwitches"));
		
		if (panel.getRange() != null) {
			HashSet<RangeEntry> rangeEntries = parseRangeEntries(panel.getRange());
			task.setRanges(rangeEntries);
			output.log("Reference Panel Ranges: " + rangeEntries);
		} else {
			output.log("Reference Panel Ranges: genome-wide");
		}

		TaskResults results = runTask(output, task);

		if (!results.isSuccess()) {
			return false;
		}

		List<String> text = new Vector<String>();

		text.add("<b>Statistics:</b>");
		if (panel.getRange() != null) {
			text.add("Ref. Panel Range: " + panel.getRange());
		}
		text.add("Alternative allele frequency > 0.5 sites: " + StringUtils.format(task.getAlternativeAlleles()));
		text.add("Reference Overlap: " + StringUtils.format(task.getFoundInLegend() / (double) (task.getFoundInLegend() + task.getNotFoundInLegend()) * 100) + " %");
		text.add("Match: " + StringUtils.format(task.getMatch()));
		text.add("Allele switch: " + StringUtils.format(task.getAlleleSwitch()));
		text.add("Strand flip: " + StringUtils.format(task.getStrandFlipSimple()));
		text.add("Strand flip and allele switch: " + StringUtils.format(task.getStrandFlipAndAlleleSwitch()));
		text.add("A/T, C/G genotypes: " + StringUtils.format(task.getComplicatedGenotypes()));
		text.add("<b>Filtered sites:</b>");
		text.add("Filter flag set: " + StringUtils.format(task.getFilterFlag()));
		text.add("Invalid alleles: " + StringUtils.format(task.getInvalidAlleles()));
		text.add("Multiallelic sites: " + StringUtils.format(task.getMultiallelicSites()));
		text.add("Duplicated sites: " + StringUtils.format(task.getDuplicates()));
		text.add("NonSNP sites: " + StringUtils.format(task.getNoSnps()));
		text.add("Monomorphic sites: " + StringUtils.format(task.getMonomorphic()));
		text.add("Allele mismatch: " + StringUtils.format(task.getAlleleMismatch()));
		text.add("SNPs call rate < 90%: " + StringUtils.format(task.getLowCallRate()));

		output.message(text);

		text = new Vector<String>();

		text.add("Excluded sites in total: " + StringUtils.format(task.getFiltered()));
		text.add("Remaining sites in total: " + StringUtils.format(task.getOverallSnps()));

		if (task.getFiltered() > 0) {
			text.add("See snps-excluded.txt for details");
		}

		if (task.getNotFoundInLegend() > 0) {
			text.add("Typed only sites: " + StringUtils.format(task.getNotFoundInLegend()));
			text.add("See typed-only.txt for details");
		}

		if (task.getRemovedChunksSnps() > 0) {
			text.add("\n<b>Warning:</b> " + StringUtils.format(task.getRemovedChunksSnps())
					+ " Chunk(s) excluded: < " + minSnps + " SNPs (see chunks-excluded.txt for details).");
		}

		if (task.getRemovedChunksCallRate() > 0) {
			text.add("\n<b>Warning:</b> " + StringUtils.format(task.getRemovedChunksCallRate())
					+ " Chunk(s) excluded: at least one sample has a call rate < " + (sampleCallrate * 100) + "% (see "
					+ "chunks-excluded.txt for details).");
		}

		if (task.getRemovedChunksOverlap() > 0) {
			text.add("\n<b>Warning:</b> " + StringUtils.format(task.getRemovedChunksOverlap())
					+ " Chunk(s) excluded: reference overlap < " + (referenceOverlap * 100) + "% (see "
					+ " chunks-excluded.txt for details).");
		}

		long excludedChunks = task.getRemovedChunksSnps() + task.getRemovedChunksCallRate()
				+ task.getRemovedChunksOverlap();

		long overallChunks = task.getOverallChunks();

		if (excludedChunks > 0) {
			text.add("\nRemaining chunk(s): " + StringUtils.format(overallChunks - excludedChunks));
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

	private HashSet<RangeEntry> parseRangeEntries(String ranges) {
		HashSet<RangeEntry> rangeEntries = new HashSet<RangeEntry>();
		for (String range : ranges.split(",")) {
			String chromosome = range.split(":")[0].trim();
			String region = range.split(":")[1].trim();
			int start = Integer.parseInt(region.split("-")[0].trim());
			int end = Integer.parseInt(region.split("-")[1].trim());
			RangeEntry entry = new RangeEntry();
			entry.setChromosome(chromosome);
			entry.setStart(start);
			entry.setEnd(end);
			rangeEntries.add(entry);
		}
		return rangeEntries;
	}

	protected TaskResults runTask(final OutputWriter output, ITask task) {

		try {
			TaskResults results = task.run();

			if (results.isSuccess()) {
				output.message(task.getName());
			} else {
				output.error(task.getName() + "\n" + results.getMessage());
			}
			return results;
		} catch (Exception e) {
			TaskResults result = new TaskResults();
			result.setSuccess(false);
			result.setMessage(e.getMessage());
			output.error("Task '" + task.getName(), e);
			return result;
		}

	}

	public void setBuild(String build) {
		this.build = build;
	}
}

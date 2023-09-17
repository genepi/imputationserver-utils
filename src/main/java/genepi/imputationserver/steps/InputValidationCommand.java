package genepi.imputationserver.steps;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Callable;

import org.apache.commons.lang.StringEscapeUtils;

import genepi.imputationserver.steps.vcf.VcfFile;
import genepi.imputationserver.steps.vcf.VcfFileUtil;
import genepi.imputationserver.util.RefPanel;
import genepi.imputationserver.util.RefPanelPopulation;
import genepi.imputationserver.util.importer.ImporterFactory;
import genepi.imputationserver.util.report.CloudgeneReport;
import genepi.io.FileUtil;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

@Command
public class InputValidationCommand implements Callable<Integer> {

	@Parameters(description = "VCF files")
	private List<String> files;

	@Option(names = "--population", description = "Reference Population", required = true)
	private String population;

	@Option(names = "--reference", description = "Reference Panel", required = true)
	private String reference;

	@Option(names = "--build", description = "Build", required = false)
	private String build = "hg19";

	@Option(names = "--r2Filter", description = "r2 Filter", required = false)
	private String r2Filter = "0";

	@Option(names = "--phasing", description = "Phasing", required = false)
	private String phasing = "eagle";

	@Option(names = "--mode", description = "Mode", required = false)
	private String mode = "n/a";

	@Option(names = "--chunksize", description = "Chunksize", required = false)
	private int chunksize = 20000000;

	@Option(names = "--minSamples", description = "Min Samples", required = false)
	private int minSamples = 20;

	@Option(names = "--maxSamples", description = "Max Samples", required = false)
	private int maxSamples = 50000;

	@Option(names = "--contactName", description = "Contact Name", required = false)
	private String contactName = "n/a";

	@Option(names = "--contactEmail", description = "Contact Mail", required = false)
	private String contactEmail = "n/a";

	@Option(names = "--report", description = "Cloudgene Report Output", required = false)
	private String report = "cloudgene.report.json";

	private CloudgeneReport context = new CloudgeneReport();

	private RefPanel panel = null;

	public InputValidationCommand() {
		VcfFileUtil.setTabixBinary("tabix");
	}

	public void setFiles(List<String> files) {
		this.files = files;
	}

	public void setReport(String report) {
		this.report = report;
	}

	public void setPhasing(String phasing) {
		this.phasing = phasing;
	}

	public void setPopulation(String population) {
		this.population = population;
	}

	public void setRefPanel(RefPanel panel) {
		this.panel = panel;
	}

	public void setBuild(String build) {
		this.build = build;
	}

	public void setR2Filter(String r2Filter) {
		this.r2Filter = r2Filter;
	}

	public void setMode(String mode) {
		this.mode = mode;
	}

	public void setChunksize(int chunksize) {
		this.chunksize = chunksize;
	}

	public void setMaxSamples(int maxSamples) {
		this.maxSamples = maxSamples;
	}

	public void setMinSamples(int minSamples) {
		this.minSamples = minSamples;
	}

	public void setContactName(String contactName) {
		this.contactName = contactName;
	}

	public void setContactEmail(String contactEmail) {
		this.contactEmail = contactEmail;
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

		if (!checkParameters()) {
			return -1;
		}

		if (!importVcfFiles()) {
			return -1;
		}

		if (!checkVcfFiles()) {
			return -1;
		} else {
			return 0;
		}

	}

	private boolean checkVcfFiles() throws Exception {

		List<VcfFile> validVcfFiles = new Vector<VcfFile>();

		context.beginTask("Analyze files ");
		List<String> chromosomes = new Vector<String>();

		int chunks = 0;
		int noSnps = 0;
		int noSamples = 0;

		boolean phased = true;

		Collections.sort(files);
		String infos = null;

		for (String filename : files) {

			if (infos == null) {
				// first files, no infos available
				context.updateTask(
						"Analyze file " + StringEscapeUtils.escapeHtml(FileUtil.getFilename(filename)) + "...",
						CloudgeneReport.RUNNING);

			} else {
				context.updateTask("Analyze file " + StringEscapeUtils.escapeHtml(FileUtil.getFilename(filename))
						+ "...\n\n" + infos, CloudgeneReport.RUNNING);
			}

			try {

				VcfFile vcfFile = VcfFileUtil.load(filename, chunksize, true);

				if (VcfFileUtil.isChrMT(vcfFile.getChromosome())) {
					vcfFile.setPhased(true);
				}

				if (VcfFileUtil.isValidChromosome(vcfFile.getChromosome())) {

					validVcfFiles.add(vcfFile);
					chromosomes.add(vcfFile.getChromosome());

					String chromosomeString = "";
					for (String chr : chromosomes) {
						chromosomeString += " " + chr;
					}

					// check if all files have same amount of samples
					if (noSamples != 0 && noSamples != vcfFile.getNoSamples()) {
						context.endTask(
								"Please double check, if all uploaded VCF files include the same amount of samples ("
										+ vcfFile.getNoSamples() + " vs " + noSamples + ")",
								CloudgeneReport.ERROR);
						return false;
					}

					noSamples = vcfFile.getNoSamples();
					noSnps += vcfFile.getNoSnps();
					chunks += vcfFile.getChunks().size();

					phased = phased && vcfFile.isPhased();

					if (vcfFile.isPhasedAutodetect() && !vcfFile.isPhased()) {

						context.endTask(
								"File should be phased, but also includes unphased and/or missing genotypes! Please double-check!",
								CloudgeneReport.ERROR);
						return false;
					}

					if (noSamples < minSamples && minSamples != 0) {
						context.endTask("At least " + minSamples + " samples must be uploaded.", CloudgeneReport.ERROR);
						return false;
					}

					if (noSamples > maxSamples && maxSamples != 0) {

						context.endTask("The maximum number of samples is " + maxSamples + ". Please contact "
								+ contactName + " (<a href=\"" + contactEmail + "\">" + contactEmail
								+ "</a>) to discuss this large imputation.", CloudgeneReport.ERROR);
						return false;
					}

					if (build.equals("hg19") && vcfFile.hasChrPrefix()) {
						context.endTask("Your upload data contains chromosome '" + vcfFile.getRawChromosome()
								+ "'. This is not a valid hg19 encoding. Please ensure that your input data is build hg19 and chromosome is encoded as '"
								+ vcfFile.getChromosome() + "'.", CloudgeneReport.ERROR);
						return false;
					}

					if (build.equals("hg38") && !vcfFile.hasChrPrefix()) {
						context.endTask("Your upload data contains chromosome '" + vcfFile.getRawChromosome()
								+ "'. This is not a valid hg38 encoding. Please ensure that your input data is build hg38 and chromosome is encoded as 'chr"
								+ vcfFile.getChromosome() + "'.", CloudgeneReport.ERROR);
						return false;
					}

					infos = "Samples: " + noSamples + "\n" + "Chromosomes:" + chromosomeString + "\n" + "SNPs: "
							+ noSnps + "\n" + "Chunks: " + chunks + "\n" + "Datatype: "
							+ (phased ? "phased" : "unphased") + "\n" + "Build: " + (build == null ? "hg19" : build)
							+ "\n" + "Reference Panel: " + panel.getId() + " (" + panel.getBuild() + ")" + "\n"
							+ "Population: " + population + "\n" + "Phasing: " + phasing + "\n" + "Mode: " + mode;

					if (r2Filter != null && !r2Filter.isEmpty() && !r2Filter.equals("0")) {
						infos += "\nRsq filter: " + r2Filter;
					}

				} else {
					context.endTask("No valid chromosomes found!", CloudgeneReport.ERROR);
					return false;
				}

			} catch (IOException e) {
				context.endTask(StringEscapeUtils.escapeHtml(e.getMessage())
						+ " (see <a href=\"/start.html#!pages/help\">Help</a>).", CloudgeneReport.ERROR);
				return false;

			}

		}

		if (validVcfFiles.size() > 0) {

			context.endTask(validVcfFiles.size() + " valid VCF file(s) found.\n\n" + infos, CloudgeneReport.OK);

			if (!phased && (phasing == null || phasing.isEmpty() || phasing.equals("no_phasing"))) {
				context.error("Your input data is unphased. Please select an algorithm for phasing.");
				return false;
			}

			// init counters
			context.incCounter("samples", noSamples);
			context.incCounter("genotypes", noSamples * noSnps);
			context.incCounter("chromosomes", noSamples * chromosomes.size());
			context.incCounter("runs", 1);
			context.incCounter("refpanel_" + panel.getId(), 1);
			context.incCounter("phasing_" + "eagle", 1);
			return true;

		} else {

			context.endTask("The provided files are not VCF files  (see <a href=\"/start.html#!pages/help\">Help</a>).",
					CloudgeneReport.ERROR);
			return false;
		}
	}

	private boolean checkParameters() throws Exception {

		try {

			if (!panel.supportsPopulation(population)) {
				StringBuilder report = new StringBuilder();
				report.append("Population '" + population + "' is not supported by reference panel '" + panel.getId()
						+ "'.\n");
				if (panel.getPopulations() != null) {
					report.append("Available populations:");
					for (RefPanelPopulation pop : panel.getPopulations()) {
						report.append("\n - " + pop.getId());
					}
				}
				context.error(report.toString());
				return false;
			}

		} catch (Exception e) {
			context.error("Unable to parse reference panel: " + StringEscapeUtils.escapeHtml(e.getMessage()));
			return false;
		}

		return true;
	}

	private boolean importVcfFiles() throws Exception {

		for (String input : files) {

			if (ImporterFactory.needsImport(input)) {

				context.log("URL-based uploads are no longer supported. Please use direct file uploads instead.");
				context.error("URL-based uploads are no longer supported. Please use direct file uploads instead.");
				return false;
			}

		}

		return true;

	}

}

package genepi.imputationserver.steps;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.concurrent.Callable;

import org.apache.commons.lang.StringEscapeUtils;

import cloudgene.sdk.internal.WorkflowContext;
import cloudgene.sdk.weblog.WebWorkflowContext;
import cloudgene.sdk.weblog.collectors.FileLog;
import genepi.imputationserver.steps.vcf.VcfFile;
import genepi.imputationserver.steps.vcf.VcfFileUtil;
import genepi.imputationserver.util.RefPanel;
import genepi.imputationserver.util.importer.ImporterFactory;
import genepi.io.FileUtil;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

@Command
public class InputValidationCommand implements Callable<Integer> {

	public InputValidationCommand() {
		VcfFileUtil.setTabixBinary("tabix");
	}

	@Parameters(description = "VCF files")
	List<String> files;

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

	@Option(names = "--maxSamples", description = "Max Samples", required = false)
	private int maxSamples = 50000;

	@Option(names = "--contactName", description = "Contact Name", required = false)
	private String contactName = "n/a";

	@Option(names = "--contactEmail", description = "Contact Mail", required = false)
	private String contactEmail = "n/a";

	@Option(names = "--output", description = "Log Output", required = false)
	private String output = "cloudgene.log";

	WebWorkflowContext context = new WebWorkflowContext();

	RefPanel panel = null;

	public void setFiles(List<String> files) {
		this.files = files;
	}

	public void setOutput(String output) {
		this.output = output;
	}

	public void setPhasing(String phasing) {
		this.phasing = phasing;
	}

	public void setPopulation(String population) {
		this.population = population;
	}

	public void setReference(Map<String, Object> reference) {

		try {
			panel = RefPanel.loadFromProperties(reference);
			if (panel == null) {
				context.error("Reference not found.");
			}
		} catch (Exception e) {
			context.error("Unable to parse reference panel: " + StringEscapeUtils.escapeHtml(e.getMessage()));
		}
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

	public void setContactName(String contactName) {
		this.contactName = contactName;
	}

	public void setContactEmail(String contactEmail) {
		this.contactEmail = contactEmail;
	}

	@Override
	public Integer call() throws Exception {

		context.setCollector(new FileLog(output));

		if (panel == null) {
			try {
				panel = RefPanel.loadFromJson(reference);
				if (panel == null) {
					context.error("Reference not found.");
					context.close();
					return -1;
				}
			} catch (Exception e) {
				context.error("Unable to parse reference panel: " + StringEscapeUtils.escapeHtml(e.getMessage()));
				context.close();
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

	protected void setupTabix(String path) {
		VcfFileUtil.setTabixBinary(path);
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
						WorkflowContext.RUNNING);

			} else {
				context.updateTask("Analyze file " + StringEscapeUtils.escapeHtml(FileUtil.getFilename(filename))
						+ "...\n\n" + infos, WorkflowContext.RUNNING);
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
								WorkflowContext.ERROR);
						context.close();
						return false;
					}

					noSamples = vcfFile.getNoSamples();
					noSnps += vcfFile.getNoSnps();
					chunks += vcfFile.getChunks().size();

					phased = phased && vcfFile.isPhased();

					if (vcfFile.isPhasedAutodetect() && !vcfFile.isPhased()) {

						context.endTask(
								"File should be phased, but also includes unphased and/or missing genotypes! Please double-check!",
								WorkflowContext.ERROR);
						context.close();
						return false;
					}

					if (noSamples > maxSamples && maxSamples != 0) {

						context.endTask("The maximum number of samples is " + maxSamples + ". Please contact "
								+ contactName + " (<a href=\"" + contactEmail + "\">" + contactEmail
								+ "</a>) to discuss this large imputation.", WorkflowContext.ERROR);
						context.close();
						return false;
					}

					if (build.equals("hg19") && vcfFile.hasChrPrefix()) {
						context.endTask("Your upload data contains chromosome '" + vcfFile.getRawChromosome()
								+ "'. This is not a valid hg19 encoding. Please ensure that your input data is build hg19 and chromosome is encoded as '"
								+ vcfFile.getChromosome() + "'.", WorkflowContext.ERROR);
						context.close();
						return false;
					}

					if (build.equals("hg38") && !vcfFile.hasChrPrefix()) {
						context.endTask("Your upload data contains chromosome '" + vcfFile.getRawChromosome()
								+ "'. This is not a valid hg38 encoding. Please ensure that your input data is build hg38 and chromosome is encoded as 'chr"
								+ vcfFile.getChromosome() + "'.", WorkflowContext.ERROR);
						context.close();
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
					context.endTask("No valid chromosomes found!", WorkflowContext.ERROR);
					context.close();
					return false;
				}

			} catch (IOException e) {
				context.endTask(StringEscapeUtils.escapeHtml(e.getMessage())
						+ " (see <a href=\"/start.html#!pages/help\">Help</a>).", WorkflowContext.ERROR);
				context.close();
				return false;

			}

		}

		if (validVcfFiles.size() > 0) {

			context.endTask(validVcfFiles.size() + " valid VCF file(s) found.\n\n" + infos, WorkflowContext.OK);

			if (!phased && (phasing == null || phasing.isEmpty() || phasing.equals("no_phasing"))) {
				context.error("Your input data is unphased. Please select an algorithm for phasing.");
				context.close();
				return false;
			}

			// init counters
			context.incCounter("samples", noSamples);
			context.incCounter("genotypes", noSamples * noSnps);
			context.incCounter("chromosomes", noSamples * chromosomes.size());
			context.incCounter("runs", 1);
			context.incCounter("refpanel_" + panel.getId(), 1);
			context.incCounter("phasing_" + "eagle", 1);
			context.close();
			return true;

		} else {

			context.endTask("The provided files are not VCF files  (see <a href=\"/start.html#!pages/help\">Help</a>).",
					WorkflowContext.ERROR);
			context.close();
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
					for (String pop : panel.getPopulations().values()) {
						report.append("\n - " + pop);
					}
				}
				context.error(report.toString());
				context.close();
				return false;
			}

		} catch (Exception e) {
			context.error("Unable to parse reference panel: " + StringEscapeUtils.escapeHtml(e.getMessage()));
			context.close();
			return false;
		}

		return true;
	}

	private boolean importVcfFiles() throws Exception {

		for (String input : files) {

			if (ImporterFactory.needsImport(input)) {

				context.log("URL-based uploads are no longer supported. Please use direct file uploads instead.");
				context.error("URL-based uploads are no longer supported. Please use direct file uploads instead.");
				context.close();
				return false;
			}

		}

		return true;

	}

	public String getFolder(Class clazz) {
		return new File(clazz.getProtectionDomain().getCodeSource().getLocation().getPath()).getParent();
	}

}

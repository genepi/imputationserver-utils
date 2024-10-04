package genepi.imputationserver.steps;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Callable;

import genepi.imputationserver.util.OutputWriter;

import genepi.imputationserver.steps.vcf.VcfFile;
import genepi.imputationserver.steps.vcf.VcfFileUtil;
import genepi.imputationserver.util.RefPanel;
import genepi.imputationserver.util.RefPanelPopulation;
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

	@Option(names = "--no-index", description = "Create no tabix index during validation", required = false)
	private boolean noIndex = false;

	@Option(names = "--report", description = "Cloudgene Report Output", required = false)
	private String report = null;

	private RefPanel panel = null;

	private OutputWriter output = null;


	public InputValidationCommand() {

	}

	public void setFiles(List<String> files) {
		this.files = files;
	}

	public void setReport(String report) {
		this.report = report;
	}

	public void setReference(String reference) {
		this.reference = reference;
	}

	public void setPhasing(String phasing) {
		this.phasing = phasing;
	}

	public void setPopulation(String population) {
		this.population = population;
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
			output.error("Unable to parse reference panel:" ,e);
			return 1;
		}

		if (!checkParameters()) {
			return 1;
		}

		if (!checkVcfFiles()) {
			return 1;
		} else {
			return 0;
		}

	}

	private boolean checkVcfFiles() throws Exception {

		List<VcfFile> validVcfFiles = new Vector<VcfFile>();
		List<String> chromosomes = new Vector<String>();

		int chunks = 0;
		int noSnps = 0;
		int noSamples = 0;

		boolean phased = true;

		Collections.sort(files);

		for (String filename : files) {

			try {

				VcfFile vcfFile = VcfFileUtil.load(filename, chunksize, !noIndex);

				if (VcfFileUtil.isChrMT(vcfFile.getChromosome())) {
					vcfFile.setPhased(true);
				}

				if (!VcfFileUtil.isValidChromosome(vcfFile.getChromosome())) {
					output.error("No valid chromosomes found!");
					return false;
				}

				validVcfFiles.add(vcfFile);
				chromosomes.add(vcfFile.getChromosome());

				// check if all files have same amount of samples
				if (noSamples != 0 && noSamples != vcfFile.getNoSamples()) {
					output.error("Please double check, if all uploaded VCF files include the same amount of samples ("
									+ vcfFile.getNoSamples() + " vs " + noSamples + ")");
					return false;
				}

				noSamples = vcfFile.getNoSamples();
				noSnps += vcfFile.getNoSnps();
				chunks += vcfFile.getChunks().size();

				phased = phased && vcfFile.isPhased();

				if (vcfFile.isPhasedAutodetect() && !vcfFile.isPhased()) {
					output.error("File should be phased, but also includes unphased and/or missing genotypes! Please double-check!");
					return false;
				}

				if (noSamples < minSamples && minSamples != 0) {
					output.error("At least " + minSamples + " samples must be uploaded.");
					return false;
				}

				if (noSamples > maxSamples && maxSamples != 0) {

					output.error("The maximum number of samples is " + maxSamples + ". Please contact "
							+ contactName + " (<a href=\"" + contactEmail + "\">" + contactEmail
							+ "</a>) to discuss this large imputation.");
					return false;
				}

				if (build.equals("hg19") && vcfFile.hasChrPrefix()) {
					output.error("Your upload data contains chromosome '" + vcfFile.getRawChromosome()
							+ "'. This is not a valid hg19 encoding. Please ensure that your input data is build hg19 and chromosome is encoded as '"
							+ vcfFile.getChromosome() + "'.");
					return false;
				}

				if (build.equals("hg38") && !vcfFile.hasChrPrefix()) {
					output.error("Your upload data contains chromosome '" + vcfFile.getRawChromosome()
							+ "'. This is not a valid hg38 encoding. Please ensure that your input data is build hg38 and chromosome is encoded as 'chr"
							+ vcfFile.getChromosome() + "'.");
					return false;
				}


			} catch (IOException e) {
				output.error(e);
				return false;
			}

		}

		if (validVcfFiles.isEmpty()) {
			output.error("The provided files are not VCF files (see <a href=\"/start.html#!pages/help\">Help</a>).");
			return false;
		}

		if (!phased && (phasing == null || phasing.isEmpty() || phasing.equals("no_phasing"))) {
			output.error("Your input data is unphased. Please select an algorithm for phasing.");
			return false;
		}

		List<String> summary = new Vector<String>();
		summary.add(validVcfFiles.size() + " valid VCF file(s) found.");
		summary.add("");
		summary.add("Samples: " + noSamples);
		summary.add("Chromosomes: " + String.join(" ", chromosomes));
		summary.add("SNPs: " + noSnps);
		summary.add("Chunks: " + chunks);
		summary.add("Datatype: " + (phased ? "phased" : "unphased"));
		summary.add("Build: " + (build == null ? "hg19" : build));
		summary.add("Reference Panel: " + panel.getId() + " (" + panel.getBuild() + ")" );
		summary.add("Population: " + population);
		summary.add("Phasing: " + phasing);
		summary.add("Mode: " + mode);
		if (r2Filter != null && !r2Filter.isEmpty() && !r2Filter.equals("0")) {
			summary.add("Rsq filter: " + r2Filter);
		}
		output.message(summary);

		// init counters
		output.print("");
		output.setCounter("samples", noSamples);
		output.setCounter("variants",  noSnps);
		output.setCounter("chromosomes", noSamples * chromosomes.size());
		output.setCounter("runs", 1);
		return true;

	}

	private boolean checkParameters() throws Exception {

		try {

			if (!panel.supportsPopulation(population)) {
				List<String> messages = new Vector<String>();
				messages.add("Population '" + population + "' is not supported by reference panel '" + panel.getId() + "'.");
				if (panel.getPopulations() != null) {
					messages.add("Available populations:");
					for (RefPanelPopulation pop : panel.getPopulations()) {
						messages.add(" - " + pop.getId());
					}
				}
				output.error(messages);
				return false;
			}

		} catch (Exception e) {
			output.error("Unable to parse reference panel. ", e);
			return false;
		}

		return true;
	}

	public void setNoIndex(boolean noIndex) {
		this.noIndex = noIndex;
	}
}

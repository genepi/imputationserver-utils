package genepi.imputationserver.steps;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import genepi.imputationserver.steps.ancestry.TraceInputValidation;
import genepi.imputationserver.steps.ancestry.TraceInputValidation.TraceInputValidationResult;
import genepi.imputationserver.steps.fastqc.ITask;
import genepi.imputationserver.steps.fastqc.ITaskProgressListener;
import genepi.imputationserver.steps.fastqc.LiftOverTask;
import genepi.imputationserver.steps.fastqc.TaskResults;
import genepi.imputationserver.steps.vcf.VcfFileUtil;
import genepi.imputationserver.util.report.CloudgeneReport;
import genepi.io.FileUtil;
import genepi.io.text.LineWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

@Command
public class PrepareTraceCommand implements Callable<Integer> {

	public static String BUILD = "hg19";

	public static final int MIN_VARIANTS = 100;

	@Parameters(description = "VCF files")
	private String[] files;

	@Option(names = "--build", description = "Reference Population", required = true)
	private String buildGwas;

	@Option(names = "--batch-size", description = "Reference Population", required = true)
	private int batchSize;

	@Option(names = "--output", description = "Output", required = true)
	private String output;

	@Option(names = "--report", description = "Report file", required = true)
	private String reportFilename;

	@Option(names = "--reference-sites", description = "Reference Sites", required = true)
	private String referenceSites;

	public PrepareTraceCommand() {
		VcfFileUtil.setTabixBinary("tabix");
	}

	protected void setupTabix(String path) {
		VcfFileUtil.setTabixBinary(path);
	}

	@Override
	public Integer call() throws Exception {
		return prepareTraceJobs() ? 0 : -1;
	}

	public boolean prepareTraceJobs() {

		CloudgeneReport report = new CloudgeneReport();
		report.setFilename(reportFilename);

		try {

			// Check if liftover is needed
			if (!buildGwas.equals(BUILD)) {

				report.warning("Uploaded data is " + buildGwas + " and reference is " + BUILD + ".");
				String chainFile = "";//TODO!!
				/*if (chainFile == null) {
					report.error("Currently we do not support liftOver from " + buildGwas + " to " + BUILD);
					return false;
				}*/

				String fullPathChainFile = "";
				if (!new File(fullPathChainFile).exists()) {
					report.error("Chain file " + fullPathChainFile + " not found.");
					return false;
				}

				String chunksDir = FileUtil.path("liftover");
				FileUtil.createDirectory(chunksDir);

				LiftOverTask task = new LiftOverTask();
				task.setVcfFilenames(files);
				task.setChainFile(fullPathChainFile);
				task.setChunksDir(chunksDir);
				task.setExcludedSnpsWriter(null);

				TaskResults results = runTask(report, task);

				if (results.isSuccess()) {
					files = task.getNewVcfFilenames();
				} else {
					return false;
				}

			}

			String mergedFile = FileUtil.path(output, "study.merged.vcf.gz");

			if (!checkDataAndMerge(report, files, mergedFile)) {
				return false;
			}

			// read number of samples from first vcf file
			VCFFileReader reader = new VCFFileReader(new File(mergedFile), false);
			VCFHeader header = reader.getFileHeader();
			reader.close();
			ArrayList<String> sampels = header.getSampleNamesInOrder();
			int nIndividuals = sampels.size();
			int batch = 0;
			int start = 1;
			int end;

			while (start <= nIndividuals) {
				end = start + batchSize - 1;
				if (end > nIndividuals) {
					end = nIndividuals;
				}

				LineWriter lineWriter = new LineWriter(FileUtil.path(output, batch + ".batch"));
				for (int i = start - 1; i <= end - 1; i++) {
					lineWriter.write(sampels.get(i));
				}
				lineWriter.close();

				start = end + 1;
				batch++;

				report.log("Created batch No. " + batch);
			}

			report.ok("Prepared " + batch + " batch job" + ((batch > 1) ? "s." : "."));

			return true;

		} catch (IOException e) {
			report.error("An internal server error occured.");
			e.printStackTrace();
		}

		report.error("Execution failed. Please, contact administrator.");

		return false;
	}

	public boolean checkDataAndMerge(CloudgeneReport report, String[] files, String mergedFile) {

		try {

			TraceInputValidation validation = new TraceInputValidation();

			TraceInputValidationResult result = validation.mergeAndCheckSites(files, referenceSites, mergedFile);

			String message = "Loaded " + result.getTotal() + " variants" + "\n" + "Variants with different alleles: "
					+ result.getAlleleMissmatch() + "\n" + "Variants with allele switches: "
					+ result.getAlleleSwitches() + "\n" + "Variants not found in reference: " + result.getNotFound()
					+ "\n" + "Overlapping variants used by LASER: " + result.getFound();

			report.ok(message);

			if (result.getFound() <= MIN_VARIANTS) {
				report.error("Number of variants shared with reference is too small (&le;" + MIN_VARIANTS
						+ ").\nPlease, check if input data are correct or try to use another ancestry reference panel.");
				return false;
			}

			return true;

		} catch (IOException e) {
			report.error("Input Validation failed: " + e);
			e.printStackTrace();
			return false;
		}
	}

	protected TaskResults runTask(final CloudgeneReport report, ITask task) {
		report.beginTask("Running " + task.getName() + "...");
		TaskResults results;
		try {
			results = task.run(new ITaskProgressListener() {

				@Override
				public void progress(String message) {
					report.updateTask(message, CloudgeneReport.RUNNING);

				}
			});

			if (results.isSuccess()) {
				report.endTask(task.getName(), CloudgeneReport.OK);
			} else {
				report.endTask(task.getName() + "\n" + results.getMessage(), CloudgeneReport.ERROR);
			}
			return results;
		} catch (Exception | Error e) {
			e.printStackTrace();
			TaskResults result = new TaskResults();
			result.setSuccess(false);
			result.setMessage(e.getMessage());
			StringWriter s = new StringWriter();
			e.printStackTrace(new PrintWriter(s));
			report.println("Task '" + task.getName() + "' failed.\nException:" + s.toString());
			report.endTask(e.getMessage(), CloudgeneReport.ERROR);
			return result;
		}

	}

}

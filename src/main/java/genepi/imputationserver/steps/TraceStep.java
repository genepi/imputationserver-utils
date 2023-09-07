package genepi.imputationserver.steps;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;

import cloudgene.sdk.internal.WorkflowContext;
import cloudgene.sdk.internal.WorkflowStep;
import genepi.imputationserver.steps.ancestry.TraceInputValidation;
import genepi.imputationserver.steps.ancestry.TraceInputValidation.TraceInputValidationResult;
import genepi.imputationserver.steps.fastqc.ITask;
import genepi.imputationserver.steps.fastqc.ITaskProgressListener;
import genepi.imputationserver.steps.fastqc.LiftOverTask;
import genepi.imputationserver.steps.fastqc.TaskResults;
import genepi.imputationserver.steps.vcf.VcfFileUtil;
import genepi.imputationserver.util.DefaultPreferenceStore;
import genepi.io.FileUtil;
import genepi.io.text.LineWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class TraceStep extends WorkflowStep {

	public static String BUILD = "hg19";

	public static final int MIN_VARIANTS = 100;

	@Override
	public boolean run(WorkflowContext context) {
		prepareTraceJobs(context);
		return true;
	}

	protected void setupTabix(String folder) {
		VcfFileUtil.setTabixBinary(FileUtil.path(folder, "bin", "tabix"));
	}

	public boolean prepareTraceJobs(WorkflowContext context) {

		String folder = getFolder(QualityControlCommand.class);
		setupTabix(folder);

		String genotypes = context.get("files");
		String buildGwas = context.get("build");

		int batchSize = Integer.parseInt(context.get("batch_size"));
		String output = context.get("output");

		String[] files = FileUtil.getFiles(genotypes, "*.vcf.gz$|*.vcf$");

		try {

			// load job.config
			File jobConfig = new File(FileUtil.path(folder, "job.config"));
			DefaultPreferenceStore store = new DefaultPreferenceStore();
			if (jobConfig.exists()) {
				store.load(jobConfig);
			} else {
				context.log(
						"Configuration file '" + jobConfig.getAbsolutePath() + "' not available. Use default values.");
			}

			// Check if liftover is needed
			if (!buildGwas.equals(BUILD)) {

				context.warning("Uploaded data is " + buildGwas + " and reference is " + BUILD + ".");
				String chainFile = store.getString(buildGwas + "To" + BUILD);
				if (chainFile == null) {
					context.error("Currently we do not support liftOver from " + buildGwas + " to " + BUILD);
					return false;
				}

				String fullPathChainFile = FileUtil.path(folder, chainFile);
				if (!new File(fullPathChainFile).exists()) {
					context.error("Chain file " + fullPathChainFile + " not found.");
					return false;
				}

				String chunksDir = FileUtil.path(context.getLocalTemp(), "liftover");
				FileUtil.createDirectory(chunksDir);

				LiftOverTask task = new LiftOverTask();
				task.setVcfFilenames(files);
				task.setChainFile(fullPathChainFile);
				task.setChunksDir(chunksDir);
				task.setExcludedSnpsWriter(null);

				TaskResults results = runTask(context, task);

				if (results.isSuccess()) {
					files = task.getNewVcfFilenames();
				} else {
					return false;
				}

			}

			String mergedFile = FileUtil.path(output, "study.merged.vcf.gz");

			if (!checkDataAndMerge(context, files, mergedFile)) {
				return false;
			}

			context.beginTask("Preparing TRACE jobs");

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

				context.log("\t- Created batch No. " + batch);
			}

			context.endTask("Prepared " + batch + " batch job" + ((batch > 1) ? "s." : "."), WorkflowContext.OK);

			return true;

		} catch (IOException e) {
			context.error("An internal server error occured.");
			e.printStackTrace();
		}

		context.error("Execution failed. Please, contact administrator.");

		return false;
	}

	public boolean checkDataAndMerge(WorkflowContext context, String[] files, String mergedFile) {

		try {

			context.beginTask("Input Validation");

			String referenceSites = context.get("reference_sites");

			TraceInputValidation validation = new TraceInputValidation();

			TraceInputValidationResult result = validation.mergeAndCheckSites(files, referenceSites, mergedFile);

			String message = "Loaded " + result.getTotal() + " variants" + "\n" + "Variants with different alleles: "
					+ result.getAlleleMissmatch() + "\n" + "Variants with allele switches: "
					+ result.getAlleleSwitches() + "\n" + "Variants not found in reference: " + result.getNotFound()
					+ "\n" + "Overlapping variants used by LASER: " + result.getFound();

			context.endTask(message, WorkflowContext.OK);

			if (result.getFound() <= MIN_VARIANTS) {
				context.error("Number of variants shared with reference is too small (&le;" + MIN_VARIANTS
						+ ").\nPlease, check if input data are correct or try to use another ancestry reference panel.");
				return false;
			}

			return true;

		} catch (IOException e) {
			context.error("Input Validation failed: " + e);
			e.printStackTrace();
			return false;
		}
	}

	protected TaskResults runTask(final WorkflowContext context, ITask task) {
		context.beginTask("Running " + task.getName() + "...");
		TaskResults results;
		try {
			results = task.run(new ITaskProgressListener() {

				@Override
				public void progress(String message) {
					context.updateTask(message, WorkflowContext.RUNNING);

				}
			});

			if (results.isSuccess()) {
				context.endTask(task.getName(), WorkflowContext.OK);
			} else {
				context.endTask(task.getName() + "\n" + results.getMessage(), WorkflowContext.ERROR);
			}
			return results;
		} catch (InterruptedException e) {
			e.printStackTrace();
			TaskResults result = new TaskResults();
			result.setSuccess(false);
			result.setMessage(e.getMessage());
			StringWriter s = new StringWriter();
			e.printStackTrace(new PrintWriter(s));
			context.println("Task '" + task.getName() + "' failed.\nException:" + s.toString());
			context.endTask(e.getMessage(), WorkflowContext.ERROR);
			return result;
		} catch (Exception e) {
			e.printStackTrace();
			TaskResults result = new TaskResults();
			result.setSuccess(false);
			result.setMessage(e.getMessage());
			StringWriter s = new StringWriter();
			e.printStackTrace(new PrintWriter(s));
			context.println("Task '" + task.getName() + "' failed.\nException:" + s.toString());
			context.endTask(task.getName() + " failed.", WorkflowContext.ERROR);
			return result;
		} catch (Error e) {
			e.printStackTrace();
			TaskResults result = new TaskResults();
			result.setSuccess(false);
			result.setMessage(e.getMessage());
			StringWriter s = new StringWriter();
			e.printStackTrace(new PrintWriter(s));
			context.println("Task '" + task.getName() + "' failed.\nException:" + s.toString());
			context.endTask(task.getName() + " failed.", WorkflowContext.ERROR);
			return result;
		}

	}

}

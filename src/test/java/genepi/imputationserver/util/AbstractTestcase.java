package genepi.imputationserver.util;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import cloudgene.sdk.internal.WorkflowStep;
import genepi.imputationserver.steps.CompressionEncryption;
import genepi.imputationserver.steps.FastQualityControl;
import genepi.imputationserver.steps.InputValidationCommand;
import genepi.imputationserver.steps.vcf.VcfFileUtil;
import genepi.io.FileUtil;
import genepi.io.text.LineReader;

public class AbstractTestcase {

	public static final boolean VERBOSE = true;

	public static final String BINARIES_HDFS = "binaries";

	public static final String PASSWORD = "random-pwd";
	
	protected boolean run(WorkflowTestContext context, WorkflowStep step) {
		step.setup(context);
		return step.run(context);
	}

	protected WorkflowTestContext buildContext(String folder, String configFolder, String refpanel) throws IOException {
		WorkflowTestContext context = new WorkflowTestContext();

		File file = new File("test-data/tmp");
		if (file.exists()) {
			FileUtil.deleteDirectory(file);
		}
		file.mkdirs();
		
		context.setInput("project", "test-job");

		context.setVerbose(VERBOSE);
		context.setInput("files", folder);
		context.setInput("population", "eur");
		context.setInput("password", PASSWORD);
		context.setInput("phasing", "eagle");

		// load reference panel form file
		context.setInput("refpanel", "apps@" + refpanel);
		Map<String, Object> panel = RefPanelUtil.loadFromFile(configFolder + "/panels.txt", refpanel);
		resolveEnvVariable(panel, configFolder);
		context.setData("refpanel", panel);


		
		context.setOutput("mafFile", file.getAbsolutePath() + "/mafFile/mafFile.txt");
		FileUtil.createDirectory(file.getAbsolutePath() + "/mafFile");

		context.setOutput("chunkFileDir", file.getAbsolutePath() + "/chunkFileDir");
		FileUtil.createDirectory(file.getAbsolutePath() + "/chunkFileDir");

		context.setOutput("statisticDir", file.getAbsolutePath() + "/statisticDir");
		FileUtil.createDirectory(file.getAbsolutePath() + "/statisticDir");

		context.setOutput("chunksDir", file.getAbsolutePath() + "/chunksDir");
		FileUtil.createDirectory(file.getAbsolutePath() + "/chunksDir");

		context.setOutput("local", file.getAbsolutePath() + "/local");
		FileUtil.createDirectory(file.getAbsolutePath() + "/local");

		context.setOutput("logfile", file.getAbsolutePath() + "/logfile");
		FileUtil.createDirectory(file.getAbsolutePath() + "/logfile");

		context.setOutput("outputimputation", "cloudgene-hdfs");

		context.setOutput("hadooplogs", file.getAbsolutePath() + "/hadooplogs");
		FileUtil.deleteDirectory(file.getAbsolutePath() + "/hadooplogs");
		FileUtil.createDirectory(file.getAbsolutePath() + "/hadooplogs");

		context.setLocalTemp("local-temp");
		FileUtil.deleteDirectory("local-temp");
		FileUtil.createDirectory("local-temp");

		return context;

	}

	public class CompressionEncryptionMock extends CompressionEncryption {

		private String folder;

		public CompressionEncryptionMock(String folder) {
			super();
			this.folder = folder;
		}

		@Override
		public String getFolder(Class clazz) {
			// override folder with static folder instead of jar location
			return folder;
		}

	}

	public class QcStatisticsMock extends FastQualityControl {

		private String folder;

		public QcStatisticsMock(String folder) {
			super();
			this.folder = folder;
		}

		@Override
		public String getFolder(Class clazz) {
			// override folder with static folder instead of jar location
			return folder;
		}

		@Override
		protected void setupTabix(String folder) {
			VcfFileUtil.setTabixBinary("files/bin/tabix");
		}

	}

	public class InputValidationMock extends InputValidationCommand {

		private String folder;

		public InputValidationMock(String folder) {
			super();
			this.folder = folder;
		}

		@Override
		public String getFolder(Class clazz) {
			// override folder with static folder instead of jar location
			return folder;
		}

		@Override
		protected void setupTabix(String folder) {
			VcfFileUtil.setTabixBinary("files/bin/tabix");
		}

	}
	
	protected int getLineCount(String filename) throws IOException {
		LineReader reader = new LineReader(filename);
		int lines = 0;
		while (reader.next()) {
			lines++;
		}
		return lines;
	}

	protected boolean checkAmountOfColumns(String filename, int tabs) throws IOException {
		LineReader reader = new LineReader(filename);
		while (reader.next()) {

			String line = reader.get();

			if (line.split("\t").length > tabs) {
				return false;
			}

		}

		return true;
	}

	protected boolean checkSortPositionInfo(String filename) throws IOException {

		LineReader reader = new LineReader(filename);
		int pos = -1;
		while (reader.next()) {

			String line = reader.get();

			if (!line.startsWith("SNP")) {
				String snp = line.split("\t")[0];
				if (Integer.valueOf(snp.split(":")[1]) <= pos) {
					return false;
				}
				pos = Integer.valueOf(snp.split(":")[1]);
			}

		}

		return true;
	}

	
	public void resolveEnvVariable(Map<String, Object> properties, String folder) {
		for (String key: properties.keySet()) {
			Object value = properties.get(key);
			if (value instanceof String) {
				String valueString = value.toString().replaceAll("\\$\\{app_local_folder\\}", folder);
				properties.put(key, valueString);
			}
		}
	}
	
	
}

package genepi.imputationserver.steps.fastqc;

import java.io.IOException;
import java.util.Vector;

import genepi.imputationserver.steps.fastqc.io.ExcludedSnpsWriter;
import genepi.imputationserver.steps.vcf.VcfFileUtil;
import genepi.imputationserver.steps.vcf.VcfLiftOverFast;
import genepi.io.FileUtil;
import genepi.io.text.LineWriter;

public class LiftOverTask implements ITask {

	public static final String X_PAR1 = "X.PAR1";
	public static final String X_PAR2 = "X.PAR2";
	public static final String X_NON_PAR = "X.nonPAR";

	// input variables
	private String chainFile;
	private String chunksDir;
	private String[] vcfFilenames;
	private String[] newVcfFilenames;

	private String statDir;

	@Override
	public String getName() {
		return "Lift Over";
	}

	public TaskResults run() throws IOException {

		String excludedSnpsFile = FileUtil.path(statDir, "lift-over.txt");
		LineWriter excludedSnpsWriter = new LineWriter(excludedSnpsFile);

		newVcfFilenames = new String[vcfFilenames.length];
		for (int i = 0; i < vcfFilenames.length; i++) {
			String filename = vcfFilenames[i];
			String name = FileUtil.getFilename(filename);
			String output = FileUtil.path(chunksDir, name + ".lifted.vcf.gz");
			String temp = FileUtil.path(chunksDir, "vcf.sorted");
			FileUtil.createDirectory(temp);
			Vector<String> errors = VcfLiftOverFast.liftOver(filename, output, chainFile, temp);
			
			// create tabix index
			VcfFileUtil.createIndex(output, true);
			
			FileUtil.deleteDirectory(temp);
			for (String error : errors) {
				excludedSnpsWriter.write(error);
			}
			newVcfFilenames[i] = output;
		}

		TaskResults result = new TaskResults();
		result.setMessage("");
		result.setSuccess(true);

		excludedSnpsWriter.close();

		return result;

	}

	public void setVcfFilenames(String[] vcfFilenames) {
		this.vcfFilenames = vcfFilenames;
	}

	public void setChunksDir(String chunksDir) {
		this.chunksDir = chunksDir;
	}

	public void setChainFile(String chainFile) {
		this.chainFile = chainFile;
	}

	public String[] getNewVcfFilenames() {
		return newVcfFilenames;
	}

	public void setStatDir(String statDir) {
		this.statDir = statDir;
	}
}

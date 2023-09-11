package genepi.imputationserver.util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;

import genepi.io.text.LineReader;

public class AbstractTestcase {

	public static final boolean VERBOSE = true;

	public static final String PASSWORD = "random-pwd";
	
	protected ArrayList<String> getFiles(String inputFolder) {
		File folder = new File(inputFolder);
		ArrayList<String> files = new ArrayList<String>();
		
		for (File file : folder.listFiles()){
			if(file.getName().endsWith("vcf.gz")) {
			files.add(file.getAbsolutePath());
			}
		}
		return files;
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

package genepi.imputationserver.steps.vcf;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.List;

import org.apache.commons.io.IOUtils;

import genepi.imputationserver.util.report.CloudgeneReport;
import genepi.io.FileUtil;
import genepi.io.text.LineReader;
import htsjdk.samtools.util.BlockCompressedStreamConstants;

public class MergedVcfFile {

	private FileOutputStream output;

	public MergedVcfFile(String filename) throws FileNotFoundException {
		output = new FileOutputStream(filename);
	}

	public void addFile(InputStream input) throws IOException {
		IOUtils.copy(input, output);
		input.close();
	}

	public void close() throws IOException {
		output.write(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK);
		output.close();
	}

	public void addHeader(CloudgeneReport report, List<String> files) throws Exception {
		// simple header check
		String headerLine = null;
		for (String file : files) {

			report.println("Read header file " + file);
			LineReader reader = null;
			try {
				reader = new LineReader(file);
				while (reader.next()) {
					String line = reader.get();
					if (line.startsWith("#CHROM")) {
						if (headerLine != null) {
							if (headerLine.equals(line)) {
								report.println("  Header is the same as header of first file.");
							} else {
								report.println("  ERROR: Header is different as header of first file.");
								report.println(headerLine);
								report.println(line);
								throw new Exception("Different sample order in chunks.");
							}
						} else {
							headerLine = line;
							addFile(new FileInputStream(file));
							report.println("  Keep this header as first header.");
						}
					}

				}
				if (reader != null) {
					reader.close();
				}
				FileUtil.deleteFile(file);
			} catch (Exception e) {
				if (reader != null) {
					reader.close();
				}
				StringWriter errors = new StringWriter();
				e.printStackTrace(new PrintWriter(errors));
				report.println("Error reading header file: " + errors.toString());
			}
		}

		if (headerLine == null || headerLine.trim().isEmpty()) {
			throw new Exception("No valid header file found");
		}
	}

}

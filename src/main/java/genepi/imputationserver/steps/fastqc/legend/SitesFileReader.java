package genepi.imputationserver.steps.fastqc.legend;

import genepi.io.FileUtil;
import genepi.io.text.LineReader;
import htsjdk.tribble.readers.TabixReader;

import java.io.*;

public class SitesFileReader {

	private SitesEntry entry = new SitesEntry();

	private String population;

	private TabixReader tabixReader;

	private int idColumn = -1;

	private int refColumn = -1;

	private int altColumn = -1;

	private int popColumn = -1;

	public static String COLUMN_ID = "ID";

	public static String COLUMN_AAF_PREFIX = "AAF_";

	public static String COLUMN_REF = "REF";

	public static String COLUMN_ALT = "ALT";


	public SitesFileReader(String filename, String population) throws IOException {
		this.population = population;
		if (!new File(filename).exists()) {
			throw new IOException("File '" + filename + "' not found.");
		}

		LineReader reader = new LineReader(openTxtOrGzipStream(filename));
		if (!reader.next()) {
			throw new IOException("File '" + filename + "' is empty.");
		}

		String header = reader.get();
		parseHeader(header);
		reader.close();

		tabixReader = new TabixReader(filename);
	}

	protected void parseHeader(String line) throws IOException {
		// parse header
		String[] tiles = line.split("\t");
		int i = 0;
		for (String tile : tiles) {
			if (tile.equalsIgnoreCase(COLUMN_ID)) {
				idColumn = i;
			}
			if (tile.equalsIgnoreCase(COLUMN_REF)) {
				refColumn = i;
			}
			if (tile.equalsIgnoreCase(COLUMN_ALT)) {
				altColumn = i;
			}
			if (tile.equalsIgnoreCase(COLUMN_AAF_PREFIX + population)) {
				popColumn = i;
			}
			i++;
		}

		// Validation: Check if all required columns have been found
		if (idColumn == -1) {
			throw new IOException("Column '" + COLUMN_ID + "' not found in file.");
		}
		if (refColumn == -1) {
			throw new IOException("Column '" + COLUMN_REF + "' not found in file.");
		}
		if (altColumn == -1) {
			throw new IOException("Column '" + COLUMN_ALT + "' not found in file.");
		}

	}

	public SitesEntry findByPosition(String chromosome, int position) throws IOException {
		TabixReader.Iterator iterator = tabixReader.query(chromosome, position - 1, position);
		String line = iterator.next();
		if (line == null) {
			return null;
		}

		//TODO: handle multiple matches with alleles... This implementation: always first, before always last.
		return parseLine(line);
	}

	public void close() {
		tabixReader.close();
	}

	protected SitesEntry parseLine(String line) {
		String[] tiles = line.split("\t");

		entry.setRsId(tiles[idColumn]);
		entry.setRefAllele(tiles[refColumn].charAt(0));
		entry.setAltAllele(tiles[altColumn].charAt(0));
		entry.setType("-");

		float aaf = 0;

		if (popColumn != -1) {
			if (!tiles[popColumn].equals(".")) {
				aaf = Float.parseFloat(tiles[popColumn]);
				entry.setFrequencies(true);
			} else {
				entry.setFrequencies(false);
			}
		} else {
			entry.setFrequencies(false);
		}

		entry.setRefFrequency(1 - aaf);
		entry.setAltFrequency(aaf);

		return entry;
	}

	private static DataInputStream openTxtOrGzipStream(String filename) throws IOException {
		FileInputStream inputStream = new FileInputStream(filename);
		InputStream in2 = FileUtil.decompressStream(inputStream);
		return new DataInputStream(in2);
	}

}

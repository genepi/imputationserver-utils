package genepi.imputationserver.steps.fastqc.legend;

import genepi.io.FileUtil;
import genepi.io.text.LineReader;
import htsjdk.tribble.readers.TabixReader;

import java.io.*;

public class LegendFileReader {

	private LegendEntry entry = new LegendEntry();

	private String population;

	private TabixReader tabixReader;

	private int idCol = -1;
	private int a0Col = -1;
	private int a1Col = -1;
	private int popCol = -1;

	public LegendFileReader(String filename, String population) throws IOException {
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
			if (tile.equals("id")) {
				idCol = i;
			}
			if (tile.equals("a0")) {
				a0Col = i;
			}
			if (tile.equals("a1")) {
				a1Col = i;
			}
			if (tile.equals(population + ".aaf")) {
				popCol = i;
			}
			i++;
		}

		// Validation: Check if all required columns have been found
		if (idCol == -1) {
			throw new IOException("Column 'id' not found in file.");
		}
		if (a0Col == -1) {
			throw new IOException("Column 'a0' not found in file.");
		}
		if (a1Col == -1) {
			throw new IOException("Column 'a1' not found in file.");
		}

	}

	public LegendEntry findByPosition(String chromosome, int position) throws IOException {
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

	protected LegendEntry parseLine(String line) {
		String[] tiles = line.split("\t");

		entry.setRsId(tiles[idCol]);
		entry.setAlleleA(tiles[a0Col].charAt(0));
		entry.setAlleleB(tiles[a1Col].charAt(0));
		entry.setType("-");

		float aaf = 0;

		if (popCol != -1) {
			if (!tiles[popCol].equals(".")) {
				aaf = Float.parseFloat(tiles[popCol]);
				entry.setFrequencies(true);
			} else {
				entry.setFrequencies(false);
			}
		} else {
			entry.setFrequencies(false);
		}

		entry.setFrequencyA(1 - aaf);
		entry.setFrequencyB(aaf);

		return entry;
	}

	private static DataInputStream openTxtOrGzipStream(String filename) throws IOException {
		FileInputStream inputStream = new FileInputStream(filename);
		InputStream in2 = FileUtil.decompressStream(inputStream);
		return new DataInputStream(in2);
	}

}

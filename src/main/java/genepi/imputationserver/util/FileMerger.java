package genepi.imputationserver.util;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import genepi.io.text.LineReader;

public class FileMerger {

	public static void mergeAndGzInfo(List<String> hdfs, String local) throws IOException {

		GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(local));

		boolean firstFile = true;

		for (String file : hdfs) {

			DataInputStream in = new DataInputStream(new FileInputStream(file));

			LineReader reader = new LineReader(in);

			boolean header = true;

			while (reader.next()) {

				String line = reader.get();

				if (header) {
					if (firstFile) {
						out.write(line.toString().getBytes());
						firstFile = false;
					}
					header = false;
				} else {
					out.write('\n');
					out.write(line.toString().getBytes());
				}
			}

			in.close();

		}

		out.close();
	}

	public static boolean keepVcfLineByInfo(String info, String field, double value) {
		String[] tilesInfo = info.split(";");
		for (String tile : tilesInfo) {
			String[] tilesRsq = tile.split("=");
			if (tilesRsq.length == 2) {
				String id = tilesRsq[0];
				String stringValue = tilesRsq[1];
				if (id.equals(field)) {
					double doubleValue = Double.parseDouble(stringValue);
					return (doubleValue > value);
				}
			}
		}
		return true;
	}

	public static String parseInfo(String line) {
		// rsq set. parse line and check rsq in info
		String[] tiles = line.split("\t", 9);
		if (tiles.length == 9) {
			String info = tiles[7];
			return info;
		} else {
			return null;
		}
	}

}

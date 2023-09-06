package genepi.imputationserver.util;

import java.io.IOException;
import java.util.ArrayList;

import genepi.io.text.LineReader;

public class CloudgeneContext {

	ArrayList<String> lines = new ArrayList<String>();

	public CloudgeneContext(String filename) {
		LineReader reader;
		try {
			reader = new LineReader(filename);
			while (reader.next()) {
				lines.add(reader.get());
			}
			reader.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public boolean hasInMemory(String content) {
		for (String str : lines) {
			if (str.contains(content)) {
				return true;
			}
		}

		return false;
	}

	public void print() {
		for (String str : lines) {
			System.out.println(str);
		}
	}

}

package genepi.imputationserver.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

public class CloudgeneContext {

	ArrayList<String> lines = new ArrayList<String>();

	public CloudgeneContext(String filename) {
		Scanner scanner;
		try {
			scanner = new Scanner(new File(filename)).useDelimiter(",\\s*");
			while (scanner.hasNext()) {
				lines.add(scanner.next());
			}
			scanner.close();
		} catch (FileNotFoundException e) {
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

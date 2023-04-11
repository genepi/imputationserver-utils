package genepi.imputationserver.steps.vcf;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import genepi.command.Command;
import genepi.io.FileUtil;
import genepi.io.text.LineReader;
import htsjdk.variant.vcf.VCFFileReader;

public class VcfFileUtil {

	public static String TABIX_PATH = "bin/";

	public static void setTabixBinary(String binaries) {
		TABIX_PATH = binaries;
	}
	
	public static String getBinary(){
		return TABIX_PATH;
	}
	
	public static VcfFile load(String vcfFilename, int chunksize, boolean createIndex) throws IOException {

		Set<Integer> chunks = new HashSet<Integer>();
		Set<String> chromosomes = new HashSet<String>();
		Set<String> rawChromosomes = new HashSet<String>();
		int noSnps = 0;
		int noSamples = 0;

		try {

			VCFFileReader reader = new VCFFileReader(new File(vcfFilename), false);

			noSamples = reader.getFileHeader().getGenotypeSamples().size();

			reader.close();

			LineReader lineReader = new LineReader(vcfFilename);

			boolean phased = true;
			boolean phasedAutodetect = true;
			boolean firstLine = true;
			while (lineReader.next()) {

				String line = lineReader.get();

				if (!line.startsWith("#")) {

					String tiles[] = line.split("\t", 10);

					if (tiles.length < 3) {
						throw new IOException("The provided VCF file is not tab-delimited");
					}

					String chromosome = tiles[0];
					rawChromosomes.add(chromosome);
					chromosome = chromosome.replaceAll("chr", "");
					int position = Integer.parseInt(tiles[1]);

					if (phased) {
						boolean containsSymbol = tiles[9].contains("/");

						if (containsSymbol) {
							phased = false;
						}

					}

					if (firstLine) {
						boolean containsSymbol = tiles[9].contains("/") || tiles[9].contains(".");

						if (!containsSymbol) {
							phasedAutodetect = true;
						} else {
							phasedAutodetect = false;
						}
						firstLine = false;

					}

					// TODO: check that all are phased
					// context.getGenotypes().get(0).isPhased();
					chromosomes.add(chromosome);
					if (chromosomes.size() > 1) {
						throw new IOException(
								"The provided VCF file contains more than one chromosome. Please split your input VCF file by chromosome");
					}

					String ref = tiles[3];
					String alt = tiles[4];

					if (ref.equals(alt)) {
						throw new IOException("The provided VCF file is malformed at variation " + tiles[2]
								+ ": reference allele (" + ref + ") and alternate allele  (" + alt + ") are the same.");
					}

					int chunk = position / chunksize;
					if (position % chunksize == 0) {
						chunk = chunk - 1;
					}
					chunks.add(chunk);
					noSnps++;

				} else {

					if (line.startsWith("#CHROM")) {

						String[] tiles = line.split("\t");

						// check sample names, stop when not unique
						HashSet<String> samples = new HashSet<>();

						for (int i = 0; i < tiles.length; i++) {

							String sample = tiles[i];

							if (samples.contains(sample)) {
								reader.close();
								throw new IOException("Two individuals or more have the following ID: " + sample);
							}
							samples.add(sample);
						}
					}

				}

			}
			lineReader.close();

			// create index
			if (createIndex && !new File(vcfFilename + ".tbi").exists()) {

				File command = new File(TABIX_PATH);
				if (!command.canExecute()){
					command.setExecutable(true, false);
				}
				
				Command tabix = new Command(TABIX_PATH);
				tabix.setParams("-f", "-p", "vcf", vcfFilename);
				tabix.saveStdErr("tabix.output");
				int returnCode = tabix.execute();

				if (returnCode != 0) {
					throw new IOException(
							"The provided VCF file is malformed. Error during index creation: "
									+ FileUtil.readFileAsString("tabix.output"));
				}

			}

			VcfFile pair = new VcfFile();
			pair.setVcfFilename(vcfFilename);
			pair.setIndexFilename(vcfFilename + ".tbi");
			pair.setNoSnps(noSnps);
			pair.setNoSamples(noSamples);
			pair.setChunks(chunks);
			pair.setChromosomes(chromosomes);
			
			boolean hasChrPrefix = false;
			for (String chromosome: rawChromosomes){
				if (chromosome.startsWith("chr")){
					hasChrPrefix = true;
				}
			}
			pair.setRawChromosomes(rawChromosomes);
			pair.setChrPrefix(hasChrPrefix);
			pair.setPhased(phased);
			pair.setPhasedAutodetect(phasedAutodetect);
			pair.setChunkSize(chunksize);
			return pair;

		} catch (Exception e) {
			throw new IOException(e.getMessage());
		}

	}

	public static Set<String> validChromosomes = new HashSet<String>();

	static {

		validChromosomes.add("1");
		validChromosomes.add("2");
		validChromosomes.add("3");
		validChromosomes.add("4");
		validChromosomes.add("5");
		validChromosomes.add("6");
		validChromosomes.add("7");
		validChromosomes.add("8");
		validChromosomes.add("9");
		validChromosomes.add("10");
		validChromosomes.add("11");
		validChromosomes.add("12");
		validChromosomes.add("13");
		validChromosomes.add("14");
		validChromosomes.add("15");
		validChromosomes.add("16");
		validChromosomes.add("17");
		validChromosomes.add("18");
		validChromosomes.add("19");
		validChromosomes.add("20");
		validChromosomes.add("21");
		validChromosomes.add("22");
		validChromosomes.add("23");
		validChromosomes.add("X");
		validChromosomes.add("MT");

		validChromosomes.add("chr1");
		validChromosomes.add("chr2");
		validChromosomes.add("chr3");
		validChromosomes.add("chr4");
		validChromosomes.add("chr5");
		validChromosomes.add("chr6");
		validChromosomes.add("chr7");
		validChromosomes.add("chr8");
		validChromosomes.add("chr9");
		validChromosomes.add("chr10");
		validChromosomes.add("chr11");
		validChromosomes.add("chr12");
		validChromosomes.add("chr13");
		validChromosomes.add("chr14");
		validChromosomes.add("chr15");
		validChromosomes.add("chr16");
		validChromosomes.add("chr17");
		validChromosomes.add("chr18");
		validChromosomes.add("chr19");
		validChromosomes.add("chr20");
		validChromosomes.add("chr21");
		validChromosomes.add("chr22");
		validChromosomes.add("chr23");
		validChromosomes.add("chrX");
		validChromosomes.add("chrMT");

	
	}

	public static boolean isValidChromosome(String chromosome) {
		return validChromosomes.contains(chromosome);
	}

	public static boolean isChrX(String chromosome) {
		return chromosome.equals("X") || chromosome.equals("23") || chromosome.equals("chrX") || chromosome.equals("chr23");
	}
	
	public static boolean isChrMT(String chromosome) {
		return chromosome.equals("MT") || chromosome.equals("chrMT");
	}

	public static void createIndex(String vcfFilename, boolean force) throws IOException{
		if (force){
			if (new File(vcfFilename + ".tbi").exists()){
				new File(vcfFilename + ".tbi").delete();
			}
		}
		if (!new File(vcfFilename + ".tbi").exists()) {
			createIndex(vcfFilename);
		}
	}


	public static void createIndex(String vcfFilename) throws IOException{
		Command tabix = new Command(TABIX_PATH);
		tabix.setParams("-f", "-p", "vcf", vcfFilename);
		tabix.saveStdErr("tabix.output");
		int returnCode = tabix.execute();

		if (returnCode != 0) {
			throw new IOException(
					"The provided VCF file is malformed. Error during index creation: "
							+ FileUtil.readFileAsString("tabix.output"));
		}
		
	}
}
package genepi.imputationserver.steps;

import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.Callable;

import genepi.imputationserver.steps.vcf.MergedVcfFile;
import genepi.imputationserver.util.FileChecksum;
import genepi.imputationserver.util.FileMerger;
import genepi.imputationserver.util.ImputationResults;
import genepi.imputationserver.util.ImputedChromosome;
import genepi.imputationserver.util.report.CloudgeneReport;
import genepi.io.FileUtil;
import genepi.io.text.LineWriter;
import net.lingala.zip4j.ZipFile;
import net.lingala.zip4j.exception.ZipException;
import net.lingala.zip4j.model.ZipParameters;
import net.lingala.zip4j.model.enums.AesKeyStrength;
import net.lingala.zip4j.model.enums.CompressionLevel;
import net.lingala.zip4j.model.enums.CompressionMethod;
import net.lingala.zip4j.model.enums.EncryptionMethod;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

@Command
public class EncryptionCommand implements Callable<Integer> {

	@Option(names = "--input", description = "Folder with imputed chunks", required = true)
	private String input;

	@Option(names = "--output", description = "Output Directory", required = true)
	private String output;

	@Option(names = "--password", description = "Password", required = true)
	private String password;

	@Option(names = "--meta", description = "Merge meta files", required = true)
	private String meta;

	@Option(names = "--mode", description = "Mode", required = true)
	private String mode;

	@Option(names = "--aesEncryption", description = "Use aesEncryption", required = true)
	private String aesEncryptionValue;

	@Option(names = "--reference", description = "Reference Panel ID", required = true)
	private String reference;

	@Option(names = "--phasing", description = "Phasing method", required = true)
	private String phasing;

	@Option(names = "--report", description = "Cloudgene Report Output", required = false)
	private String reportFileanme = "cloudgene.report.json";

	@Override
	public Integer call() throws Exception {

		CloudgeneReport report = new CloudgeneReport();
		report.setFilename(reportFileanme);

		boolean phasingOnly = false;
		if (mode != null && mode.equals("phasing")) {
			phasingOnly = true;
		}

		boolean mergeMetaFiles = !phasingOnly && (meta != null && meta.equals("yes"));

		boolean aesEncryption = (aesEncryptionValue != null && aesEncryptionValue.equals("yes"));

		try {

			// get sorted directories (
			List<String> folders = getDirectories(input);

			ImputationResults imputationResults = new ImputationResults(folders, phasingOnly);
			Map<String, ImputedChromosome> imputedChromosomes = imputationResults.getChromosomes();

			Set<String> chromosomes = imputedChromosomes.keySet();

			for (String name : chromosomes) {

				ImputedChromosome imputedChromosome = imputedChromosomes.get(name);

				report.println("Export and merge chromosome " + name);

				// create temp dir
				String temp = FileUtil.path(output, "temp");
				FileUtil.createDirectory(temp);

				// output files

				ArrayList<File> files = new ArrayList<File>();

				// merge info files
				if (!phasingOnly) {
					String infoOutput = FileUtil.path(temp, "chr" + name + ".info.gz");
					FileMerger.mergeAndGzInfo(imputedChromosome.getInfoFiles(), infoOutput);
					files.add(new File(infoOutput));
				}

				// merge all dosage files

				String dosageOutput;
				if (phasingOnly) {
					dosageOutput = FileUtil.path(temp, "chr" + name + ".phased.vcf.gz");
				} else {
					dosageOutput = FileUtil.path(temp, "chr" + name + ".dose.vcf.gz");
				}
				files.add(new File(dosageOutput));

				MergedVcfFile vcfFile = new MergedVcfFile(dosageOutput);
				vcfFile.addHeader(report, imputedChromosome.getHeaderFiles());

				for (String file : imputedChromosome.getDataFiles()) {
					report.println("Read file " + file);
					vcfFile.addFile(new FileInputStream(file));
					FileUtil.deleteFile(file);
				}

				vcfFile.close();

				// merge all meta files
				if (mergeMetaFiles) {

					report.println("Merging meta files...");

					String dosageMetaOutput = FileUtil.path(temp, "chr" + name + ".empiricalDose.vcf.gz");
					MergedVcfFile vcfFileMeta = new MergedVcfFile(dosageMetaOutput);

					String headerMetaFile = imputedChromosome.getHeaderMetaFiles().get(0);
					report.println("Use header from file " + headerMetaFile);

					vcfFileMeta.addFile(new FileInputStream(headerMetaFile));

					for (String file : imputedChromosome.getDataMetaFiles()) {
						report.println("Read file " + file);
						vcfFileMeta.addFile(new FileInputStream(file));
						FileUtil.deleteFile(file);
					}
					vcfFileMeta.close();

					report.println("Meta files merged.");

					files.add(new File(dosageMetaOutput));
				}

				// create zip file
				String fileName = "chr_" + name + ".zip";
				String filePath = FileUtil.path(output, fileName);
				File file = new File(filePath);
				createEncryptedZipFile(file, files, password, aesEncryption);

				// add checksum to hash file
				String checksumFilename = filePath + ".md5";
				LineWriter writer = new LineWriter(checksumFilename);
				report.log("Creating file checksum for " + filePath);
				long checksumStart = System.currentTimeMillis();
				String checksum = FileChecksum.HashFile(new File(filePath), FileChecksum.Algorithm.MD5);
				writer.write(checksum + " " + fileName);
				long checksumEnd = (System.currentTimeMillis() - checksumStart) / 1000;
				report.log("File checksum for " + filePath + " created in " + checksumEnd + " seconds.");
				writer.close();

				// delete temp dir
				FileUtil.deleteDirectory(temp);

				report.ok("Chromosome " + name + " encrypted.");

			}

		} catch (Exception e) {
			e.printStackTrace();
			report.error("Data export failed: " + e.getMessage());
			return -1;
		}

		// submit counters!
		report.submitCounter("samples");
		report.submitCounter("genotypes");
		report.submitCounter("chromosomes");
		report.submitCounter("runs");

		// submit panel and phasing method counters
		report.submitCounter("refpanel_" + reference);
		report.submitCounter("phasing_" + phasing);
		report.submitCounter("23andme-input");

		return 0;

	}

	private void createEncryptedZipFile(File file, List<File> files, String password, boolean aesEncryption)
			throws ZipException {
		ZipParameters param = new ZipParameters();
		param.setEncryptFiles(true);
		param.setEncryptionMethod(EncryptionMethod.ZIP_STANDARD);

		if (aesEncryption) {
			param.setEncryptionMethod(EncryptionMethod.AES);
			param.setAesKeyStrength(AesKeyStrength.KEY_STRENGTH_256);
			param.setCompressionMethod(CompressionMethod.DEFLATE);
			param.setCompressionLevel(CompressionLevel.NORMAL);
		}

		ZipFile zipFile = new ZipFile(file, password.toCharArray());
		zipFile.addFiles(files, param);

	}

	private List<String> getDirectories(String directory) {

		List<String> result = new Vector<>();
		File[] files = new File(directory).listFiles();
		for (File file : files) {
			if (file.isDirectory()) {
				result.add(file.getAbsolutePath());
			}
		}
		Collections.sort(result);
		return result;
	}

}

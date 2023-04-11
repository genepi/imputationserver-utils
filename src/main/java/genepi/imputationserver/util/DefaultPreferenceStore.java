package genepi.imputationserver.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.Properties;
import java.util.Set;

public class DefaultPreferenceStore {

	private Properties properties = new Properties(defaults());

	public void load(File file) {
		try {
			properties.load(new FileInputStream(file));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public DefaultPreferenceStore() {

	}

	public String getString(String key) {
		return properties.getProperty(key);
	}

	public void setString(String key, String value) {
		properties.setProperty(key, value);
	}

	public Set<Object> getKeys() {
		return new HashSet<Object>(Collections.list(properties.propertyNames()));
	}

	public static Properties defaults() {

		Properties defaults = new Properties();
		defaults.setProperty("chunksize", "20000000");
		defaults.setProperty("phasing.window", "5000000");
		defaults.setProperty("minimac.window", "500000");
		defaults.setProperty("minimac.sendmail", "no");
		defaults.setProperty("server.url", "https://imputationserver.sph.umich.edu");
		defaults.setProperty("minimac.tmp", "/tmp");
		defaults.setProperty("minimac.command",
				"--refHaps ${ref} --haps ${vcf} --start ${start} --end ${end} --window ${window} --prefix ${prefix} --chr ${chr} --cpus 1 --noPhoneHome --format GT,DS,GP --allTypedSites --meta --minRatio 0.00001 ${chr =='MT' ? '--myChromosome ' + chr : ''} ${unphased ? '--unphasedOutput' : ''} ${mapMinimac != null ? '--referenceEstimates --map ' + mapMinimac : ''}");
		defaults.setProperty("eagle.command",
				"--vcfRef ${ref} --vcfTarget ${vcf} --geneticMapFile ${map} --outPrefix ${prefix} --bpStart ${start} --bpEnd ${end} --allowRefAltSwap --vcfOutFormat z --keepMissingPloidyX");
		defaults.setProperty("beagle.command",
				"-jar ${beagle} ref=${ref} gt=${vcf} out=${prefix} nthreads=1 chrom=${chr}:${start}-${end} map=${map} impute=false");
		defaults.setProperty("ref.fasta", "v37");
		defaults.setProperty("contact.name", "Christian Fuchsberger");
		defaults.setProperty("contact.email", "cfuchsb@umich.edu");
		defaults.setProperty("hg38Tohg19", "chains/hg38ToHg19.over.chain.gz");
		defaults.setProperty("hg19Tohg38", "chains/hg19ToHg38.over.chain.gz");
		defaults.setProperty("sanitycheck", "yes");

		return defaults;
	}

}

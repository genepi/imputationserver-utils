package genepi.imputationserver.util;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.esotericsoftware.yamlbeans.YamlReader;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonIOException;
import com.google.gson.JsonSyntaxException;

public class RefPanel {

	public static final String ALLELE_SWITCHES = String.valueOf(Integer.MAX_VALUE);
	public static final String STRAND_FLIPS = "100";
	public static final String SAMPLE_CALL_RATE = "0.5";
	public static final String MIN_SNPS = "3";
	public static final String OVERLAP = "0.5";
	public static final String CHR_X_MIXED_GENOTYPES = "0.1";

	private String id;

	private String genotypes;

	private String legend;

	private String mapMinimac;

	private String build = "hg19";

	private String mapEagle;

	private String refEagle;

	private String refBeagle;

	private String mapBeagle;

	private List<RefPanelPopulation> populations;

	private Map<String, String> defaultQcFilter;

	private Map<String, String> qcFilter;

	private String range;

	/**
	 * 
	 */
	public RefPanel() {
		defaultQcFilter = new HashMap<String, String>();
		defaultQcFilter.put("overlap", OVERLAP);
		defaultQcFilter.put("minSnps", MIN_SNPS);
		defaultQcFilter.put("sampleCallrate", SAMPLE_CALL_RATE);
		defaultQcFilter.put("mixedGenotypeschrX", CHR_X_MIXED_GENOTYPES);
		defaultQcFilter.put("strandFlips", STRAND_FLIPS);
		defaultQcFilter.put("alleleSwitches", ALLELE_SWITCHES);
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getGenotypes() {
		return genotypes;
	}

	public void setGenotypes(String genotypes) {
		this.genotypes = genotypes;
	}

	public String getLegend() {
		return legend;
	}

	public void setLegend(String legend) {
		this.legend = legend;
	}

	public void setBuild(String build) {
		this.build = build;
	}

	public String getBuild() {
		return build;
	}

	public void setMapMinimac(String mapMinimac) {
		this.mapMinimac = mapMinimac;
	}

	public String getMapMinimac() {
		return mapMinimac;
	}

	public void setMapEagle(String mapEagle) {
		this.mapEagle = mapEagle;
	}

	public String getMapEagle() {
		return mapEagle;
	}

	public String getRefBeagle() {
		return refBeagle;
	}

	public void setRefBeagle(String refBeagle) {
		this.refBeagle = refBeagle;
	}

	public void setRefEagle(String refEagle) {
		this.refEagle = refEagle;
	}

	public String getRefEagle() {
		return refEagle;
	}

	public int getSamplesByPopulation(String population) {
		if (population == null) {
			return 0;
		}

		RefPanelPopulation panelPopulation = getPopulation(population);
		if (panelPopulation != null) {
			return panelPopulation.getSamples();
		} else {
			return 0;
		}
	}

	public void setPopulations(List<RefPanelPopulation> populations) {
		this.populations = populations;
	}

	public List<RefPanelPopulation> getPopulations() {
		return populations;
	}

	public boolean supportsPopulation(String population) {
		if (population == null || population.equals("")) {
			return true;
		}
		return (getPopulation(population) != null);
	}

	public RefPanelPopulation getPopulation(String id) {
		if (populations == null) {
			return null;
		}

		for (RefPanelPopulation population : populations) {
			if (population.getId().equalsIgnoreCase(id)) {
				return population;
			}
		}
		return null;
	}

	public Map<String, String> getQcFilter() {
		return qcFilter;
	}

	public double getQcFilterByKey(String key) {

		if (qcFilter == null) {
			qcFilter = defaultQcFilter;
		}
		Object n = qcFilter.get(key);

		if (n == null) {
			return Double.parseDouble(defaultQcFilter.get(key));
		}

		if (n instanceof String) {
			return Double.parseDouble((String) n);
		} else if (n instanceof Number) {
			return ((Number) n).doubleValue();
		} else {
			throw new IllegalArgumentException("Value associated with key is neither a String nor a Number: " + n);
		}

	}

	public void setQcFilter(Map<String, String> qcFilter) {
		this.qcFilter = qcFilter;
	}

	public void setRange(String range) {
		this.range = range;
	}

	public String getRange() {
		return range;
	}

	public String getMapBeagle() {
		return mapBeagle;
	}

	public void setMapBeagle(String mapBeagle) {
		this.mapBeagle = mapBeagle;
	}

	public static RefPanel fromProperties(Object properties) throws IOException {

		if (properties == null) {
			throw new IOException("Propertie map not set.");
		}

		RefPanel panel = new RefPanel();
		Map<String, Object> map = (Map<String, Object>) properties;

		if (map.get("genotypes") != null) {
			panel.setGenotypes(map.get("genotypes").toString());
		} else {
			throw new IOException("Property 'genotypes' not found in cloudgene.yaml.");
		}

		if (map.get("id") != null) {
			panel.setId(map.get("id").toString());
		} else {
			throw new IOException("Property 'id' not found in cloudgene.yaml.");
		}

		if (map.get("legend") != null) {
			panel.setLegend(map.get("legend").toString());
		} else {
			throw new IOException("Property 'legend' not found in cloudgene.yaml.");
		}

		if (map.get("mapEagle") != null) {
			panel.setMapEagle(map.get("mapEagle").toString());
		}

		if (map.get("refEagle") != null) {
			panel.setRefEagle(map.get("refEagle").toString());
		}

		if (map.get("mapBeagle") != null) {
			panel.setMapBeagle(map.get("mapBeagle").toString());
		}

		if (map.get("refBeagle") != null) {
			panel.setRefBeagle(map.get("refBeagle").toString());
		}

		if (map.get("populations") != null) {
			panel.setPopulations(RefPanelPopulation.fromProperties((List<Map<String, Object>>) map.get("populations")));
		} else {
			throw new IOException("Property 'populations' not found in cloudgene.yaml.");
		}

		if (map.get("qcFilter") != null) {
			panel.setQcFilter((Map<String, String>) map.get("qcFilter"));
		}

		// optional parameters
		if (map.get("reference_build") != null) {
			panel.setBuild(map.get("reference_build").toString());
		}

		if (map.get("build") != null) {
			panel.setBuild(map.get("build").toString());
		}

		if (map.get("range") != null) {
			panel.setRange(map.get("range").toString());
		} else {
			panel.setRange(null);
		}

		if (map.get("mapMinimac") != null) {
			panel.setMapMinimac(map.get("mapMinimac").toString());
		} else {
			panel.setMapMinimac(null);
		}

		return panel;

	}

	public static void resolveEnvVariable(Map<String, Object> properties, String folder) {
		for (String key : properties.keySet()) {
			Object value = properties.get(key);
			if (value instanceof String) {
				String valueString = value.toString().replaceAll("\\$\\{app_local_folder\\}", folder);
				properties.put(key, valueString);
			}
		}
	}

	public static RefPanel loadFromJson(String filename) throws JsonSyntaxException, JsonIOException, IOException {
		Gson gson = (new GsonBuilder()).create();
		Map<String, Object> panel = gson.fromJson(new FileReader(filename), Map.class);
		return fromProperties(panel);
	}

	public static RefPanel loadFromYamlFile(String filename, String id) throws IOException {

		YamlReader reader = new YamlReader(new FileReader(filename));
		List<Map<String, Object>> panels = reader.read(List.class);
		for (Map<String, Object> panel : panels) {
			if (panel.get("id").equals(id)) {
				resolveEnvVariable(panel, (new File(filename).getParent()));
				return fromProperties(panel);
			}
		}

		throw new IOException("Reference panel '" + id + "' not found in file '" + filename + "'.");

	}

}

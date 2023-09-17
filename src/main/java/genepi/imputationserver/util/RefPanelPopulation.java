package genepi.imputationserver.util;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Vector;

public class RefPanelPopulation {

	private String id;

	private String name;

	private int samples;

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int getSamples() {
		return samples;
	}

	public void setSamples(int samples) {
		this.samples = samples;
	}

	public static List<RefPanelPopulation> fromProperties(List<Map<String, Object>> properties) throws IOException {
		List<RefPanelPopulation> result = new Vector<RefPanelPopulation>();
		for (Map<String, Object> property : properties) {
			RefPanelPopulation population = new RefPanelPopulation();
			if (property.containsKey("id")) {
				population.setId(property.get("id").toString());
			} else {
				throw new IOException("Property 'id' not found in population list.");
			}
			if (property.containsKey("name")) {
				population.setName(property.get("name").toString());
			} else {
				throw new IOException("Property 'name' not found in population list.");
			}
			if (property.containsKey("samples")) {
				Double samples = Double.parseDouble(property.get("samples").toString());
				population.setSamples(samples.intValue());
			} else {
				throw new IOException("Property 'samples' not found in population list.");
			}
			result.add(population);
		}
		return result;

	}

}

package genepi.imputationserver.steps.fastqc.legend;

public class SitesEntry {

	private char refAllele;

	private char altAllele;

	private String rsId;

	private float refFrequency;

	private float altFrequency;

	private boolean frequencies = false;;

	private String type;
	
	private String genotype = null;

	public char getRefAllele() {
		return refAllele;
	}

	public void setRefAllele(char refAllele) {
		this.refAllele = refAllele;
		this.genotype = null;
	}

	public char getAltAllele() {
		return altAllele;
	}

	public void setAltAllele(char altAllele) {
		this.altAllele = altAllele;
		this.genotype = null;
	}

	public float getRefFrequency() {
		return refFrequency;
	}

	public void setRefFrequency(float refFrequency) {
		this.refFrequency = refFrequency;
	}

	public float getAltFrequency() {
		return altFrequency;
	}

	public void setAltFrequency(float altFrequency) {
		this.altFrequency = altFrequency;
	}

	public String getRsId() {
		return rsId;
	}

	public void setRsId(String rsId) {
		this.rsId = rsId;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public void setFrequencies(boolean frequencies) {
		this.frequencies = frequencies;
	}

	public boolean hasFrequencies() {
		return frequencies;
	}
	
	public String getGenotype() {
		if (genotype == null) {
			StringBuilder builder = new StringBuilder(2);
			builder.append(refAllele);
			builder.append(altAllele);
			genotype = builder.toString();
		}
		return genotype;
	}

}

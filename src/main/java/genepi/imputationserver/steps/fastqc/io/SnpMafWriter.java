package genepi.imputationserver.steps.fastqc.io;

import genepi.imputationserver.steps.fastqc.SnpStats;
import genepi.imputationserver.steps.vcf.MinimalVariantContext;
import genepi.io.table.writer.CsvTableWriter;
import genepi.io.text.LineWriter;

import java.io.IOException;

public class SnpMafWriter {

    private CsvTableWriter writer;

    private String filename;

    public SnpMafWriter(String filename) throws IOException {
       this.filename = filename;
    }

    public void write(MinimalVariantContext snp, SnpStats stats) throws IOException {
        if (writer == null) {
            writer = new CsvTableWriter(filename, '\t', false);
            writer.setColumns(new String[]{"ID","CHROM","POS", "REF", "ALT",
                    "AAF", "REFERENCE_REF",
                    "REFERENCE_ALT", "REFERENCE_AAF",
                    "CHISQ", "OVERLAP_WITH_REFERENCE", "TYPE"});
        }

        writer.setString("ID", snp.getId());
        writer.setString("CHROM", stats.getChromosome());
        writer.setInteger("POS", stats.getPosition());
        writer.setString("REF", String.valueOf(stats.getAlleleA()));
        writer.setString("ALT", String.valueOf(stats.getAlleleB()));
        writer.setDouble("AAF", stats.getFrequencyB());
        writer.setString("REFERENCE_REF", stats.getRefAlleleA() != Byte.MAX_VALUE ? String.valueOf(stats.getRefAlleleA()) : "NA");
        writer.setString("REFERENCE_ALT", stats.getRefAlleleB() != Byte.MAX_VALUE ? String.valueOf(stats.getRefAlleleB()) : "NA");
        writer.setString("REFERENCE_AAF", Float.isNaN(stats.getRefFrequencyB()) ? "NA" : Float.toString(stats.getRefFrequencyB()));
        writer.setString("CHISQ", Double.isNaN(stats.getChisq()) ? "NA" : Double.toString(stats.getChisq()));
        writer.setString("OVERLAP_WITH_REFERENCE", String.valueOf(stats.isOverlapWithReference()));
        writer.setString("TYPE", stats.getType());
        writer.next();

    }

    public void close() throws IOException {
        if (writer != null) {
            writer.close();
        }
    }
}

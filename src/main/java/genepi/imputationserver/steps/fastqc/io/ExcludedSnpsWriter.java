package genepi.imputationserver.steps.fastqc.io;

import genepi.imputationserver.steps.vcf.MinimalVariantContext;
import genepi.io.table.writer.CsvTableWriter;
import genepi.io.text.LineWriter;

import java.io.IOException;

public class ExcludedSnpsWriter {

    private CsvTableWriter writer;

    private String filename;

    public ExcludedSnpsWriter(String filename) throws IOException {
        this.filename = filename;
    }

    public void write(MinimalVariantContext variant) throws IOException {
        write(variant, "");
    }

    public void write(MinimalVariantContext variant, String details) throws IOException {
        if (writer == null) {
            writer = new CsvTableWriter(filename, '\t', false);
            writer.setColumns(new String[]{"ID", "CHROM", "POS", "REF", "ALT", "INFO"});
        }
        writer.setString("ID", variant.getId());
        writer.setString("CHROM", variant.getContig());
        writer.setInteger("POS", variant.getStart());
        writer.setString("REF", variant.getReferenceAllele());
        writer.setString("ALT", variant.getAlternateAllele());
        writer.setString("INFO", details);
        writer.next();
    }

    public void close() throws IOException {
        if (writer != null) {
            writer.close();
        }
    }

}

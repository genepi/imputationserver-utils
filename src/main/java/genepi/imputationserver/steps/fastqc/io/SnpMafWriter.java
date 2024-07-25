package genepi.imputationserver.steps.fastqc.io;

import genepi.imputationserver.steps.fastqc.SnpStats;
import genepi.imputationserver.steps.vcf.MinimalVariantContext;
import genepi.io.text.LineWriter;

import java.io.IOException;

public class SnpMafWriter {

    private LineWriter writer;

    private String filename;

    public SnpMafWriter(String filename) throws IOException {
       this.filename = filename;
    }

    public void write(MinimalVariantContext snp, SnpStats statistics) throws IOException {
        if (writer == null) {
            writer = new LineWriter(filename);

        }
        writer.write(snp.getId() + "\t" + statistics.toString());
    }

    public void close() throws IOException {
        if (writer != null) {
            writer.close();
        }
    }
}

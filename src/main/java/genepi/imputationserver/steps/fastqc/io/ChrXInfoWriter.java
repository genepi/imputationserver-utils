package genepi.imputationserver.steps.fastqc.io;

import genepi.imputationserver.steps.vcf.MinimalVariantContext;
import genepi.io.text.LineWriter;

import java.io.IOException;

public class ChrXInfoWriter {

    private LineWriter writer;

    private String filename;

    public ChrXInfoWriter(String filename) throws IOException {
        this.filename = filename;
    }

    public void write(String name, String details) throws IOException {
        if (writer == null) {
            writer = new LineWriter(filename);
            writer.write("#Sample\tPosition", false);
        }
        writer.write(name + "\t" + details);
    }

    public void close() throws IOException {
        if (writer != null) {
            writer.close();
        }
    }

}

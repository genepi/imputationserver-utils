package genepi.imputationserver.steps.fastqc.io;

import genepi.imputationserver.steps.vcf.MinimalVariantContext;
import genepi.io.table.writer.CsvTableWriter;
import genepi.io.text.LineWriter;

import java.io.IOException;

public class ChrXInfoWriter {

    private CsvTableWriter writer;

    private String filename;

    public ChrXInfoWriter(String filename) throws IOException {
        this.filename = filename;
    }

    public void write(String sample, String position) throws IOException {
        if (writer == null) {
            writer = new CsvTableWriter(filename, '\t', false);
            writer.setColumns(new String[]{"SAMPLE", "POSITION"});
        }
        writer.setString("SAMPLE", sample);
        writer.setString("POSITION", position);
        writer.next();
    }

    public void close() throws IOException {
        if (writer != null) {
            writer.close();
        }
    }

}

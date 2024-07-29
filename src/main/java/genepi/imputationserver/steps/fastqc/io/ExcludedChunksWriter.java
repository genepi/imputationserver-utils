package genepi.imputationserver.steps.fastqc.io;

import genepi.imputationserver.steps.vcf.MinimalVariantContext;
import genepi.imputationserver.steps.vcf.VcfChunk;
import genepi.io.table.writer.CsvTableWriter;
import genepi.io.text.LineWriter;

import java.io.IOException;

public class ExcludedChunksWriter {

    private CsvTableWriter writer;

    private String filename;

    public ExcludedChunksWriter(String filename) throws IOException {
        this.filename = filename;
    }

    public void write(VcfChunk chunk, double overlap, int countLowSamples) throws IOException {
        if (writer == null) {
            writer = new CsvTableWriter(filename, '\t', false);
            writer.setColumns(new String[]{"CHUNK", "SNPS", "REFERENCE_OVERLAP", "SAMPLES_LOW_CALL_RATE"});
        }
        writer.setString("CHUNK", chunk.toString());
        writer.setInteger("SNPS", chunk.overallSnpsChunk);
        writer.setDouble("REFERENCE_OVERLAP", overlap);
        writer.setInteger("SAMPLES_LOW_CALL_RATE", countLowSamples);
        writer.next();
    }

    public void close() throws IOException {
        if (writer != null) {
            writer.close();
        }
    }

}

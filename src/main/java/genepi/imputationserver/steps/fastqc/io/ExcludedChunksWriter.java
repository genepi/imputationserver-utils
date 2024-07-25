package genepi.imputationserver.steps.fastqc.io;

import genepi.imputationserver.steps.vcf.MinimalVariantContext;
import genepi.imputationserver.steps.vcf.VcfChunk;
import genepi.io.text.LineWriter;

import java.io.IOException;

public class ExcludedChunksWriter {

    private LineWriter writer;

    private String filename;

    public ExcludedChunksWriter(String filename) throws IOException {
        this.filename = filename;
    }

    public void write(VcfChunk chunk, double overlap, int countLowSamples) throws IOException {
        if (writer == null) {
            writer = new LineWriter(filename);
            writer.write(
                    "#Chunk" + "\t" + "SNPs (#)" + "\t" + "Reference Overlap (%)" + "\t" + "Low Sample Call Rates (#)",
                    false);
        }
        writer.write(chunk.toString() + "\t" + chunk.overallSnpsChunk + "\t" + overlap + "\t" + countLowSamples);
    }

    public void close() throws IOException {
        if (writer != null) {
            writer.close();
        }
    }

}

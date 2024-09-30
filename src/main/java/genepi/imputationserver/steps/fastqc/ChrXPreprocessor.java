package genepi.imputationserver.steps.fastqc;

import genepi.imputationserver.steps.fastqc.io.ChrXInfoWriter;
import genepi.imputationserver.steps.fastqc.io.ExcludedSnpsWriter;
import genepi.imputationserver.steps.vcf.MinimalVariantContext;
import genepi.imputationserver.util.GenomicTools;
import genepi.io.FileUtil;
import genepi.io.text.LineReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Vector;

public class ChrXPreprocessor {

    public static final String X_PAR1 = "X.PAR1";

    public static final String X_PAR2 = "X.PAR2";

    public static final String X_NON_PAR = "X.nonPAR";

    private final String filename;

    private final String build;

    private boolean chrXPloidyError = false;

    private boolean chrXMissingRate = false;

    private int invalidAlleles = 0;

    private int filtered = 0;

    public ChrXPreprocessor(String filename, String build) {
        this.filename = filename;
        this.build = build;
    }

    public List<String> splitVcf(double mixedGenotypesChrX, String output, String chrXInfoFile, ExcludedSnpsWriter excludedSnpsWriter) throws IOException {

        ChrXInfoWriter chrXInfoWriter = new ChrXInfoWriter(chrXInfoFile);

        HashSet<String> hapSamples = new HashSet<>();

        List<String> paths = new Vector<String>();
        String nonPar = FileUtil.path(output, X_NON_PAR + ".vcf.gz");
        VariantContextWriter vcfChunkWriterNonPar = new VariantContextWriterBuilder().setOutputFile(nonPar)
                .setOption(Options.INDEX_ON_THE_FLY).setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF).build();

        String par1 = FileUtil.path(output, X_PAR1 + ".vcf.gz");
        VariantContextWriter vcfChunkWriterPar1 = new VariantContextWriterBuilder().setOutputFile(par1)
                .setOption(Options.INDEX_ON_THE_FLY).setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF).build();

        String par2 = FileUtil.path(output, X_PAR2 + ".vcf.gz");
        VariantContextWriter vcfChunkWriterPar2 = new VariantContextWriterBuilder().setOutputFile(par2)
                .setOption(Options.INDEX_ON_THE_FLY).setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF).build();

        VCFFileReader vcfReader = new VCFFileReader(new File(filename), true);

        VCFHeader header = vcfReader.getFileHeader();
        vcfChunkWriterNonPar.writeHeader(header);
        vcfChunkWriterPar1.writeHeader(header);
        vcfChunkWriterPar2.writeHeader(header);

        int mixedGenotypes[] = null;
        int count = 0;

        int nonParStart = 2699520;
        int nonParEnd = 154931044;

        if (build.equals("hg38")) {
            nonParStart = 2781479;
            nonParEnd = 155701383;
        }

        VCFCodec codec = new VCFCodec();
        codec.setVCFHeader(vcfReader.getFileHeader(), VCFHeaderVersion.VCF4_1);
        LineReader reader = new LineReader(filename);

        while (reader.next()) {

            String lineString = reader.get();

            if (!lineString.startsWith("#")) {

                String tiles[] = lineString.split("\t", 6);
                String ref = tiles[3];
                String alt = tiles[4];

                // filter invalid alleles
                //TODO: could be deleted and handled central in StatisticTask. ATM needed for codec.decode
                if (!GenomicTools.isValid(ref) || !GenomicTools.isValid(alt)) {
                    MinimalVariantContext variantContext = new MinimalVariantContext(1);
                    variantContext.setId(tiles[2]);
                    variantContext.setContig(tiles[0]);
                    variantContext.setStart(Integer.parseInt(tiles[1]));
                    variantContext.setReferenceAllele(ref);
                    variantContext.setAlternateAllele(alt);
                    excludedSnpsWriter.write(variantContext, "Invalid Alleles");
                    invalidAlleles++;
                    filtered++;
                    continue;
                }

                // now decode, since it's a valid VCF
                VariantContext line = codec.decode(lineString);

                if (line.getContig().equals("23")) {
                    line = new VariantContextBuilder(line).chr("X").make();
                } else if (line.getContig().equals("chr23")) {
                    line = new VariantContextBuilder(line).chr("chrX").make();
                }

                if (line.getStart() < nonParStart) {

                    vcfChunkWriterPar1.add(line);

                    if (!paths.contains(par1)) {
                        paths.add(par1);
                    }

                }

                else if (line.getStart() >= nonParStart && line.getStart() <= nonParEnd) {

                    count++;

                    checkPloidy(header.getGenotypeSamples(), line, hapSamples, chrXInfoWriter);

                    mixedGenotypes = checkMixedGenotypes(mixedGenotypes, line);

                    vcfChunkWriterNonPar.add(line);

                    if (!paths.contains(nonPar)) {
                        paths.add(nonPar);
                    }

                }

                else {

                    vcfChunkWriterPar2.add(line);

                    if (!paths.contains(par2)) {
                        paths.add(par2);
                    }

                }

            }

        }

        if (mixedGenotypes != null) {
            for (int i = 0; i < mixedGenotypes.length; i++) {
                double missingRate = mixedGenotypes[i] / (double) count;
                if (missingRate > mixedGenotypesChrX) {
                    //TODO: throw exceptions
                    this.chrXMissingRate = true;
                    break;
                }

            }
        }

        vcfReader.close();
        reader.close();

        vcfChunkWriterPar1.close();
        vcfChunkWriterPar2.close();
        vcfChunkWriterNonPar.close();

        chrXInfoWriter.close();

        return paths;
    }

    // mixed genotype: ./1; 1/.;
    private int[] checkMixedGenotypes(int[] mixedGenotypes, VariantContext line) {

        //TODO: init outside at beginning?
        if (mixedGenotypes == null) {

            mixedGenotypes = new int[line.getNSamples()];
            for (int i = 0; i < line.getNSamples(); i++) {
                mixedGenotypes[i] = 0;
            }
        }

        for (int i = 0; i < line.getNSamples(); i++) {
            Genotype genotype = line.getGenotype(i);

            if (genotype.isMixed()) {
                mixedGenotypes[i] += 1;
            }

        }
        return mixedGenotypes;
    }

    private void checkPloidy(List<String> samples, VariantContext snp, HashSet<String> hapSamples, ChrXInfoWriter chrXInfoWriter) throws IOException {

        for (final String name : samples) {

            Genotype genotype = snp.getGenotype(name);

            if (hapSamples.contains(name) && genotype.getPloidy() != 1) {
                chrXInfoWriter.write(name, snp.getContig() + ":" + snp.getStart());
                //TODO: throw exceptions
                this.chrXPloidyError = true;
            }

            if (genotype.getPloidy() == 1) {
                hapSamples.add(name);
            }

        }
    }

    public boolean isChrXMissingRate() {
        return chrXMissingRate;
    }

    public boolean isChrXPloidyError() {
        return chrXPloidyError;
    }

    public int getInvalidAlleles() {
        return invalidAlleles;
    }

    public int getFiltered() {
        return filtered;
    }

    public static String getContig(String filename) {
        if (filename.contains(X_PAR1)) {
            return X_PAR1;
        } else if (filename.contains(X_PAR2)) {
            return X_PAR2;
        } else {
            return X_NON_PAR;
        }
    }

}

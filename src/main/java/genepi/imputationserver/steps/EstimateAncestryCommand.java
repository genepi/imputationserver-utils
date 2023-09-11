package genepi.imputationserver.steps;

import java.util.concurrent.Callable;

import genepi.imputationserver.steps.ancestry.PopulationPredictor;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

@Command
public class EstimateAncestryCommand implements Callable<Integer> {

	@Option(names = "--samples", description = "Reference Samples File", required = true)
	private String samples;

	@Option(names = "--reference-pc", description = "Reference PCs File", required = true)
	private String reference_pc;

	@Option(names = "--study-pc", description = "Study PCs File", required = true)
	private String study_pc;

	@Option(names = "--max-pcs", description = "Max PCs", required = true)
	private int max_pcs;

	@Option(names = "--k", description = "k Nearest Neighbors", required = true)
	private int k;

	@Option(names = "--threshold", description = "Threshold", required = true)
	private Double threshold;

	@Option(names = "--output", description = "Output File", required = true)
	private String output;

	@Override
	public Integer call() throws Exception {

		PopulationPredictor predictor = new PopulationPredictor();

		predictor.setSamplesFile(samples);
		predictor.setReferenceFile(reference_pc);
		predictor.setStudyFile(study_pc);
		predictor.setMaxPcs(max_pcs);
		predictor.setK(k);
		predictor.setWeightThreshold(threshold);

		predictor.predictPopulation(output);
		return 0;

	}

}

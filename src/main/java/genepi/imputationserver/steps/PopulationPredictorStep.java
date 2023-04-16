package genepi.imputationserver.steps;

import cloudgene.sdk.internal.WorkflowContext;
import cloudgene.sdk.internal.WorkflowStep;
import genepi.imputationserver.steps.ancestry.PopulationPredictor;

public class PopulationPredictorStep extends WorkflowStep {

	@Override
	public boolean run(WorkflowContext context) {

		PopulationPredictor predictor = new PopulationPredictor();
		predictor.setSamplesFile(context.get("samples"));
		predictor.setReferenceFile(context.get("reference_pc"));
		predictor.setStudyFile(context.get("study_pc"));
		int pcs = Integer.parseInt(context.get("max_pcs"));
		predictor.setMaxPcs(pcs);
		int k = Integer.parseInt(context.get("k"));
		predictor.setK(k);
		double threshold = Double.parseDouble(context.get("threshold"));
		predictor.setWeightThreshold(threshold);
		predictor.predictPopulation(context.get("output"));
		return true;

	}

}

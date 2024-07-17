package genepi.imputationserver.steps.fastqc;

public interface ITask {

	public String getName();
	
	public TaskResults run() throws Exception;

}

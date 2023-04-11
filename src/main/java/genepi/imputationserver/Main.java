package genepi.imputationserver;

import java.lang.reflect.InvocationTargetException;

import genepi.base.Toolbox;
import genepi.imputationserver.steps.imputation.ImputationPipeline;
import genepi.imputationserver.tools.VersionTool;

public class Main extends Toolbox {

	public Main(String command, String[] args) {
		super(command, args);
		
		System.out.println();
		System.out.println(ImputationPipeline.PIPELINE_VERSION);
		System.out.println("https://imputationserver.sph.umich.edu");
		System.out.println("Lukas Forer, Sebastian Schoenherr and Christian Fuchsberger");
		
	}

	public static void main(String[] args) throws InstantiationException, IllegalAccessException, SecurityException,
			NoSuchMethodException, IllegalArgumentException, InvocationTargetException {

		Main main = new Main("imputationserver.jar", args);
		main.addTool("version", VersionTool.class);

		main.start();
	}
}

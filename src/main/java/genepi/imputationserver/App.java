package genepi.imputationserver;

import genepi.imputationserver.steps.EstimateAncestryCommand;
import genepi.imputationserver.steps.InputValidationCommand;
import genepi.imputationserver.steps.QualityControlCommand;
import genepi.imputationserver.steps.PrepareTraceCommand;
import picocli.CommandLine;
import picocli.CommandLine.Command;

@Command(name = App.APP, version = App.VERSION)
public class App {

	public static final String APP = "imputationserver-utils";

	public static final String VERSION = "1.5.1";

	public static final String URL = "https://github.com/genepi/imputationserver-utils";

	public static final String COPYRIGHT = "(c) 2023-2024 Lukas Forer and Sebastian Schönherr";

	public static String[] ARGS = new String[0];

	public static void main(String[] args) {

		System.err.println();
		System.err.println(APP + " " + VERSION);
		if (URL != null && !URL.isEmpty()) {
			System.err.println(URL);
		}
		if (COPYRIGHT != null && !COPYRIGHT.isEmpty()) {
			System.err.println(COPYRIGHT);
		}
		System.err.println();

		ARGS = args;

		CommandLine commandLine = new CommandLine(new App());
		commandLine.addSubcommand("validate", new InputValidationCommand());
		commandLine.addSubcommand("run-qc", new QualityControlCommand());
		commandLine.addSubcommand("estimate-ancestry", new EstimateAncestryCommand());
		commandLine.addSubcommand("prepare-trace", new PrepareTraceCommand());
		commandLine.setExecutionStrategy(new CommandLine.RunLast());
		int result = commandLine.execute(args);
		System.exit(result);

	}

}

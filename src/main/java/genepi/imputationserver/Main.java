package genepi.imputationserver;

import java.lang.reflect.InvocationTargetException;

import genepi.base.Toolbox;

public class Main extends Toolbox {

	public static String APP = "imputationserver-utils";
	
	public static String VERSION = "1.1.2";
	
	public Main(String command, String[] args) {
		super(command, args);
		
		System.out.println();
		System.out.println(Main.APP + " " + Main.VERSION);
		System.out.println("https://imputationserver.sph.umich.edu");
		System.out.println("Lukas Forer, Sebastian Schoenherr and Christian Fuchsberger");
		
	}

	public static void main(String[] args) throws InstantiationException, IllegalAccessException, SecurityException,
			NoSuchMethodException, IllegalArgumentException, InvocationTargetException {

		Main main = new Main("imputationserver-utils.jar", args);
		main.start();
	}
}

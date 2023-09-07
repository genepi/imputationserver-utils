package genepi.imputationserver.util.importer;

public class ImporterFactory {

	public static boolean needsImport(String url) {
		return url.startsWith("sftp://") || url.startsWith("http://") || url.startsWith("https://")
				|| url.startsWith("ftp://") || url.startsWith("s3://");
	}

}

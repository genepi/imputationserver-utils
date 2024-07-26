package genepi.imputationserver.util;

public class StringUtils {

    public static String resolveVariable(String text, String variable, String value) {
        return text.replaceAll("\\$" + variable, value).replaceAll("\\$\\{" + variable + "\\}", value);
    }

}

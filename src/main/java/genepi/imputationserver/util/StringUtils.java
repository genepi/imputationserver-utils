package genepi.imputationserver.util;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

public class StringUtils {

    private static  DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.US);

    private static DecimalFormat formatterInteger = new DecimalFormat("###,###.###", symbols);

    private static DecimalFormat formatterDouble = new DecimalFormat("#.00", symbols);

    public static String resolveVariable(String text, String variable, String value) {
        return text.replaceAll("\\$" + variable, value).replaceAll("\\$\\{" + variable + "\\}", value);
    }

    public static String format(int number) {
        return formatterInteger.format(number);
    }

    public static String format(long number) {
        return formatterInteger.format(number);
    }

    public static String format(double number) {
        return formatterDouble.format(number);
    }
}

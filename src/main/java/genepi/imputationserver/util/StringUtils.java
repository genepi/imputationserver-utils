package genepi.imputationserver.util;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

public class StringUtils {

    public static String resolveVariable(String text, String variable, String value) {
        return text.replaceAll("\\$" + variable, value).replaceAll("\\$\\{" + variable + "\\}", value);
    }

    public static DecimalFormat getFormatter() {
        DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.US);
        return new DecimalFormat("###,###.###", symbols);
    }

    public static DecimalFormat getFormatterDecimal() {
        DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.US);
        return new DecimalFormat("#.00", symbols);
    }
}

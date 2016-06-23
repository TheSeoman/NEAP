import java.util.Map;

/**
 * Created by Seoman on 22.06.2016.
 */
public class ExpressionData {
    private String type;
    private double[][] values;
    private Map<String, Integer> idMap;
    private String[] samples;

    public ExpressionData(String type, double[][] values, Map<String, Integer> idMap, String[] samples){
        this.type = type;
        this.values = values;
        this.idMap = idMap;
        this.samples = samples;
    }

    public String getType() {
        return type;
    }

    public double[][] getValues() {
        return values;
    }

    public Map<String, Integer> getIdMap() {
        return idMap;
    }

    public String[] getSamples() {
        return samples;
    }
}

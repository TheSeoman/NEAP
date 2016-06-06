import java.util.Map;

/**
 * Created by seoman on 5/19/16.
 */
public class MotifList {
    private Map<Integer, double[]> motifs;
    private Map<Integer, double[]> pValues;
    private Map<String, Integer> motifIdMap;
    private Map<Integer, Integer> geneIdMap;

    double[][] pccs;
    double[][] zTrans;

    public MotifList(Map<Integer, double[]> motifs, Map<Integer, double[]> pValues, Map<String, Integer> motifIdMap) {
        this.motifIdMap = motifIdMap;
        this.motifs = motifs;
        this.pValues = pValues;
    }

    public double[][] getPccs() {
        return pccs;
    }

    public void setPccs(double[][] pccs) {
        this.pccs = pccs;
    }

    public Map<Integer, Integer> getGeneIdMap() {
        return geneIdMap;
    }

    public void setGeneIdMap(Map<Integer, Integer> geneIdMap) {
        this.geneIdMap = geneIdMap;
    }

    public Map<Integer, double[]> getMotifs() {
        return motifs;
    }

    public Map<Integer, double[]> getpValues() {
        return pValues;
    }

    public Map<String, Integer> getMotifIdMap() {
        return motifIdMap;
    }

    public double[][] getzTrans() {
        return zTrans;
    }

    public void setzTrans(double[][] zTrans) {
        this.zTrans = zTrans;
    }
}

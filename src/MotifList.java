import java.util.Map;

/**
 * Created by seoman on 5/19/16.
 */
public class MotifList {
    private Map<Integer, int[]> motifs;
    private Map<String, Integer> motifIdMap;
    private Map<Integer, Integer> geneIdMap;

    double[][] pccs;
    double[][] zTrans;

    public MotifList(Map<Integer, int[]> motifs, Map<String, Integer> motifIdMap) {
        this.motifIdMap = motifIdMap;
        this.motifs = motifs;
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

    public Map<Integer, int[]> getMotifs() {
        return motifs;
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

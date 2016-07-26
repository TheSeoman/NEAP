import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Created by seoman on 7/24/16.
 */
public class AbberrantNeighbor {
    Network net;
    Set<Integer> genesOfInterest;
    int scoringMethod;
    double edgeWeightThreshhold;
    double cumulCutOff;

    public AbberrantNeighbor(Network net, Set<Integer> genesOfInterest, int scoringMethod, double param1) {
        this.net = net;
        this.genesOfInterest = genesOfInterest;
        this.scoringMethod = scoringMethod;
        this.edgeWeightThreshhold = param1;

    }

    public double calculateAssociation(Map<Integer, Double> abberrantGenes) {
        double totalAssociation = 0;
        int c = 0;
        for (int id1 : abberrantGenes.keySet()) {
            if (net.genMap.containsKey(id1)) {
                double geneAssociation = 0.0;
                for (int id2 : genesOfInterest) {
                    if (net.genMap.containsKey(id2)) {
                        double edge = net.getEdge(id1, id2);
                        geneAssociation = score(geneAssociation, edge, abberrantGenes.get(id1));
                    }
                }
                totalAssociation += geneAssociation;
                c++;
            }
        }
        return totalAssociation / c;
    }

    private double score(double geneAssociation, double edge, double fc) {
        if (scoringMethod == 0) {
            if (edge > edgeWeightThreshhold)
                return 1;
            else
                return geneAssociation;
        }
        if (scoringMethod == 1) {
            if (edge > edgeWeightThreshhold)
                return geneAssociation + (1 - geneAssociation) * edge;
            else
                return geneAssociation;
        }
        if (scoringMethod == 2) {
            return geneAssociation + (1 - geneAssociation) * Math.pow(edge, 2);
        }
        if (scoringMethod == 3){
            if (edge > edgeWeightThreshhold)
                return geneAssociation + edge;
            else
                return geneAssociation;
        }
        if (scoringMethod == 4){
            if (edge > edgeWeightThreshhold)
                return geneAssociation + edge * Math.abs(fc);
            else
                return geneAssociation;
        }
        return 0;
    }

    public void runOnPatientSet(String fcPath, String outPath, Map<String, String> patientMap, String disease) {
        Map<String, Double> associations = new HashMap<>();
        File[] files = new File(fcPath).listFiles();
        for (File file : files) {
            if (file.isFile()) {
                Map<Integer, Double> aberrant = ExpressionParser.getAberrantGenes(file.getAbsolutePath(), 1);
                String patientId = file.getName().substring(0, file.getName().indexOf("."));
                associations.put(patientId, calculateAssociation(aberrant));
            }
        }
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(outPath));
            for(String patient : associations.keySet()){
                String label = patientMap.get(patient).equals(disease) ? "0" : "1";
                out.write(patient + "\t" + associations.get(patient) + "\t" + label + "\n");
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void runOnTCGAData(String fcPath, String outPath) {
        Map<String, Double> associations = new HashMap<>();
        Map<String, String> labels = new HashMap<>();
        String[] tissues = new String[]{"prostate", "thyroid", "lung", "breast", "kidney"};
        for(String tissue : tissues) {
            File[] files = new File(fcPath + "/" + tissue).listFiles();
            for(File file : files) {
                if(!file.getName().equals("total.fc.tsv")) {
                    Map<Integer, Double> aberrant = ExpressionParser.getAberrantGenes(file.getAbsolutePath(), 1);
                    String patientId = file.getName().substring(0, file.getName().indexOf("."));
                    associations.put(patientId, calculateAssociation(aberrant));
                    labels.put(patientId, tissue);
                }
            }

        }
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(outPath));
            for(String patient : associations.keySet()){
                String label = labels.get(patient).equals("prostate") ? "0" : "1";
                out.write(patient + "\t" + associations.get(patient) + "\t" + label + "\n");
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}

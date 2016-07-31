import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Created by seoman on 7/24/16.
 */
public class AbberrantNeighbor {
    Network net;
    Set<Integer> genesOfInterest;
    double edgeWeightThreshold;

    public AbberrantNeighbor(Network net, Set<Integer> genesOfInterest, double param1) {
        this.net = net;
        this.genesOfInterest = genesOfInterest;
        this.edgeWeightThreshold = param1;

    }

    public double calculateAssociation(Map<Integer, Double> aberrantGenes) {
        double totalAssociation = 0;
        int c = 0;
        for (int id1 : aberrantGenes.keySet()) {
            if (net.genMap.containsKey(id1)) {
                double geneAssociation = 0.0;
//                if (genesOfInterest.contains(id1)) {
//                    geneAssociation = 1 * Math.abs(aberrantGenes.get(id1));
//                } else {
                    for (int id2 : genesOfInterest) {
                        if (net.genMap.containsKey(id2)) {
                            double edge = net.getEdge(id1, id2);
                            if (edge >= edgeWeightThreshold)
                                geneAssociation += edge * Math.abs(aberrantGenes.get(id1));
                        }
                    }
//                    if (geneAssociation > 0) {
//                        int neighborCount = 0;
//                        for (int id2 : net.genMap.keySet()) {
//                            if (net.getEdge(id1, id2) >= edgeWeightThreshold) {
//                                neighborCount++;
//                            }
//                        }
//                        geneAssociation /= neighborCount;
//                    }
//                }
                totalAssociation += geneAssociation;
                c++;
            }
        }
        return totalAssociation / c;
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
            for (String patient : associations.keySet()) {
                String label = patientMap.get(patient).equals(disease) ? "0" : "1";
                out.write(patient + "\t" + associations.get(patient) + "\t" + label + "\n");
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void runOnTCGAData(String[] fcPaths, String[] tissues, String disease, String outPath) {
        Map<String, Double> associations = new HashMap<>();
        Map<String, String> labels = new HashMap<>();
        int c = 0;
        for (int i = 0; i < fcPaths.length; i++) {
            List<Map<Integer, Double>> aberrant = ExpressionParser.getTotalAberrantGenes(fcPaths[i], 1);
            for (int j = 0; j < aberrant.size(); j++) {
                String patientId = String.valueOf(c);
                associations.put(patientId, calculateAssociation(aberrant.get(j)));
                labels.put(patientId, tissues[i]);
                c++;
            }
        }
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(outPath));
            for (String patient : associations.keySet()) {
                String label = labels.get(patient).equals(disease) ? "0" : "1";
                out.write(patient + "\t" + associations.get(patient) + "\t" + label + "\n");
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by seoman on 7/26/16.
 */
public class MCSubnet {
    Network net;
    Map<Integer, boolean[]> aberrantVectors;
    Map<Integer, Set<Integer>> neighbors;

    double edgeThreshhold;
    double fcThreshold;
    double consistencyThreshold;

    private String fcDir = "/home/seoman/Documents/NEAP/Prostate Cancer/FoldChange/prostate";
    private String networkPath = "/media/seoman/9CBA3874BA384CD0/Users/User/Documents/Networks/Maria/prostate_gland";
    private String allGenesPath = "/home/seoman/Documents/NEAP/all_genes.txt";

    public MCSubnet(double edgeThreshold, double fcThreshold, double consistencyThreshold){
        this.edgeThreshhold = edgeThreshold;
        this.fcThreshold = fcThreshold;
        this.consistencyThreshold = consistencyThreshold;
        net = NetworkParser.readBinaryNetwork(networkPath, allGenesPath);
        File[] files = new File(fcDir).listFiles();
        Map<String, Map<Integer, Double>> patients = new HashMap<>();
        for(File file : files) {
            if(!file.getName().equals("total.fc.tsv")) {
                Map<Integer, Double> aberrant = ExpressionParser.getAberrantGenes(file.getAbsolutePath(), fcThreshold);
                String patientId = file.getName().substring(0, file.getName().indexOf("."));
                patients.put(patientId, aberrant);
            }
        }

        aberrantVectors = new HashMap<>();

        for(int geneId : net.genMap.keySet()){
            boolean[] vector = new boolean[patients.size()];
            int c = 0;
            for(Map<Integer, Double> patient : patients.values()){
                if(patient.containsKey(geneId)){
                    vector[c] = true;
                }
                c++;
            }
            if(calcConsistency(vector) >= consistencyThreshold){
                aberrantVectors.put(geneId, vector);
            }
        }

        System.out.println("Aberrant genes: " + aberrantVectors.size());

        neighbors = new HashMap<>();

        int connectedGenes = 0;
        int edges = 0;
        for(int geneId : aberrantVectors.keySet()){
            neighbors.put(geneId, new HashSet<Integer>());
            for(int neighbor : aberrantVectors.keySet()){
                if(net.getEdge(geneId, neighbor) >= edgeThreshold ){
                    neighbors.get(geneId).add(neighbor);
                    edges++;
                }
            }
            if(neighbors.get(geneId).size() > 0){
                connectedGenes++;
            }
        }
        System.out.println("Conected aberrant Genes: " + connectedGenes);
        System.out.println("Edges found: " + edges/2);

        findGreedySubnet(6425);
    }

    public void findGreedySubnet(int startGene){
        Set<Integer> currentSet = new HashSet<>();
        currentSet.add(startGene);
        Set<Integer> removedGenes = new HashSet<>();
        double currentConsistency = calcConsistency(aberrantVectors.get(startGene));
        Set<Integer> activeNeighbors = new HashSet<>(neighbors.get(startGene));
        boolean[] currentVector = aberrantVectors.get(startGene);

        while(currentConsistency >= consistencyThreshold && activeNeighbors.size() > 0){
            int bestNeighbor = -1;
            double maxConsistency = -1;
            Set<Integer> toRemove = new HashSet<>();
            for(int neighbor : activeNeighbors){
                boolean[] joinedVector = intersect(currentVector, aberrantVectors.get(neighbor));
                double joinedConsistency = calcConsistency(joinedVector);
                if(joinedConsistency < consistencyThreshold){
                    toRemove.add(neighbor);
                    removedGenes.add(neighbor);
                } else {
                    if(joinedConsistency > maxConsistency){
                        maxConsistency = joinedConsistency;
                        bestNeighbor = neighbor;
                    }
                }
            }
            activeNeighbors.removeAll(toRemove);

            if(bestNeighbor == -1){
                break;
            } else {
                activeNeighbors.remove(bestNeighbor);
                currentSet.add(bestNeighbor);
                currentVector = intersect(currentVector, aberrantVectors.get(bestNeighbor));
                for(int gene : neighbors.get(bestNeighbor)){
                    if(!currentSet.contains(gene) && !removedGenes.contains(gene)){
                        activeNeighbors.add(gene);
                    }
                }
                currentConsistency = calcConsistency(currentVector);
            }
        }
        System.out.println(currentSet);
        System.out.println(currentConsistency);
        System.out.println(currentVector);
    }

    public double calcConsistency(boolean[] vector){
        int aberrant = 0;
        for(int i = 0; i < vector.length; i++){
            if(vector[i])
                aberrant++;
        }
        return (double) aberrant / vector.length;
    }

    public boolean[] intersect(boolean[] v1, boolean[] v2){
        boolean[] v = new boolean[v1.length];
        for(int i = 0; i < v1.length; i++){
            v[i] = v1[i] && v2[i];
        }
        return v;
    }
}

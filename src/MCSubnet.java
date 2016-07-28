import java.io.File;
import java.util.*;

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

    public MCSubnet(double edgeThreshold, double fcThreshold, double consistencyThreshold) {
        this.edgeThreshhold = edgeThreshold;
        this.fcThreshold = fcThreshold;
        this.consistencyThreshold = consistencyThreshold;
        net = NetworkParser.readBinaryNetwork(networkPath, allGenesPath);
        File[] files = new File(fcDir).listFiles();
        Map<String, Map<Integer, Double>> patients = new HashMap<>();
        for (File file : files) {
            if (!file.getName().equals("total.fc.tsv")) {
                Map<Integer, Double> aberrant = ExpressionParser.getAberrantGenes(file.getAbsolutePath(), fcThreshold);
                String patientId = file.getName().substring(0, file.getName().indexOf("."));
                patients.put(patientId, aberrant);
            }
        }

        aberrantVectors = new HashMap<>();

        for (int geneId : net.genMap.keySet()) {
            boolean[] vector = new boolean[patients.size()];
            int c = 0;
            for (Map<Integer, Double> patient : patients.values()) {
                if (patient.containsKey(geneId)) {
                    vector[c] = true;
                }
                c++;
            }
            if (calcConsistency(vector) >= consistencyThreshold) {
                aberrantVectors.put(geneId, vector);
            }
        }

        System.out.println("Aberrant genes: " + aberrantVectors.size());

        neighbors = new HashMap<>();

        int connectedGenes = 0;
        int edges = 0;
        for (int geneId : aberrantVectors.keySet()) {
            neighbors.put(geneId, new HashSet<Integer>());
            for (int neighbor : aberrantVectors.keySet()) {
                if (net.getEdge(geneId, neighbor) >= edgeThreshold) {
                    neighbors.get(geneId).add(neighbor);
                    edges++;
                }
            }
            if (neighbors.get(geneId).size() > 0) {
                connectedGenes++;
            }
        }
        System.out.println("Connected aberrant Genes: " + connectedGenes);
        System.out.println("Edges found: " + edges / 2);

        findTotalKGreedySubnet(2);
    }

    public void findGreedySubnet(int startGene) {
        Set<Integer> currentSet = new HashSet<>();
        currentSet.add(startGene);
        Set<Integer> removedGenes = new HashSet<>();
        double currentConsistency = calcConsistency(aberrantVectors.get(startGene));
        Set<Integer> activeNeighbors = new HashSet<>(neighbors.get(startGene));
        boolean[] currentVector = aberrantVectors.get(startGene);

        while (currentConsistency >= consistencyThreshold && activeNeighbors.size() > 0) {
            int bestNeighbor = -1;
            double maxConsistency = -1;
            Set<Integer> toRemove = new HashSet<>();
            for (int neighbor : activeNeighbors) {
                boolean[] joinedVector = intersect(currentVector, aberrantVectors.get(neighbor));
                double joinedConsistency = calcConsistency(joinedVector);
                if (joinedConsistency < consistencyThreshold) {
                    toRemove.add(neighbor);
                    removedGenes.add(neighbor);
                } else {
                    if (joinedConsistency > maxConsistency) {
                        maxConsistency = joinedConsistency;
                        bestNeighbor = neighbor;
                    }
                }
            }
            activeNeighbors.removeAll(toRemove);

            if (bestNeighbor == -1) {
                break;
            } else {
                activeNeighbors.remove(bestNeighbor);
                currentSet.add(bestNeighbor);
                currentVector = intersect(currentVector, aberrantVectors.get(bestNeighbor));
                for (int gene : neighbors.get(bestNeighbor)) {
                    if (!currentSet.contains(gene) && !removedGenes.contains(gene)) {
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

    public void findTotalKGreedySubnet(int k){
        List<Set<Integer>> startSets = new ArrayList<>();
        for(int gene : aberrantVectors.keySet()){
            Set<Integer> temp = new HashSet<>();
            temp.add(gene);
            startSets.add(temp);
        }
        findKGreedySubnet(startSets, k);
    }

    public void findSeedKGreedySubnet(Set<Integer> seed, int k){
        List<Set<Integer>> startSets = new ArrayList<>();
        startSets.add(seed);
        findKGreedySubnet(startSets, k);
    }

    public void findKGreedySubnet(List<Set<Integer>> startSets, int k) {
        List<Set<Integer>> currentSets = startSets;
        List<Set<Integer>> activeNeighbors = new ArrayList<>();
        List<boolean[]> currentVectors = new ArrayList<>();
        List<Set<Integer>> removedGenes = new ArrayList<>();
        for(Set<Integer> set : startSets){
            Set<Integer> activeNeighborSet = new HashSet<>();
            boolean[] vector = null;
            for(int gene : set){
                activeNeighborSet.addAll(neighbors.get(gene));
                if(vector == null){
                    vector = aberrantVectors.get(gene);
                } else {
                    vector = intersect(vector, aberrantVectors.get(gene));
                }
            }
            activeNeighbors.add(activeNeighborSet);
            currentVectors.add(vector);
            removedGenes.add(new HashSet<Integer>());
        }

        List<Set<Integer>> newCurrentSets;
        List<Set<Integer>> newActiveNeighbors;
        List<boolean[]> newCurrentVectors;
        List<Set<Integer>> newRemovedGenes;

        List<Set<Integer>> terminatedSets = new ArrayList<>();
        List<boolean[]> terminatedVectors = new ArrayList<>();

        while (currentSets.size() > 0) {
            newCurrentSets = new ArrayList<>();
            newActiveNeighbors = new ArrayList<>();
            newCurrentVectors = new ArrayList<>();
            newRemovedGenes = new ArrayList<>();

            for (int i = 0; i < currentSets.size(); i++) {


                Map<Integer, boolean[]> joinedVectors = new HashMap<>();
                Map<Integer, Double> joinedConsistencies = new HashMap<>();
                List<Integer> order = new LinkedList<>();
                Set<Integer> toRemove = new HashSet<>();
                for (int neighbor : activeNeighbors.get(i)) {
                    boolean[] joinedVector = intersect(currentVectors.get(i), aberrantVectors.get(neighbor));
                    double joinedConsistency = calcConsistency(joinedVector);
                    if (joinedConsistency < consistencyThreshold) {
                        toRemove.add(neighbor);
                        removedGenes.get(i).add(neighbor);
                    } else {
                        joinedVectors.put(neighbor, joinedVector);
                        joinedConsistencies.put(neighbor, joinedConsistency);
                        boolean added = false;
                        for (int j = 0; j < order.size(); j++) {
                            if (joinedConsistency > joinedConsistencies.get(order.get(j))) {
                                order.add(j, neighbor);
                                added = true;
                                break;
                            }
                        }
                        if (!added) {
                            order.add(neighbor);
                        }
                    }
                }
                activeNeighbors.get(i).removeAll(toRemove);

                if (order.size() == 0) {
                    System.out.println("No more genes can be added to: " + currentSets.get(i));
                    terminatedSets.add(currentSets.get(i));
                    terminatedVectors.add(currentVectors.get(i));
                } else {
                    for (int j = 0; j < Math.min(order.size(), k); j++) {
                        int geneToAdd = order.get(j);
                        System.out.println(geneToAdd + "\t" + currentSets.get(i));
                        Set<Integer> newCurrentSet = new HashSet<>(currentSets.get(i));
                        newCurrentSet.add(geneToAdd);
                        boolean duplicate = false;
                        for (Set<Integer> currentSet : newCurrentSets) {
                            if (currentSet.containsAll(newCurrentSet)) {
                                System.out.println("Set " + currentSet + " already in newCurrentSets");
                                duplicate = true;
                                break;
                            }
                        }
                        if (!duplicate) {
                            newCurrentSets.add(newCurrentSet);
                            newCurrentVectors.add(joinedVectors.get(geneToAdd));
                            Set<Integer> newActiveNeighbor = new HashSet<>(activeNeighbors.get(i));
                            newActiveNeighbor.remove(geneToAdd);
                            for (int gene : neighbors.get(geneToAdd)) {
                                if (!newCurrentSet.contains(gene) && !removedGenes.get(i).contains(gene)) {
                                    newActiveNeighbor.add(gene);
                                }
                            }
                            newActiveNeighbors.add(newActiveNeighbor);
                            Set<Integer> newRemovedGeneSet = new HashSet<>(removedGenes.get(i));
                            newRemovedGenes.add(newRemovedGeneSet);
                        }
                    }
                }
            }
            currentSets = newCurrentSets;
            currentVectors = newCurrentVectors;
            activeNeighbors = newActiveNeighbors;
            removedGenes = newRemovedGenes;
        }
        System.out.println(terminatedSets.size());
    }

    private double calcConsistency(boolean[] vector) {
        int aberrant = 0;
        for (int i = 0; i < vector.length; i++) {
            if (vector[i])
                aberrant++;
        }
        return (double) aberrant / vector.length;
    }

    private boolean[] intersect(boolean[] v1, boolean[] v2) {
        boolean[] v = new boolean[v1.length];
        for (int i = 0; i < v1.length; i++) {
            v[i] = v1[i] && v2[i];
        }
        return v;
    }

}

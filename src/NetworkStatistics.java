import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created by schmidtju on 25.04.16.
 */
public class NetworkStatistics {

    public static void main(String[] args) {
        if (args.length != 3) {
            System.out.println("Usage: networkstats.jar <network file><weight threshold><output dir>");
        } else {
            String path = args[0];
            Double threshold = Double.parseDouble(args[1]);
            String dir = args[2];
            if (!dir.endsWith("/")) {
                dir += "/";
            }
            System.out.println("Parsing " + path + "...");
            Net net = NetworkParser.parseNetworkFileGZ(path, threshold, 26000);
            String[] split = path.split("/");
            String tissue = split[split.length - 1];
            tissue = tissue.substring(0, tissue.length() - 3);
            System.out.println("Calculating statistics for " + tissue + "...");
//            calcBasicStats(net, tissue, dir);
            saveIds(net, tissue, dir);
        }
    }

    public static void calcBasicStats(Net network, String tissue, String dir) {
        int nodes = 0;
        int edges = 0;
        byte[][] net = network.getNet();
        Map<Integer, Integer> idMap = network.getIdMap();
        Map<Integer, Integer> hist = new HashMap<Integer, Integer>();
        for (int i = 0; i < net.length; i++) {
            boolean present = false;
            for (int j = 0; j < net.length; j++) {
                if (net[i][j] == 1) {
                    edges++;
                    present = true;
                }
            }
            if (present) {
                nodes++;
                if (!hist.containsKey(edges)) {
                    hist.put(edges, 1);
                } else {
                    hist.put(edges, hist.get(edges) + 1);
                }
            }
        }
        double graphDensity = edges / (nodes * (nodes - 1.0));

        try {
            BufferedWriter out_stats = new BufferedWriter(new FileWriter(dir + tissue + "_stats"));
            out_stats.write("#nodes:\t" + nodes + "\n");
            out_stats.write("#edges:\t" + edges + "\n");
            out_stats.write("Graph Density:\t" + graphDensity + "\n");
            out_stats.close();
            BufferedWriter out_hist = new BufferedWriter(new FileWriter(dir + tissue + "_hist"));
            for (int degree : hist.keySet()) {
                out_hist.write(degree + "\t" + hist.get(degree) + "\n");
            }
            out_hist.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void saveIds(Net net, String tissue, String dir) {
        Map<Integer, Integer> idMap = net.getIdMap();
        try {
            BufferedWriter out_ids = new BufferedWriter(new FileWriter(dir + tissue + "_ids"));
            for (int gene : idMap.keySet()) {
                out_ids.write(gene + "\n");
            }
            out_ids.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void calcClusteringCoefficient(byte[][] net) {
        double[] localCC = new double[net.length];
        double globalCC = 0;
        for (int i = 0; i < net.length; i++) {
            localCC[i] = 0;
        }
        for (int i = 0; i < net.length; i++) {
            int edges = 0;
            for (int j = i; j < net.length; j++) {
                for (int h = j + 1; h < net.length; h++) {
                    if (net[i][j] == 1)
                        edges++;
                    if (net[i][h] == 1 && net[j][h] == 1) {
                        localCC[i] += 1;
                        localCC[j] += 1;
                    }
                }
            }
            localCC[i] /= edges;
            globalCC += localCC[i];
        }
        System.out.println("Global Clustering Coefficient: " + globalCC);
    }

    public static void calcClusteringCoefficient(Map<Integer, Map<Integer, Double>> net) {
        Map<Integer, Double> localCC = new HashMap<Integer, Double>(net.size());
        double globalCC = 0;
        for (int gene1 : net.keySet()) {
            localCC.put(gene1, 0.0);
            for (int gene2 : net.get(gene1).keySet()) {
                for (int gene3 : net.get(gene1).keySet()) {
                    if (net.containsKey(gene2) && net.get(gene2).containsKey(gene3))
                        localCC.put(gene1, localCC.get(gene1) + 1);
                }
            }
            localCC.put(gene1, localCC.get(gene1) / (net.get(gene1).size() * (net.get(gene1).size() - 1)));
            globalCC += localCC.get(gene1);
        }
        System.out.println("Global Clustering Coefficient: " + globalCC);
    }

    public static double calcTStatistic(Network network, Set<Integer> geneSet) {
        int cSet = 0;
        int cOut = 0;
        List<Double> inSet = new LinkedList<>();
        List<Double> adjacent = new LinkedList<>();
        for (int id1 : geneSet) {
            if(!network.genMap.containsKey(id1)){
                System.out.println("Entrez Gene Id: " + id1 + " not in network.");
                continue;
            }
            for (int id2 : network.genMap.keySet()) {
                if (geneSet.contains(id2)) {
                    if (id1 > id2) {
                        inSet.add(network.getEdge(id1, id2));
                    }
                } else {
                    if (id1 > id2) {
                        adjacent.add(network.getEdge(id1, id2));
                    }
                }
            }
        }

        int nw = inSet.size();
        int nb = adjacent.size();

        double Xw = 0; //mean within
        for(double d : inSet){
            Xw += d;
        }
        Xw /= nw;

        double Xb = 0; //mean adjacent
        for(double d : adjacent){
            Xb += d;
        }
        Xb /= nb;

        double sw = 0; //standard deviation within
        for(double d : inSet){
            sw += Math.pow(d - Xw, 2);
        }
        sw = Math.sqrt(sw/nw);

        double sb = 0; //standard deviation adjacent
        for(double d : inSet){
            sb += Math.pow(d - Xb, 2);
        }
        sb = Math.sqrt(sw/nb);

        double sx = Math.sqrt(Math.pow(sw, 2)/nw + Math.pow(sb, 2)/nb);

        double t = (Xw - Xb) / sx;

        return t;
    }

    public static double calcZStatistic(Network network, Set<Integer> geneSet, int iterations) {
        geneSet.retainAll(network.genMap.keySet());
        double t_up = NetworkStatistics.calcTStatistic(network, geneSet);
        double[] randomTs = new double[iterations];
        double randomSum = 0;
        for(int i = 0; i < iterations; i++){
            Set<Integer> genesRand = GeneIdParser.generateRandomGeneSet(network.genMap.keySet(), geneSet.size());
            double randomT = NetworkStatistics.calcTStatistic(network, genesRand);
            randomSum += randomT;
            randomTs[i] = randomT;
        }
        double randomMean = randomSum/iterations;
        double randomDev = 0;
        for(double randomT : randomTs){
            randomDev += Math.pow(randomT - randomMean, 2);
        }
        randomDev = Math.sqrt(randomDev/iterations);
        System.out.println("T: " + t_up);
        System.out.println("randomMean: " + randomMean);
        System.out.println("randomStandardDev: " + randomDev);

        double z = (t_up - randomMean)/randomDev;
        System.out.println("Z: " + z);

        return z;
    }

    public static void calculateMcSubnet(Network network, HashMap<Integer, Integer> aberrantCount, double edgeWeightThreshhold, int aberrantThreshhold){
        Set<Integer> visited = new HashSet<>();
        List<Set<Integer>> subsets = new ArrayList<Set<Integer>>();

        for(int id1 : network.genMap.keySet()){
            if(visited.contains(id1)){
                continue;
            }
            Set<Integer> subset = new HashSet<>();
            subset.add(id1);
            for(int id2 : network.genMap.keySet()){
                if(visited.contains(id2)){
                    break;
                }
                if(network.getEdge(id1, id2) >= edgeWeightThreshhold && aberrantCount.get(id2) >= aberrantThreshhold){

                }
            }

        }

    }

}

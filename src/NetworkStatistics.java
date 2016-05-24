import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

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
            if(!dir.endsWith("/")){
                dir += "/";
            }
            System.out.println("Parsing " + path + "...");
            Network network = NetworkParser.parseNetworkFileGZ(path, threshold, 26000);
            String[] split = path.split("/");
            String tissue = split[split.length - 1];
            tissue = tissue.substring(0, tissue.length() - 3);
            System.out.println("Calculating statistics for " + tissue + "...");
//            calcBasicStats(network, tissue, dir);
            saveIds(network, tissue, dir);
        }
    }

    public static void calcBasicStats(Map<Integer, Map<Integer, Double>> net, String dir, String tissue) {
        int nodes = net.size();
        int edges = 0;
        Map<Integer, Integer> hist = new HashMap<Integer, Integer>();
        for (int gene : net.keySet()) {
            edges += net.get(gene).size();
            if (!hist.containsKey(net.get(gene).size()))
                hist.put(net.get(gene).size(), 1);
            else
                hist.put(net.get(gene).size(), hist.get(net.get(gene).size()) + 1);
        }
        double graphDensity = edges / (nodes * (nodes - 1.0) / 2);

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

    public static void calcBasicStats(Network network, String tissue, String dir) {
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

    public static void saveIds(Network network, String tissue, String dir) {
        Map<Integer, Integer> idMap = network.getIdMap();
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
                    if (net[i][h] == 1 && net[j][h] == 1){
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

}

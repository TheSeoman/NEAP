import java.io.*;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;
/**
 * Created by schmidtju on 25.04.16.
 */
public class NetworkParser {
    public static Net parseNetworkFile(String path, double threshold, int size) {
        HashMap<Integer, Integer> idMap = new HashMap<Integer, Integer>(size);
        byte[][] net = new byte[size][size];
        int id1 = 0, id2 = 0, idcount = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            while((line = br.readLine()) != null){
                String[] split = line.split("\t");
                if(!split[0].matches("^\\d+$"))
                    continue;
                int gene1 = Integer.parseInt(split[0]);
                if(idMap.containsKey(gene1)){
                    id1 = idMap.get(gene1);
                } else {
                    id1 = idcount;
                    idcount++;
                    idMap.put(gene1, id1);
                }
                if(!split[1].matches("^\\d+$"))
                    continue;
                int gene2 = Integer.parseInt(split[1]);
                if(idMap.containsKey(gene2)){
                    id2 = idMap.get(gene2);
                } else {
                    id2 = idcount;
                    idcount++;
                    idMap.put(gene2, id2);
                }
                if(split.length > 2) {
                    double weight = Double.parseDouble(split[2]);
                    if (weight > threshold) {
                        net[id1][id2] = 1;
                        net[id2][id1] = 1;
                    }
                } else {
                    net[id1][id2] = 1;
                    net[id2][id1] = 1;
                }

            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return new Net(net, idMap);
    }

    public static Network readBinaryNetwork(String networkPath, String allGenePath){
        try {
            return new Network(networkPath, allGenePath);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    public static Net parseNetworkFileGZ(String path, double threshold, int size) {
        HashMap<Integer, Integer> idMap = new HashMap<Integer, Integer>(size);
        byte[][] net = new byte[size][size];
        int id1 = 0, id2 = 0, idcount = 0;
        try {
            InputStream fileStream = new FileInputStream(path);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
            BufferedReader br = new BufferedReader(decoder);
            String line;
            while((line = br.readLine()) != null){
                String[] split = line.split("\t");
                int gene1 = Integer.parseInt(split[0]);
                if(idMap.containsKey(gene1)){
                    id1 = idMap.get(gene1);
                } else {
                    id1 = idcount;
                    idcount++;
                    idMap.put(gene1, id1);
                }
                int gene2 = Integer.parseInt(split[1]);
                if(idMap.containsKey(gene2)){
                    id2 = idMap.get(gene2);
                } else {
                    id2 = idcount;
                    idcount++;
                    idMap.put(gene2, id2);
                }
                double weight = Double.parseDouble(split[2]);
                if(weight > threshold){
                    net[id1][id2] = 1;
                    net[id2][id1] = 1;
                }

            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return new Net(net, idMap);
    }

}

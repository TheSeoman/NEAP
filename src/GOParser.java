import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by Seoman on 06.06.2016.
 */
public class GOParser {

    public static Map<Integer, Set<Integer>> readPositiveNegativePairs(String path) {
        Map<Integer, Set<Integer>> pairs = new HashMap<>();
        try {
            String line;
            BufferedReader br = new BufferedReader(new FileReader(path));
            while ((line = br.readLine()) != null) {
                String[] split = line.split("\t");
                if(split.length != 2 || split[0].equals("null") || split[1].equals("null"))
                    continue;
                int id1 = Integer.parseInt(split[0]);
                int id2 = Integer.parseInt(split[1]);
                if(!pairs.containsKey(id1)){
                    pairs.put(id1, new HashSet<Integer>(id2));
                } else {
                    pairs.get(id1).add(id2);
                }
            }
            br.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
        return pairs;
    }
}
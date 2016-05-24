import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by seoman on 5/24/16.
 */
public class MotifParser {
    public static MotifList readMotifs(String path, double threshold, int motifCount) {
        Map<String, Integer> idMap = new HashMap<String, Integer>();
        Map<Integer, int[]> motifs = new HashMap<Integer, int[]>();
        int c = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line = br.readLine();
            while ((line = br.readLine()) != null) {
                String[] split = line.split("\t");
                int geneId = Integer.parseInt(split[1]);
                if (!motifs.containsKey(geneId)) {
                    motifs.put(geneId, new int[motifCount]);
                }
                if (!idMap.containsKey(split[0])) {
                    idMap.put(split[0], c);
                    c++;
                }
                if (Double.parseDouble(split[6]) < 0.001)
                    motifs.get(geneId)[idMap.get(split[0])] = 1;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return new MotifList(motifs, idMap);
    }
}

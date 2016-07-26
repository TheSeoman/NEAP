import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by schmidtju on 27.06.16.
 */
public class GeneIdParser {
    public static Set<Integer> parseGeneIds(String path) {
        Set<Integer> ids = new HashSet<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            while ((line = br.readLine()) != null) {
                int id = Integer.parseInt(line.trim());
                ids.add(id);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return ids;
    }

    public static Set<String> parseGeneIds(String path, int col) {
        Set<String> ids = new HashSet<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            while ((line = br.readLine()) != null) {
                String[] split = line.split("\t");
                if(split.length > col) {
                    String id = split[col];
                    ids.add(id);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return ids;
    }

    public static Set<Integer> parseEntrezIds(String path, int col) {
        Set<Integer> ids = new HashSet<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            while ((line = br.readLine()) != null) {
                String[] split = line.split("\t");
                if(split.length > col) {
                    int id = Integer.parseInt(split[col]);
                    ids.add(id);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return ids;
    }

    public static Set<Integer> generateRandomGeneSet(Set<Integer> allIds, int size) {
        List<Integer> all = new ArrayList<>(allIds);
        Set<Integer> ids = new HashSet<>();
        Random rand = new Random(System.currentTimeMillis());

        for (int i = 0; i < size; i++) {
            ids.add(all.remove(rand.nextInt(all.size())));
        }
        return ids;
    }

    public static Map<String, String> parseMappingFile(String mappingPath, int keyCol, int valueCol) {
        Map<String, String> mapping = new HashMap<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(mappingPath));
            String line;
            String[] split;
            int expectedCols = Math.max(keyCol, valueCol);
            while ((line = br.readLine()) != null) {
                split = line.split("\t");
                if (split.length > expectedCols && !split[valueCol].equals("") && !split[keyCol].equals("")) {
                    if (!mapping.containsKey(split[keyCol])) {
                        mapping.put(split[keyCol], split[valueCol]);
                    } else {
                        System.out.println("Ignore multiple occurrence of key: " + split[keyCol]);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return mapping;
    }
}
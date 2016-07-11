import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by schmidtju on 27.06.16.
 */
public class GeneIdParser {
    public static Set<Integer> readGeneIds(String path) {
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

    public static Set<Integer> generateRandomGeneSet(Set<Integer> allIds, int size) {
        List<Integer> all = new ArrayList<>(allIds);
        Set<Integer> ids = new HashSet<>();
        Random rand = new Random(System.currentTimeMillis());

        for (int i = 0; i < size; i++) {
            ids.add(all.remove(rand.nextInt(all.size())));
        }
        return ids;
    }

    public static Map<String, String> parseMappingFile(String mappingPath, int entrezCol, int ensemblCol) {
        Map<String, String> ensembl2Entrez = new HashMap<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(mappingPath));
            String line;
            String[] split;
            int expectedCols = Math.max(entrezCol, ensemblCol);
            while ((line = br.readLine()) != null) {
                split = line.split("\t");
                if (split.length > expectedCols && !split[ensemblCol].equals("") && !split[entrezCol].equals("")) {
                    if (!ensembl2Entrez.containsKey(split[ensemblCol])) {
                        ensembl2Entrez.put(split[ensemblCol], split[entrezCol]);
                    } else {
                        System.out.println("Ignore multiple occurrence of: " + split[ensemblCol]);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return ensembl2Entrez;
    }
}
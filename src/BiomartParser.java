import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * Created by schmidtju on 09.05.16.
 */
public class BiomartParser {
    public static void changeToEntrezId(String pathIdMapping, String pathSeqFile, String pathOut) {
        HashMap<String, String> idMap = new HashMap<String, String>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(pathIdMapping));
            String line;
            while ((line = br.readLine()) != null) {
                String[] split = line.split("\t");
                String entrezId = split[1];
                String ensemblId = split[0];
                idMap.put(ensemblId, entrezId);
            }
            br.close();
            br = new BufferedReader(new FileReader(pathSeqFile));
            BufferedWriter out = new BufferedWriter(new FileWriter(pathOut));
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    String ensemblId = line.substring(1, 16);
                    String entrezId = idMap.get(ensemblId);
                    out.write(">" + entrezId + "\n");
                } else {
                    out.write(line + "\n");
                }

            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}

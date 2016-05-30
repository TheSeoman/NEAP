import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
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

    public static void shuffleSequences(String pathIn, String pathOut){
        try {
            String line;
            BufferedReader br = new BufferedReader(new FileReader(pathIn));
            BufferedWriter out = new BufferedWriter(new FileWriter(pathOut));
            String seq = "";
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    if(seq != ""){
                        out.write(line + "_shuffled\n");
                        seq = shuffle(seq);
                        int i;
                        for(i = 0; i < seq.length() - 60; i += 60){
                            out.write(seq.substring(i, i + 60) + "\n");
                        }
                        out.write(seq.substring(i, seq.length()) + "\n");
                    }
                    seq = "";
                } else {
                    seq += line;
                }
            }
            br.close();
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static String shuffle(String input){
        List<Character> characters = new ArrayList<Character>();
        for(char c:input.toCharArray()){
            characters.add(c);
        }
        StringBuilder output = new StringBuilder(input.length());
        while(characters.size()!=0){
            int randPicker = (int)(Math.random()*characters.size());
            output.append(characters.remove(randPicker));
        }
        return output.toString();
    }
}

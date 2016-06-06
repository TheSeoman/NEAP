import java.io.*;
import java.util.*;

/**
 * Created by seoman on 5/24/16.
 */
public class MotifParser {
    public static MotifList readMotifs(String path, double threshold, int motifCount) {
        Map<String, Integer> idMap = new HashMap<String, Integer>();
        Map<Integer, double[]> motifs = new HashMap<Integer, double[]>();
        Map<Integer, double[]> pValues = new HashMap<Integer, double[]>();

        int c = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line = br.readLine();
            while ((line = br.readLine()) != null) {
                String[] split = line.split("\t");
                int geneId = Integer.parseInt(split[1]);
                if (!motifs.containsKey(geneId)) {
                    motifs.put(geneId, new double[motifCount]);
                    pValues.put(geneId, new double[motifCount]);
                }
                if (!idMap.containsKey(split[0])) {
                    idMap.put(split[0], c);
                    c++;
                }
                double pValue = Double.parseDouble(split[6]);
                pValues.get(geneId)[idMap.get(split[0])] = pValue;
                if (pValue < 0.001)
                    motifs.get(geneId)[idMap.get(split[0])] = 1;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return new MotifList(motifs, pValues, idMap);
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

    public static void filterForGo(Map<Integer, Set<Integer>> positivePairs, Map<Integer, Set<Integer>> negativePairs, String pathMotifs){
        try {
            String line;
            BufferedReader br = new BufferedReader(new FileReader(pathMotifs));
            BufferedWriter outPositive = new BufferedWriter(new FileWriter(pathMotifs.substring(0, pathMotifs.length() - 4) + "_positive.tsv"));
            BufferedWriter outNegative = new BufferedWriter(new FileWriter(pathMotifs.substring(0, pathMotifs.length() - 4) + "_negative.tsv"));
            BufferedWriter outRest = new BufferedWriter(new FileWriter(pathMotifs.substring(0, pathMotifs.length() - 4) + "_rest.tsv"));

            String seq = "";
            while ((line = br.readLine()) != null) {
                String[] split = line.split("\t");
                int id1 = Integer.parseInt(split[0]);
                int id2 = Integer.parseInt(split[1]);
                if((positivePairs.containsKey(id1) && positivePairs.get(id1).contains(id2)) || (positivePairs.containsKey(id2) && positivePairs.get(id2).contains(id1))){
                    outPositive.write(line + "\n");
                } else if((negativePairs.containsKey(id1) && negativePairs.get(id1).contains(id2)) || (negativePairs.containsKey(id2) && negativePairs.get(id2).contains(id1))){
                    outNegative.write(line + "\n");
                } else {
                    outRest.write(line + "\n");
                }
            }
            br.close();
            outPositive.close();
            outNegative.close();
            outRest.close();
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

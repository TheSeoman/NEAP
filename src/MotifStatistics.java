import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by seoman on 5/24/16.
 */
public class MotifStatistics {
    public static double[][] calculatePCCs(MotifList mlist, String pathOut, boolean usepValues, boolean ignoreNull) {
        Map<Integer, double[]> motifs = usepValues ? mlist.getpValues() : mlist.getMotifs();
        Map<Integer, Integer> idMap = new HashMap<>();
        double[][] pccs = new double[motifs.size()][motifs.size()];
        double[][] zTrans;
        int c = 0;
        for (int gene1 : motifs.keySet()) {
            idMap.put(gene1, c);
            c++;
        }
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(pathOut));
            double i = 0.0;
            for (int gene1 : motifs.keySet()) {
                for (int gene2 : motifs.keySet()) {
                    if (gene1 > gene2) {
                        int id1 = idMap.get(gene1);
                        int id2 = idMap.get(gene2);
                        double pcc = ignoreNull ? Statistics.calculatePCCIgnoreZero(motifs.get(gene1), motifs.get(gene2)) : Statistics.calculatePCC(motifs.get(gene1), motifs.get(gene2));
                        pccs[id1][id2] = pcc;
                        pccs[id2][id1] = pcc;

                        i++;
                        if(i / (motifs.size() * (motifs.size() - 1)/200) % 1 == 0)
                            System.out.println("Progress calculationg pccs: " + (i / (motifs.size() * (motifs.size() - 1)/200)) + "%");
                    }
                }
            }
            zTrans = normalize(pccs);
            System.out.println("Pccs normalized.");
            for (int gene1 : motifs.keySet()) {
                for (int gene2 : motifs.keySet()) {
                    if (gene1 > gene2) {
                        int id1 = idMap.get(gene1);
                        int id2 = idMap.get(gene2);
                        writer.write(gene1 + "\t" + gene2 + "\t" + Math.round(pccs[id1][id2] * 10000) / 10000.0 + "\t" + Math.round(zTrans[id1][id2] * 10000) / 10000.0 + "\n");
                    }
                }
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return pccs;
    }

    public static double[][] normalize(double[][] pccs) {
        double[][] z = new double[pccs.length][];
        for (int i = 0; i < pccs.length; i++) {
            z[i] = Arrays.copyOf(pccs[i], pccs[i].length);
        }
        double mean = 0;
        int c = 0;
        for (int i = 0; i < pccs.length; i++) {
            for (int j = i + 1; j < pccs.length; j++) {
                if(Double.isNaN(z[i][j])){
                    continue;
                }
//                z[i][j] = 0.5 * Math.log((1 + pccs[i][j]) / (1 - pccs[i][j]));
                mean += z[i][j];
                c++;
            }
        }
        mean /= c;
        double sdev = 0;
        for (int i = 0; i < pccs.length; i++) {
            for (int j = i + 1; j < pccs.length; j++) {
                if(Double.isNaN(z[i][j])){
                    continue;
                }
                z[i][j] -= mean;
                sdev += Math.pow(z[i][j], 2);
            }
        }
        sdev = Math.sqrt(sdev/(c - 1));
        for (int i = 0; i < pccs.length; i++) {
            for (int j = i + 1; j < pccs.length; j++) {
                if(Double.isNaN(z[i][j])){
                    continue;
                }
                z[i][j] /= sdev;
                z[j][i] = z[i][j];
            }
        }
        return z;
    }
}

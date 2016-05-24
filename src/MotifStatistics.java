import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by seoman on 5/24/16.
 */
public class MotifStatistics {
    public static double[][] calculatePCCs(MotifList mlist, String pathOut, boolean ignoreNull) {
        Map<Integer, int[]> motifs = mlist.getMotifs();
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
            for (int gene1 : motifs.keySet()) {
                for (int gene2 : motifs.keySet()) {
                    if (gene1 >= gene2) {
                        int id1 = idMap.get(gene1);
                        int id2 = idMap.get(gene2);
                        double pcc = calculatePCC(motifs.get(gene1), motifs.get(gene2));
                        pccs[id1][id2] = pcc;
                        pccs[id2][id1] = pcc;
                    }
                }
            }
            zTrans = zTransform(pccs);
            for (int gene1 : motifs.keySet()) {
                for (int gene2 : motifs.keySet()) {
                    if (gene1 >= gene2) {
                        int id1 = idMap.get(gene1);
                        int id2 = idMap.get(gene2);
                        writer.write(gene1 + "\t" + gene2 + "\t" + Math.round(pccs[id1][id2] * 10000) / 10000.0 + "\t" + zTrans[id1][id2] + "\n");
                    }
                }
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return pccs;
    }

    public static double[][] zTransform(double[][] pccs) {
        double[][] z = new double[pccs.length][pccs.length];
        double mean = 0;
        int c = 0;
        for (int i = 0; i < pccs.length; i++) {
            for (int j = i + 1; j < pccs.length; j++) {
                z[i][j] = 0.5 * Math.log((1 + pccs[i][j]) / (1 - pccs[i][j]));
                mean += z[i][j];
                c++;
            }
        }
        mean /= c;
        double sdev = 0;
        for (int i = 0; i < pccs.length; i++) {
            for (int j = i + 1; j < pccs.length; j++) {
                z[i][j] -= mean;
                sdev += Math.pow(z[i][j], 2);
            }
        }
        sdev /= c - 1;
        for (int i = 0; i < pccs.length; i++) {
            for (int j = i + 1; j < pccs.length; j++) {
                z[i][j] /= sdev;
            }
        }
        return z;
    }

    public static double calculatePCC(int[] v1, int[] v2) {
        double v1sum = 0;
        double v2sum = 0;
        double pcc = 0;
        for (int i = 0; i < v1.length; i++) {
            v1sum += v1[i];
            v2sum += v2[i];
        }
        double v1q = v1sum / v1.length;
        double v2q = v2sum / v1.length;
        double q = 0;
        double d1 = 0;
        double d2 = 0;

        for (int i = 0; i < v1.length; i++) {
            q += (v1[i] - v1q) * (v2[i] - v2q);
            d1 += Math.pow(v1[i] - v1q, 2);
            d2 += Math.pow(v2[i] - v2q, 2);
        }

        pcc = q / (Math.sqrt(d1) * Math.sqrt(d2));
        return pcc;
    }

    public static double calculatePCCIngoreZero(int[] v1, int[] v2) {
        double v1sum = 0;
        double v2sum = 0;
        double pcc = 0;
        int c = 0;
        for (int i = 0; i < v1.length; i++) {
            if (v1[i] != 0 || v2[i] != 0) {
                v1sum += v1[i];
                v2sum += v2[i];
                c++;
            }
        }
        double v1q = v1sum / c;
        double v2q = v2sum / c;
        double q = 0;
        double d1 = 0;
        double d2 = 0;

        for (int i = 0; i < v1.length; i++) {
            if (v1[i] != 0 || v2[i] != 0) {
                q += (v1[i] - v1q) * (v2[i] - v2q);
                d1 += Math.pow(v1[i] - v1q, 2);
                d2 += Math.pow(v2[i] - v2q, 2);
            }
        }

        pcc = q / (Math.sqrt(d1) * Math.sqrt(d2));
        return pcc;
    }

}

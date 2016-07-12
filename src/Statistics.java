import java.util.Arrays;

/**
 * Created by seoman on 6/27/16.
 */
public class Statistics {
    public static double calculatePCC(double[] v1, double[] v2) {
        double v1sum = 0;
        double v2sum = 0;
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

        double pcc = q / (Math.sqrt(d1) * Math.sqrt(d2));
        return pcc;
    }

    public static double calculatePCCIgnoreZero(double[] v1, double[] v2) {
        if (v1.length != v2.length) {
            System.out.println("Different Vector lengths");
            return Double.NaN;
        }
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

    public static double[] normalize(double[] pccs) {
        double[] z = new double[pccs.length];
        z = Arrays.copyOf(pccs, pccs.length);

        double mean = 0;
        int c = 0;
        for (int i = 0; i < pccs.length; i++) {
            if (Double.isNaN(z[i])) {
                continue;
            }
            mean += z[i];
            c++;
        }

        mean /= c;
        double sdev = 0;
        for (int i = 0; i < pccs.length; i++) {
            if (Double.isNaN(z[i])) {
                continue;
            }
            z[i] -= mean;
            sdev += Math.pow(z[i], 2);

        }
        sdev = Math.sqrt(sdev / (c - 1));
        for (int i = 0; i < pccs.length; i++) {
            if (Double.isNaN(z[i])) {
                continue;
            }
            z[i] /= sdev;
        }
        return z;
    }
}

package patientsFC;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Calendar;
import java.util.HashMap;

public class Runner {

    static HashMap<Integer, Integer> nodes = new HashMap<Integer, Integer>();

    public static void main(String[] args) throws IOException {
        String startTime = Calendar.getInstance().getTime().toString().split("\\s+")[3];
        System.out.println("Start at "+startTime);

        String output = "C:/Users/Stefan/Desktop/BLOCKPHASE/NEAP/Prostate Cancer/subnetwork_prostate_01_threshold";
//        String genes = "C:/Users/Stefan/Desktop/BLOCKPHASE/NEAP/Prostate Cancer/all_genes.txt";
//        File file = new File("C:/Users/Stefan/Desktop/BLOCKPHASE/prostate_subnetwork.txt");
        File file = new File("C:/Users/Stefan/Desktop/TissuesGIANT/MALACARD_GENES/subnetwork_genes_weights_01_threshold_updated.txt");
        String genes = "C:/Users/Stefan/Desktop/BLOCKPHASE/aberrant_genes_prostate_subnetwork.txt";

        readInGenes(genes);
        System.out.println("Read in genes successfully...");

        String fileName = file.getName();
        System.out.println(fileName+" read in...");

        new TxtToRAF(file, output, nodes);

        String endTime = Calendar.getInstance().getTime().toString().split("\\s+")[3];
        System.out.println("--------------------------------------------");
        System.out.println("End "+output+" at "+endTime);
        System.out.println("--------------------------------------------");
    }

    private static void readInGenes(String genes) {

        System.out.println("Read in genes in network..");
        InputStream is = null;
        try {
            is = new FileInputStream(genes);
            BufferedReader br = new BufferedReader(new InputStreamReader(is));

            String line = br.readLine();
            int count = 0;

            while ((line = br.readLine()) != null) {
                nodes.put(Integer.parseInt(line.split("\\t")[0]), count++);
            }

//            while (line != null) {
//
//                line = line.trim();
//                nodes.put(Integer.parseInt(line), count);
//                count++;
//
//                line = br.readLine();
//            }

            br.close();
        } catch (IOException e) {
            System.err.print("File can not be read.");
            e.printStackTrace();
        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

}

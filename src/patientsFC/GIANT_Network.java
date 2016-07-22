package patientsFC;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.HashMap;

/**
 * Created by Stefan on 22.07.2016.
 */
public class GIANT_Network {

    public static void main(String[] args) throws IOException {

        String all_genes = "C:/Users/Stefan/Desktop/MaPra/all_genes.txt";
        String pathGIANT = "C:/Users/Stefan/Desktop/TissuesGIANT/RAF/prostate_gland";

        GIANT_Network n = new GIANT_Network(pathGIANT, all_genes);

        System.out.println("!! "+genMap.size());

        String path = "C:/Users/Stefan/Desktop/BLOCKPHASE/";
        File[] patientFolder = new File(path+"NEAP/Prostate Cancer/patient fcs/PATIENT_SET1/").listFiles();
        File pradPatientsSet1 = new File(path+"PRAD_patients_set_1.txt");
        File[] tcgaPatients = new File(path+"NEAP/Prostate Cancer/FoldChange/prostate").listFiles();
        File malacardsFile = new File(path+"NEAP/Prostate Cancer/Malacards/all_unique_prad.txt");
        File allGenesGIANTNetworks = new File(path+"NEAP/Prostate cancer/all_genes.txt");

        double threshold = 1.0;

        AberrantGenes abGenes = new AberrantGenes(allGenesGIANTNetworks, pradPatientsSet1, patientFolder, tcgaPatients, malacardsFile, threshold);

        System.out.println(abGenes.aberrantGeneMap.size() + "...");
    }

    private int geneSize;
    private ByteBuffer[] input;
    static HashMap<Integer, Integer> genMap = new HashMap<Integer, Integer>();

    private boolean isInsideMap = false;

    public GIANT_Network(String path, String genes) throws IOException {
        readInGenes(genes);
        geneSize = genMap.size();
        input = new ByteBuffer[geneSize];
        readRaf(path);
    }

    	/*
     * readInGenes: read in list of genes(nodes) present in Network and generate
	 * hashmap.
	 */

    private static void readInGenes(String genes) {

        System.out.println("Read in genes in network..");
        InputStream is = null;
        try {
            is = new FileInputStream(genes);
            BufferedReader br = new BufferedReader(new InputStreamReader(is));

            String line = br.readLine();
            int count = 0;

            while (line != null) {

                line = line.trim();
                genMap.put(Integer.parseInt(line), count);
                count++;

                line = br.readLine();
            }

            br.close();
        } catch (IOException e) {
            System.err.print("File can not be read.");
            e.printStackTrace();
        }

    }

	/*
	 * readRaf: read in as many bytes as necessary and save them in ByteBuffer
	 * array array length: number of nodes in network buffer lengths:
	 * incrementing=> half matrix
	 */

    private void readRaf(String file) throws IOException {
        RandomAccessFile raf = new RandomAccessFile(file, "r");
        FileChannel channel = raf.getChannel();

        try {
            ByteBuffer bb;
            for (int i = 0; i < geneSize; i++) {
                bb = ByteBuffer.allocate((i+1) * 8);
                channel.read(bb);
                bb.flip();
                input[i] = bb;
                bb.clear();

            }
        } finally {

            channel.close();
            raf.close();
        }

    }

	/*
	 * getEdge: inputs two entrezIDs and returns the corresponding edge weight
	 */

    public double getEdge(Integer a, Integer b) {
        this.setInsideMap(false);
        a = genMap.get(a);
        b = genMap.get(b);
        if (a != null && b != null) {
            this.setInsideMap(true);
            if (a <= b) {
                return (input[b].getDouble((a) * 8));

            } else {

                return (input[a].getDouble((b) * 8));

            }

        }
        return 0.0;
    }

    public boolean isInsideMap() {
        return isInsideMap;
    }

    public void setInsideMap(boolean isInsideMap) {
        this.isInsideMap = isInsideMap;
    }


}

package patientsFC;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by Stefan on 22.07.2016.
 */
public class SubNetwork {

    public static void main(String[] args) throws IOException {

        String all_genes = "C:/Users/Stefan/Desktop/MaPra/all_genes.txt";
        String pathGIANT = "C:/Users/Stefan/Desktop/TissuesGIANT/RAF/prostate_gland";

        SubNetwork n = new SubNetwork(pathGIANT, all_genes, false);

        System.out.println("network successfully read...");

        String path = "C:/Users/Stefan/Desktop/BLOCKPHASE/";
        File[] patientFolder = new File(path+"NEAP/Prostate Cancer/patient fcs/PATIENT_SET1/").listFiles();
        File pradPatientsSet1 = new File(path+"PRAD_patients_set_1.txt");
        File[] tcgaPatients = new File(path+"NEAP/Prostate Cancer/FoldChange/prostate").listFiles();
        File malacardsFile = new File(path+"NEAP/Prostate Cancer/Malacards/all_unique_prad.txt");
        File allGenesGIANTNetworks = new File(path+"NEAP/Prostate cancer/all_genes.txt");

        double threshold = 0.75;

        boolean withPatientSet = false;

        double aberrantThreshold = .7;

        AberrantGenes abGenes = new AberrantGenes(allGenesGIANTNetworks, pradPatientsSet1, patientFolder, tcgaPatients, malacardsFile, withPatientSet, threshold, aberrantThreshold);

        System.out.println(abGenes.aberrantGeneMap.size()+"...");

        int c = 0;
        for (Integer i : abGenes.aberrantGeneMap.keySet()) {
            if (neighborsGene(n, abGenes, i, 0.7).size() > 0) {
                HashSet<Integer> cur = neighborsGene(n, abGenes, i, 0.7);
                System.out.println(i+" :"+cur.size()+cur);
                c++;
            }
        }
        System.out.println("Aberrent Genes with Aberrent neighbors = "+c);
    }

    private static int[] receiveAberrantGene(HashMap<Integer, int[]> map, int gene) {
        return map.get(gene);
    }

    private int geneSize;
    private ByteBuffer[] input;
    static HashMap<Integer, Integer> genMap = new HashMap<Integer, Integer>();
    private int numberOfTCGACases = 0;
    private boolean isInsideMap = false;

    public SubNetwork(String path, String genes, boolean subnetwork) throws IOException {
        readInGenes(genes, subnetwork);
        geneSize = genMap.size();
        input = new ByteBuffer[geneSize];
        readRaf(path);
    }

    private static HashSet<Integer> neighborsGene(SubNetwork n, AberrantGenes abGenes, Integer geneID,
                                                  double edgeWeightThreshold) {
        HashSet<Integer> neighbors = new HashSet<Integer>();
        for (Integer i : abGenes.aberrantGeneMap.keySet()) {
            if (n.getEdge(geneID, i) >= edgeWeightThreshold) {
                neighbors.add(i);
            }
        }
        return neighbors;
    }

    	/*
     * readInGenes: read in list of genes(nodes) present in Network and generate
	 * hashmap.
	 */

    private static void readInGenes(String genes, boolean subNetwork) {

        System.out.println("Read in genes in network..");
        InputStream is = null;
        try {
            is = new FileInputStream(genes);
            BufferedReader br = new BufferedReader(new InputStreamReader(is));

            String line = br.readLine();
            if (subNetwork) //header
                line = br.readLine();

            int count = 0;

            while (line != null) {

                line = line.trim();
                if (subNetwork) {
                    genMap.put(Integer.parseInt(line.split("\\t")[0]), count);
                } else {
                    genMap.put(Integer.parseInt(line), count);
                }
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

    public int getNumberOfTCGACases() {
        return numberOfTCGACases;
    }

    public void setNumberOfTCGACases(int numberOfTCGACases) {
        this.numberOfTCGACases = numberOfTCGACases;
    }


}

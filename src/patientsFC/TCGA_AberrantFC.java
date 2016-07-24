package patientsFC;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

/**
 * Created by Stefan on 24.07.2016.
 * <p>
 * Create a map with aberrant genes
 * </p>
 */
public class TCGA_AberrantFC {

    public static void main(String[] args) throws IOException {
        String path = "C:/Users/Stefan/Desktop/BLOCKPHASE/";

        File[] tcgaFolder = new File(path+"NEAP/Prostate Cancer/FoldChange/prostate/").listFiles();

        File allGenes = new File("C:/Users/Stefan/Desktop/MaPra/all_genes.txt");

        double fcThreshold = 1.0;

        double aberrantThreshold = 0.7;

        getAllGenes(allGenes);

        System.out.println(aberrantMap.size());

        aberrantGenesPatients(tcgaFolder, fcThreshold, aberrantThreshold);
    }

    public static HashMap<Integer, ArrayList<Byte>> aberrantMap = new HashMap<Integer, ArrayList<Byte>>();

    /**
     * create a map with each gene id of GIANT network stored as key
     *
     * @param allGenes file containing all genes in GIANT network
     * @return map containing all genes in GIANT network as key
     * @throws IOException
     */
    private static HashMap<Integer, ArrayList<Byte>> getAllGenes(File allGenes) throws IOException {
        BufferedReader bur = openReader(allGenes);

        String sLine = null;

        while ((sLine = bur.readLine()) != null) {
            ArrayList<Byte> list = new ArrayList<Byte>();
            aberrantMap.put(Integer.parseInt(sLine), list);
        }

        bur.close();

        return aberrantMap;
    }

    /**
     * create aberrant genes map for all genes in network
     *
     * @param folder            folder to tcga prostate cancer cases
     * @param fcThreshold       threshold for log fc value
     * @param aberrantThreshold threshold for aberrant cases in single gene
     * @return map with key geneID and value ArrayList with aberrant tcga where aberrant percentage is above fcThreshold
     * @throws IOException
     */
    private static HashMap<Integer, ArrayList<Byte>> aberrantGenesPatients(File[] folder, double fcThreshold, double aberrantThreshold) throws IOException {
        String sLine = null;
        int fileCur = 0;
        for (File file : folder) {
            BufferedReader bur = openReader(file);
            fileCur++;
            while ((sLine = bur.readLine()) != null) {

                int geneID = Integer.parseInt(sLine.split("\\t")[0]);

                // check if current gene is aberrant
                double fcValue = Double.parseDouble(sLine.split("\\t")[1]);
                int aberrantGene = (fcValue >= fcThreshold || fcValue <= -fcThreshold) ? 1 : 0;

                if (aberrantMap.containsKey(geneID)) {
                    ArrayList<Byte> l = aberrantMap.get(geneID);
                    if (l == null) {
                        ArrayList<Byte> cur = new ArrayList<Byte>();
                        cur.add((byte) aberrantGene);
                        aberrantMap.put(geneID, cur);
                    } else {
                        l.add((byte) aberrantGene);
                        aberrantMap.put(geneID, l);
                    }
                }
            }
            bur.close();

            //add not aberrants to gene if not in file
            for (Integer i : aberrantMap.keySet()) {
                ArrayList<Byte> l = aberrantMap.get(i);
                if (l == null) {
                    ArrayList<Byte> cur = new ArrayList<Byte>();
                    cur.add((byte) 0);
                    aberrantMap.put(i, cur);
                } else if (l.size() < fileCur) {
                    l.add((byte) 0);
                    aberrantMap.put(i, l);
                }
            }
        }

        // remove genes from map that are not "aberrant enough"

        HashMap<Integer, Double> percentageAberrantMap = new HashMap<Integer, Double>();

        for (Integer i : aberrantMap.keySet()) {
            ArrayList<Byte> l = aberrantMap.get(i);
            int aberrantOccurence = Collections.frequency(l, (byte) 1);

            double curPercentageAberrant = (double) aberrantOccurence / (double) l.size();

            percentageAberrantMap.put(i, curPercentageAberrant);
        }

        for (Integer i : percentageAberrantMap.keySet()) {
            if (percentageAberrantMap.get(i) < aberrantThreshold) {
                aberrantMap.remove(i);
            }
        }

        return aberrantMap;
    }


    /**
     * @param file to be opened
     * @return buffered reader for foile
     * @throws IOException
     */
    private static BufferedReader openReader(File file) throws IOException {
        try {
            return new BufferedReader(new FileReader(file));
        } catch (FileNotFoundException eFNF) {
            eFNF.printStackTrace();
        }
        return null;
    }


}

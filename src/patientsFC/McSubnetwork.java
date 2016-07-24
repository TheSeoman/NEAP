package patientsFC;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by Stefan on 23.07.2016.
 */
public class McSubnetwork {

    public static void main(String[] args) throws IOException {
        /* edge weight threshold 0.7, fcThreshold 0.8, aberrantGenesThreshold 0.6 */
        String path = "C:/Users/Stefan/Desktop/BLOCKPHASE/";
        String pathGIANT = path+"NEAP/Prostate Cancer/subnetwork_prostate";
        String allGenes = path+"aberrant_genes_prostate_subnetwork.txt";

//        String allGenes = "C:/Users/Stefan/Desktop/MaPra/all_genes.txt";
//        String pathGIANT = "C:/Users/Stefan/Desktop/TissuesGIANT/RAF/prostate_gland";

        File[] patientFolder = new File(path+"NEAP/Prostate Cancer/patient fcs/PATIENT_SET1/").listFiles();
        File pradPatientsSet1 = new File(path+"PRAD_patients_set_1.txt");
        File[] tcgaPatientsFolder = new File(path+"NEAP/Prostate Cancer/FoldChange/prostate").listFiles();
        File malacardsFile = new File(path+"NEAP/Prostate Cancer/Malacards/all_unique_prad.txt");
        File allGenesGIANTNetworks = new File(path+"NEAP/Prostate cancer/all_genes.txt");

        SubNetwork n = new SubNetwork(pathGIANT, allGenes, true);

        System.out.println("STARTING...");

        double edgeWeightThreshold = 0.7;

        double fcThreshold = 0.6;

        double aberrantGenesThreshold = 0.6;

        int tcgaPatientsNumber = 52;

        int minimalNodesSubnetwork = 5;

        new McSubnetwork(patientFolder, tcgaPatientsFolder, pradPatientsSet1, malacardsFile, false, pathGIANT, allGenes, n, edgeWeightThreshold, fcThreshold, aberrantGenesThreshold, tcgaPatientsNumber, minimalNodesSubnetwork);
    }

    public McSubnetwork(File[] patientFolder, File[] tcgaFolder, File pradPatientsFile, File malacardsFile, boolean withPatientSet,
                        String network, String genes, SubNetwork n, double edgeWeightThreshold, double fcThreshold, double aberrantGenesThreshold, int tcgaPatients, int minimalNodesMCSubnetwork) throws IOException {
        AberrantGenes abGenes = new AberrantGenes(new File(genes), pradPatientsFile, patientFolder, tcgaFolder, malacardsFile, withPatientSet, fcThreshold, aberrantGenesThreshold);

        System.out.println("aberrant genes created...");

        HashMap<Integer, int[]> aberrantMap = abGenes.aberrantGeneMap;

        System.out.println("Aberrant genes = " + aberrantMap.size());

        HashMap<Integer, HashSet<Integer>> neighborsMap = neighbors(n, aberrantMap, edgeWeightThreshold);

        System.out.println("neighbors map created..." + "\t" + neighborsMap.size());

        CreateAberrantPaths.calcUnion(aberrantMap, neighborsMap, minimalNodesMCSubnetwork, aberrantGenesThreshold, tcgaPatients, aberrantGenesThreshold);
    }

    /**
     * @param n                   GIANT or OWN network
     * @param edgeWeightThreshold threshold for edge weight of 2 genes
     * @return neighbors map
     */
    private HashMap<Integer, HashSet<Integer>> neighbors(SubNetwork n, HashMap<Integer, int[]> aberrantGenesMap, double edgeWeightThreshold) {
        HashMap<Integer, HashSet<Integer>> map = new HashMap<Integer, HashSet<Integer>>();
        for (Integer i : aberrantGenesMap.keySet()) {
            //neighbor if edge weight is above used cutoff 0.7
            for (Integer j : aberrantGenesMap.keySet()) {
                if (n.getEdge(i, j) >= edgeWeightThreshold) {
                    fillNeighbors(map, i, j);
                    fillNeighbors(map, j, i);
                }
            }
        }
        return map;
    }

    /**
     * @param map      neighbors map
     * @param gene     current gene
     * @param neighbor current neighbor
     * @return updated neighbors map
     */
    private HashMap<Integer, HashSet<Integer>> fillNeighbors(HashMap<Integer, HashSet<Integer>> map, int gene, int neighbor) {
        HashSet<Integer> set = map.get(gene);

        if (set == null) {
            HashSet<Integer> set1 = new HashSet<Integer>();
            set1.add(neighbor);
            map.put(gene, set1);
        } else {
            set.add(neighbor);
            map.put(gene, set);
        }

        return map;
    }

    /**
     * @param file current file as {@link String}
     * @return opened {@link BufferedReader} for file
     * @throws IOException
     */
    private BufferedReader openReader(String file) throws IOException {
        try {
            return new BufferedReader(new FileReader(new File(file)));
        } catch (FileNotFoundException eFNF) {
            eFNF.printStackTrace();
        }
        return null;
    }

}

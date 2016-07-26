package patientsFC;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

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

        String allGenes = "C:/Users/Stefan/Desktop/MaPra/all_genes.txt";
        String pathGIANT = "C:/Users/Stefan/Desktop/TissuesGIANT/RAF/prostate_gland";

        double fcThreshold = 1.0;

        double aberrantThreshold = .7;

        double intersectionThreshold = .5;

        int sizeOfPatients = 52;

        getAllGenes(new File(allGenes));

        aberrantGenesPatients(tcgaFolder, fcThreshold, aberrantThreshold);

        System.out.println("Aberrant Map: "+aberrantMap.size());

        System.out.println();

        createNeighborsMap(pathGIANT, allGenes, 1.0);

        System.out.println("NeighborsMap:"+neighborsMap.size());

        fillIntersectionList(outputIntersectionList, sizeOfPatients);

//        createIntersection(outputIntersectionList, 1191, intersectionThreshold);

        int startingGene = 1191;

        testIntersection(aberrantMap.get(startingGene), startingGene, intersectionThreshold, sizeOfPatients);
        System.out.println("END: "+outputIntersectionList);
        System.out.println("Visited: "+visitedGenesForIntersection);
        System.out.println();
        System.out.println(m);
        System.out.println("SIZE: "+m.size());
        for (Integer i : m.keySet()) {
            System.out.println(i+" visited");
        }

        System.out.println();
        System.out.println(neighborsMap.size());
    }

    public static HashMap<Integer, ArrayList<Byte>> aberrantMap = new HashMap<Integer, ArrayList<Byte>>();

    public static HashMap<Integer, HashSet<Integer>> neighborsMap = new HashMap<Integer, HashSet<Integer>>();

    public static ArrayList<Byte> outputIntersectionList = new ArrayList<Byte>();
    private static HashSet<Integer> visitedGenesForIntersection = new HashSet<Integer>();
    private static HashMap<Integer, Integer> m = new HashMap<Integer, Integer>();
    private static int n;

    /**
     * @param list
     * @param sizeOfPatients
     * @return list of patient size filled with (byte) 1
     */
    private static ArrayList<Byte> fillIntersectionList(ArrayList<Byte> list, int sizeOfPatients) {
        for (int i = 0; i < sizeOfPatients; i++) {
            list.add((byte) 1);
        }

        return list;
    }

    /**
     * @param out                   arraylist containing current maximal consistent intersection
     * @param gene                  current gene
     * @param intersectionThreshold threshold for "good" intersection - (0.0 : 1.0)
     * @param sizeOfPatients        number of case files (tcga, patients...)
     * @return maximal consistent intersection vector (% of aberrant above intersectionThreshold)
     * @throws IOException
     */
    private static ArrayList<Byte> testIntersection(ArrayList<Byte> out, int gene, double intersectionThreshold, int sizeOfPatients) throws IOException {
        // add gene to visited list
        if (!visitedGenesForIntersection.contains(gene)) {
            visitedGenesForIntersection.add(gene);
        }

        // get current intersection vector
        ArrayList<Byte> currentIntersection = out;

        // create temporary intersection array list
        // used to compare with intersectionThreshold
        ArrayList<Byte> tempIntersection = new ArrayList<Byte>();
        fillIntersectionList(tempIntersection, sizeOfPatients);

        // loop through gene neighbors
        for (Integer i : neighborsMap.get(gene)) {
            // check if gene has not been visited so far
            if (!visitedGenesForIntersection.contains((i))) {
                visitedGenesForIntersection.add(i);

                // get aberrant vector for current gene i
                ArrayList<Byte> currentGeneAberrantList = aberrantMap.get(i);
                System.out.println("OUT: "+gene+"-"+currentIntersection);
                System.out.println("CUR: "+i+"-"+currentGeneAberrantList);
                System.out.print("NEW: [");

                //create intersection of current gene with currentIntersection
                for (int k = 0; k < currentIntersection.size(); k++) {
                    if (currentIntersection.get(k) == (byte) 1 && currentGeneAberrantList.get(k) == (byte) 1) {
                        System.out.print("1, ");
                        tempIntersection.set(k, (byte) 1);
                    } else {
                        tempIntersection.set(k, (byte) 0);
                        System.out.print("0, ");
                    }
                }

                System.out.println();

                // check if intersection is still "good" (above intersectionThreshold)
                int aberrantOccurences = Collections.frequency(tempIntersection, (byte) 1);
                double aberrantPercantage = (double) aberrantOccurences / (double) tempIntersection.size();
                System.out.println(aberrantOccurences+" - "+tempIntersection.size()+" -> "+aberrantPercantage);
                if (aberrantPercantage >= intersectionThreshold) {
                    if (!m.containsKey(gene)) {
                        System.out.println("PUT IN MAP: "+gene+" at position: "+n);
                        m.put(gene, n++);
                    }
                    if (!m.containsKey(i)) {
                        System.out.println("PUT IN MAP: "+i+" at position: "+n);
                        m.put(i, n++);
                    }
                    outputIntersectionList = tempIntersection;
                    System.out.println("NOUT:"+outputIntersectionList);
                    System.out.println();
                    testIntersection(outputIntersectionList, i, intersectionThreshold, sizeOfPatients);
                } else {
                    continue;
                }
            }

        }

        return out;
    }

    /**
     * @param pathGIANT          GIANT or OWN network
     * @param pathAllGenes       all genes in network line by line
     * @param edgeWeighThreshold only neighbors with at least edge weight of threshold taken into account
     * @return neighbors hashmap
     * @throws IOException
     */
    private static HashMap<Integer, HashSet<Integer>> createNeighborsMap(String pathGIANT, String pathAllGenes, double edgeWeighThreshold) throws IOException {
        SubNetwork n = new SubNetwork(pathGIANT, pathAllGenes, false);

        System.out.println("Start creating neighbors map...");

        // edges count
        int edgeCounter = 0;

        // only take direct neighbors into account that have an edge weight of at least edgeWeighThreshold
        for (Integer ida : aberrantMap.keySet()) {
            for (Integer idb : aberrantMap.keySet()) {
                if (n.getEdge(ida, idb) >= edgeWeighThreshold) {
                    edgeCounter++;
                    HashSet<Integer> neighborsIDA = neighborsMap.get(ida);

                    if (neighborsIDA == null) {
                        HashSet<Integer> newList = new HashSet<Integer>();
                        newList.add(idb);
                        neighborsMap.put(ida, newList);
                    } else {
                        neighborsIDA.add(idb);
                        neighborsMap.put(ida, neighborsIDA);
                    }
                }
            }
        }

        System.out.println("Edges counted = "+edgeCounter);

        return neighborsMap;
    }

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

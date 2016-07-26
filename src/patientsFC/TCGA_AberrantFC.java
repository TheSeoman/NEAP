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

        int minimumSizeOfSubnetwork = 8;

//        double outputIntersectionThreshold = .9;

        getAllGenes(new File(allGenes));

        aberrantGenesPatients(tcgaFolder, fcThreshold, aberrantThreshold);

        createNeighborsMap(pathGIANT, allGenes, 1.0);

        System.out.println(aberrantMap.size()+" aberrant map size");
        System.out.println(neighborsMap.size()+" neighbors map size");

        naiveIntersectionGeneOfInterest(aberrantMap.get(1191), 1191, intersectionThreshold, sizeOfPatients);

        System.out.println(outputIntersectionList+"\n"+maximalConsistentMap);

//        findMaximalIntersectionNaive(intersectionThreshold, sizeOfPatients, minimumSizeOfSubnetwork); //, outputIntersectionThreshold);
//
//        for (Integer i : maximalIntersetionsGenesSubnetwork.keySet()) {
//            double curPrecentage = (100 * getPercentageOfAberrantGenesCurrentIntersection(maximalIntersectionsPerGene.get(i)));
//            if (curPrecentage >= 60.0) {
//                System.out.println("MCSubnet of genes number: "+maximalIntersetionsGenesSubnetwork.get(i).size()+" with "+curPrecentage+"% aberrant patients and containing "+maximalIntersetionsGenesSubnetwork.get(i));
//                System.out.println(maximalIntersectionsPerGene.get(i));
//            }
//        }
//
////        System.out.println();
//        System.out.println();
//
//
//        System.out.println(aberrantMap.get(6288)+" - "+6288);
//        System.out.println(aberrantMap.get(2625)+" - "+2625);
//        System.out.println(aberrantMap.get(7857)+" - "+7857);
//        System.out.println(aberrantMap.get(5443)+" - "+5443);
//        System.out.println(aberrantMap.get(5268)+" - "+5268);
//        System.out.println(aberrantMap.get(89780)+" - "+89780);
//        System.out.println(aberrantMap.get(2950)+" - "+2950);
//        System.out.println(aberrantMap.get(10568)+" - "+10568);
//        System.out.println(aberrantMap.get(5579)+" - "+5579);
//        System.out.println(aberrantMap.get(6286)+" - "+6286);
    }


    public static HashMap<Integer, ArrayList<Byte>> aberrantMap = new HashMap<Integer, ArrayList<Byte>>();

    public static HashMap<Integer, HashSet<Integer>> neighborsMap = new HashMap<Integer, HashSet<Integer>>();

    public static ArrayList<Byte> outputIntersectionList = new ArrayList<Byte>();
    private static HashMap<Integer, Boolean> visitedGenesForIntersection = new HashMap<Integer, Boolean>();
    private static HashMap<Integer, Integer> maximalConsistentMap = new HashMap<Integer, Integer>();
    private static int geneVisitedPositionMCSubnet = 0;

    private static HashMap<Integer, ArrayList<Byte>> maximalIntersectionsPerGene = new HashMap<Integer, ArrayList<Byte>>();
    private static HashMap<Integer, ArrayList<Integer>> maximalIntersetionsGenesSubnetwork = new HashMap<Integer, ArrayList<Integer>>();

    private static void findMaximalIntersectionNaive(double intersectionThreshold, int sizeOfPatients, int minimumSizeOfSubnetwork) throws IOException {
        // loop through all aberrant genes
        for (Integer curAberrantGene : neighborsMap.keySet()) {
            // create new list, set, map and set gene visited to 0
            outputIntersectionList = new ArrayList<Byte>();
            visitedGenesForIntersection = new HashMap<Integer, Boolean>();
            maximalConsistentMap = new HashMap<Integer, Integer>();
            geneVisitedPositionMCSubnet = 0;

            fillIntersectionList(outputIntersectionList, sizeOfPatients);

            // get maximal naive intersecction for current gene
            ArrayList<Byte> maximalIntersecionCurrentGene = naiveIntersectionGeneOfInterest(aberrantMap.get(curAberrantGene), curAberrantGene, intersectionThreshold, sizeOfPatients);

            // store genes of maximal naive intersection for new subnetwork
            ArrayList<Integer> currentGeneSubnetworkGenes = new ArrayList<Integer>();
            for (Integer subnetworkGene : maximalConsistentMap.keySet()) {
                currentGeneSubnetworkGenes.add(subnetworkGene);
            }

            // check if current subnetwork contains at least minimumSizeOfSubnetwork nodes
            if (currentGeneSubnetworkGenes.size() >= minimumSizeOfSubnetwork) {
                maximalIntersectionsPerGene.put(curAberrantGene, maximalIntersecionCurrentGene);
                maximalIntersetionsGenesSubnetwork.put(curAberrantGene, currentGeneSubnetworkGenes);
            }
        }
    }

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
    private static ArrayList<Byte> naiveIntersectionGeneOfInterest(ArrayList<Byte> out, int gene, double intersectionThreshold, int sizeOfPatients) throws IOException {
        // add gene to visited list
        if (!visitedGenesForIntersection.containsKey(gene)) {
            visitedGenesForIntersection.put(gene, true);
        }

        // create temporary intersection array list
        // used to compare with intersectionThreshold
        ArrayList<Byte> tempIntersection = new ArrayList<Byte>();
        fillIntersectionList(tempIntersection, sizeOfPatients);

        System.out.println(out+" start");

        // loop through gene neighbors
        for (Integer i : neighborsMap.get(gene)) {
            // check if gene has not been visited so far
//            System.out.print(i+"");
            if (!visitedGenesForIntersection.containsKey((i))) {
                visitedGenesForIntersection.put(i, true);
//                System.out.println(" CUR");

                // get aberrant vector for current gene i
                ArrayList<Byte> currentGeneAberrantList = aberrantMap.get(i);
//                System.out.println("OUT: "+gene+"\t-"+out);
//                System.out.println("CUR: "+i+"\t-"+currentGeneAberrantList);
//                System.out.print("NEW:     \t-[");

                //create intersection of current gene with currentIntersection
                for (int k = 0; k < out.size(); k++) {
                    if (out.get(k) == (byte) 1 && currentGeneAberrantList.get(k) == (byte) 1) {
//                        System.out.print("1, ");
                        tempIntersection.set(k, (byte) 1);
                    } else {
                        tempIntersection.set(k, (byte) 0);
//                        System.out.print("0, ");
                    }
                }

//                System.out.println("\n");

                // check if intersection is still "good" (above intersectionThreshold)
                double aberrantPercantage = getPercentageOfAberrantGenesCurrentIntersection(tempIntersection);

                if (aberrantPercantage >= intersectionThreshold) {
                    if (!maximalConsistentMap.containsKey(gene)) {
                        maximalConsistentMap.put(gene, geneVisitedPositionMCSubnet++);
                    }
                    if (!maximalConsistentMap.containsKey(i)) {
                        maximalConsistentMap.put(i, geneVisitedPositionMCSubnet++);
                    }
//                    System.out.println(currentGeneAberrantList+" cur");
                    System.out.println(tempIntersection+" tmp");
                    naiveIntersectionGeneOfInterest(tempIntersection, i, intersectionThreshold, sizeOfPatients);
                } else {
                    continue;
                }
            } else {
//                System.out.println(" INSIDE");
                continue;
            }

        }

        return outputIntersectionList;
    }

    /**
     * @param tempIntersection represents intersection list
     * @return percentage of aberrant patients in intersecion list
     */
    private static double getPercentageOfAberrantGenesCurrentIntersection(ArrayList<Byte> tempIntersection) {
        return (double) Collections.frequency(tempIntersection, (byte) 1) / tempIntersection.size();
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

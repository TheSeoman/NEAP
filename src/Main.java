
import java.io.File;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by seoman on 5/24/16.
 */

public class Main {
    public static void main(String[] args) {
//        generateRapidMinerInput();


//        runAberrantNeighbor();


        runMCSubnet();


//        ExpressionParser.savePatientsFoldChanges("/home/proj/biosoft/praktikum/neap-ss16/assignments/data/PATIENT_SET1/", "/home/sch/schmidtju/IntellijProjects/NEAP/Prostate Cancer/patient fcs/PATIENT_SET1/");


//        ExpressionParser.testPrediction("/home/seoman/Documents/NEAP/Prostate Cancer/GIANT_PRAD.prediction1", "/home/seoman/Documents/NEAP/Prostate Cancer/set1.patient.map");



//        runDETStat();
//    String experiment = "GDS2545";
//        ExpressionData cancer = ExpressionParser.parseSoftGz("/home/seoman/Documents/NEAP/Prostate Cancer/" + experiment + "_full.soft.gz", new String[]{"prostate tumor"}, "Gene ID", "counts");
//        ExpressionParser.saveExpressionData(cancer, "/home/seoman/Documents/NEAP/Prostate Cancer/FilteredCountFiles/" + experiment + "_cancer");
//        ExpressionData normal = ExpressionParser.parseSoftGz("/home/seoman/Documents/NEAP/Prostate Cancer/" + experiment + "_full.soft.gz", new String[]{"normal prostate tissue"}, "Gene ID", "counts");
//        ExpressionParser.saveExpressionData(normal, "/home/seoman/Documents/NEAP/Prostate Cancer/FilteredCountFiles/" + experiment + "_normal");
//        ExpressionData data = ExpressionParser.parseExpressionData("/home/seoman/Documents/NEAP/Prostate Cancer/FilteredCountFiles/" + experiment + "_cancer", "count");
//        ExpressionStatistics.saveExpressionCorrelations(data, "/home/seoman/Documents/NEAP/Prostate Cancer/Correlations/" + experiment + "_cancer_cor");
//        ExpressionData dat = ExpressionParser.parseExpressionData("/home/seoman/Documents/NEAP/Prostate Cancer/FilteredCountFiles/GDS4114_cancer", "counts");
//        ExpressionStatistics.saveExpressionCorrelations(dat, "/home/seoman/Documents/NEAP/Prostate Cancer/GDS4114_cancer_corr");
        if (args.length < 2) {
            printUsage();
        } else if (args[0].equals("-generateFasta")) {
            try {
                GtfParser.generateFastaFromGff3(args[1], args[2], args[3]);
            } catch (Exception e) {
                printUsage();
            }
        } else if (args[0].equals("-calcCorrelationzScore")) {
            try {
                MotifList mlist = MotifParser.readMotifs(args[1], Double.parseDouble(args[3]), Integer.parseInt(args[4]));
                MotifStatistics.calculatePCCs(mlist, args[2], false, Boolean.parseBoolean(args[5]));
            } catch (Exception e) {
                printUsage();
            }
        } else if (args[0].equals("-filterForGo")) {
            try {
                MotifParser.filterForGo(NetworkParser.parseNetworkFile(args[1], 0, 25000),
                        NetworkParser.parseNetworkFile(args[2], 0, 25000),
                        args[3]);
            } catch (Exception e) {
                printUsage();
            }
        } else if (args[0].equals("-countIntoBins")) {
            try {
                MotifParser.countIntoBins(args[1]);
            } catch (Exception e) {
                printUsage();
            }
        } else {
            printUsage();
        }
    }

    private static void printUsage() {
        System.out.println("Usage:");
        System.out.println("motifs.jar -generateFasta <id file> <gff3 file> <chromosome file dir>");
        System.out.println("motifs.jar -calcCorrelationzScore <fimo output file> <output file> <p-value threshhold> <number of motifs> <ignore zeroes (true|false)>");
        System.out.println("motifs.jar -filterForGo <positive GO file> <negative GO file> <correlation zScore file>");
        System.out.println("motifs.jar -countIntoBins <correlation zScore file>");
    }

    public static void generateRapidMinerInput() {
        String caseFcDir = "/home/seoman/Documents/NEAP/Prostate Cancer/FoldChangePseudo/";
        String[] patientFcDirs = new String[]{"/home/seoman/Documents/NEAP/Prostate Cancer/PatientFoldChangePseudo/PATIENT_SET1",
                "/home/seoman/Documents/NEAP/Prostate Cancer/PatientFoldChangePseudo/PATIENT_SET2",
                "/home/seoman/Documents/NEAP/Prostate Cancer/PatientFoldChangePseudo/PATIENT_SET3"};
//        String featureName = "malacard";
//        Set<Integer> features = GeneIdParser.parseEntrezIds("/home/seoman/Documents/NEAP/Prostate Cancer/Malacards/all_unique_prad.txt",1);
        double[] edgeThresholds = new double[]{0.1};
        double[] consistencyThresholds = new double[]{0.5};
        double fcThreshold = 1.0;
        for(double edgeThreshold : edgeThresholds) {
            for(double consistencyThreshold : consistencyThresholds){
                String featureName = String.valueOf(edgeThreshold) + "_" + String.valueOf(consistencyThreshold) + "_1";
                Set<Integer> features = GeneIdParser.parseEntrezIds("/home/seoman/Documents/NEAP/Prostate Cancer/mcSubnet/" + featureName + ".tsv", 0);
                String disease = "prostate";
                String rapidMinerOutDir = "/home/seoman/Documents/NEAP/Prostate Cancer/RapidMinerInput/";
                for (String patientFcDir : patientFcDirs) {
                    ExpressionParser.saveRapidMinerFoldChanges(caseFcDir, patientFcDir, features, disease,
                            rapidMinerOutDir + "training_" + featureName + ".tsv",
                            rapidMinerOutDir + patientFcDir.substring(patientFcDir.lastIndexOf("/")) + "_" + featureName + ".tsv");
                }
        }
        String featureName = "0.5_0.5_1";

        }
    }

    public static void runAberrantNeighbor() {
        String[] networks = new String[]{"prostate_gland", "thyroid_gland", "lung", "mammary_gland", "kidney"};
        String[] tissues = new String[]{"prostate", "thyroid", "lung", "breast", "kidney"};
        String networkPath = "/media/seoman/9CBA3874BA384CD0/Users/User/Documents/Networks/Maria/";
        String fcDir = "/home/seoman/Documents/NEAP/Prostate Cancer/FoldChangePseudo/";
        String outDir = "/home/seoman/Documents/NEAP/Prostate Cancer/AberrantNeighbors/";
        String[] fcPaths = new String[tissues.length];
        String disease = tissues[1];
        double edgeWeightThreshhold = 0.5;
        for (int i = 0; i < 5; i++) {
            fcPaths[i] = fcDir + tissues[i] + "/total.fc.tsv";
        }
        for (int i = 0; i < 5; i++) {
            Network network = NetworkParser.readBinaryNetwork(networkPath + networks[i], "/home/seoman/Documents/NEAP/all_genes.txt");
            Set<Integer> genesOfInterest = GeneIdParser.parseEntrezIds("/home/seoman/Documents/NEAP/Prostate Cancer/Malacards/thyroid_genes.txt", 0);
            AbberrantNeighbor n = new AbberrantNeighbor(network, genesOfInterest, edgeWeightThreshhold);
            n.runOnTCGAData(fcPaths, tissues, disease, outDir + disease + "_" + networks[i] + "_" + edgeWeightThreshhold + ".tsv");

        }
    }

    public static void runMCSubnet() {
        double[] edgeThresholds = new double[]{0.1};
        double[] consistencyThresholds = new double[]{0.5};
        int k = 2;
        double fcThreshold = 1.0;
        for(double edgeThreshold : edgeThresholds) {
            for(double consistencyThreshold : consistencyThresholds) {
                System.out.println("Genereate mcSubnet: " + edgeThreshold + ", " + fcThreshold + ", " + consistencyThreshold);
                MCSubnet mcs = new MCSubnet(edgeThreshold, fcThreshold, consistencyThreshold);
                mcs.findTotalKGreedySubnet(k);
            }
        }
    }

    public static void generateTCGACountFiles(){
        String[] tissues = new String[]{"prostate", "thyroid", "lung", "breast", "kidney"};
        for(String tissue : tissues) {
            ExpressionParser.mergeTCGACountFiles("/home/seoman/Documents/NEAP/Prostate Cancer/TCGA_" + tissue + "_healthy.json",
                    "/home/seoman/Documents/NEAP/Prostate Cancer/TCGA_" + tissue + "_tumor.json",
                    "/home/seoman/programs/gdc transfer tool/", "/home/seoman/Documents/NEAP/geneId2ensembl",
                    "/home/seoman/Documents/NEAP/Prostate Cancer/TCGAcases/" + tissue + "/",
                    "/home/seoman/Documents/NEAP/Prostate Cancer/FoldChange/" + tissue + "/",
                    "/home/seoman/Documents/NEAP/Prostate Cancer/RapidMinerInput/");
        }
    }

    public static void runDETStat() {
        String[] experiments = new String[]{"GDS2545"};
        String deDir = "/home/sch/schmidtju/IntellijProjects/NEAP/Prostate Cancer/R_out/";

        String[] networks = new String[]{"all_tissues", "prostate_gland", "thyroid_gland", "mammary_gland", "lung"};
        String networkDir = "/home/sch/schmidtju/IntellijProjects/NEAP/networks/";

        int iterations = 100;

        for (String net : networks) {
            Network network = NetworkParser.readBinaryNetwork(networkDir + net, "/home/sch/schmidtju/IntellijProjects/NEAP/all_genes.txt");
            for (String experiment : experiments) {
                System.out.println(net + "\t" + experiment + "\t" + iterations + " iterations");
                Set<Integer> genesUp = GeneIdParser.parseGeneIds(deDir + experiment + "_up.tsv");
                if (genesUp.size() == 0) {
                    System.out.println("No significantly up-regulated Genes.");
                } else {
                    System.out.println("Up-regulated genes: " + genesUp.size());
                    double z_up = NetworkStatistics.calcZStatistic(network, genesUp, iterations);
                }
                Set<Integer> genesDown = GeneIdParser.parseGeneIds(deDir + experiment + "_down.tsv");
                if (genesDown.size() == 0) {
                    System.out.println("No significantly down-regulated Genes.");
                } else {
                    double z_down = NetworkStatistics.calcZStatistic(network, genesDown, iterations);
                }
            }
        }
    }
}

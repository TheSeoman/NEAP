
import java.util.Set;

/**
 * Created by seoman on 5/24/16.
 */

public class Main {
    public static void main(String[] args) {

//        ExpressionParser.saveRapidMinerFoldChanges("/home/seoman/Documents/NEAP/Prostate Cancer/TCGAcases/",
//                "/home/seoman/Documents/NEAP/Prostate Cancer/RapidMinerInput/training.fcs.tsv",
//                GeneIdParser.parseGeneIds("/home/seoman/Documents/NEAP/Prostate Cancer/Malacards/all_unique_prad.txt", 1));
//        ExpressionParser.savePatientsFoldChanges("/home/proj/biosoft/praktikum/neap-ss16/assignments/data/PATIENT_SET1/", "/home/sch/schmidtju/IntellijProjects/NEAP/Prostate Cancer/patient fcs/PATIENT_SET1/");

//        ExpressionParser.saveRapidMinerPatientsFoldChanges("/home/proj/biosoft/praktikum/neap-ss16/assignments/data/PATIENT_SET2/",
//                "/home/sch/schmidtju/IntellijProjects/NEAP/Prostate Cancer/RapidPATIENT_SET2.fc.tsv",
//                GeneIdParser.parseGeneIds("/home/sch/schmidtju/IntellijProjects/NEAP/Prostate Cancer/Malacards/all_unique_prad.txt", 1));

//        ExpressionParser.testPrediction("/home/seoman/Documents/NEAP/Prostate Cancer/GIANT_PRAD.prediction1", "/home/seoman/Documents/NEAP/Prostate Cancer/set1.patient.map");

//        String[] tissues = new String[]{"prostate", "thyroid", "lung", "breast"};
//        for(String tissue : tissues) {
//            ExpressionParser.mergeTCGACountFiles("/home/seoman/Documents/NEAP/Prostate Cancer/TCGA_" + tissue + "_healthy.json",
//                    "/home/seoman/Documents/NEAP/Prostate Cancer/TCGA_" + tissue + "_tumor.json",
//                    "/home/seoman/programs/gdc transfer tool/", "/home/seoman/Documents/NEAP/geneId2ensembl",
//                    "/home/seoman/Documents/NEAP/Prostate Cancer/TCGAcases/" + tissue + "/",
//                    "/home/seoman/Documents/NEAP/Prostate Cancer/FoldChange/" + tissue + "/",
//                    "/home/seoman/Documents/NEAP/Prostate Cancer/RapidMinerInput/");
//        }

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
        if(args.length < 2) {
            printUsage();
        } else if(args[0].equals("-generateFasta")){
            try{
                GtfParser.generateFastaFromGff3(args[1], args[2], args[3]);
            } catch (Exception e){
                printUsage();
            }
        } else if(args[0].equals("-calcCorrelationzScore")){
            try{
                MotifList mlist = MotifParser.readMotifs(args[1], Double.parseDouble(args[3]), Integer.parseInt(args[4]));
                MotifStatistics.calculatePCCs(mlist, args[2], false, Boolean.parseBoolean(args[5]));
            } catch (Exception e){
                printUsage();
            }
        } else if(args[0].equals("-filterForGo")){
            try{
                MotifParser.filterForGo(NetworkParser.parseNetworkFile(args[1], 0, 25000),
                NetworkParser.parseNetworkFile(args[2], 0, 25000),
                args[3]);
            } catch (Exception e){
                printUsage();
            }
        } else if(args[0].equals("-countIntoBins")) {
            try{
                MotifParser.countIntoBins(args[1]);
            } catch (Exception e){
                printUsage();
            }
        } else {
            printUsage();
        }
    }

    private static void printUsage(){
        System.out.println("Usage:");
        System.out.println("motifs.jar -generateFasta <id file> <gff3 file> <chromosome file dir>");
        System.out.println("motifs.jar -calcCorrelationzScore <fimo output file> <output file> <p-value threshhold> <number of motifs> <ignore zeroes (true|false)>");
        System.out.println("motifs.jar -filterForGo <positive GO file> <negative GO file> <correlation zScore file>");
        System.out.println("motifs.jar -countIntoBins <correlation zScore file>");
    }

    public static void runDETStat(){
        String[] experiments = new String[]{"GDS2545"};
        String deDir = "/home/sch/schmidtju/IntellijProjects/NEAP/Prostate Cancer/R_out/";

        String[] networks = new String[]{"all_tissues", "prostate_gland", "thyroid_gland", "mammary_gland", "lung"};
        String networkDir = "/home/sch/schmidtju/IntellijProjects/NEAP/networks/";

        int iterations = 100;

        for(String net : networks){
            Network network = NetworkParser.readBinaryNetwork(networkDir + net, "/home/sch/schmidtju/IntellijProjects/NEAP/all_genes.txt");
            for(String experiment : experiments){
                System.out.println(net + "\t" + experiment + "\t" + iterations + " iterations");
                Set<Integer> genesUp = GeneIdParser.parseGeneIds(deDir + experiment + "_up.tsv");
                if(genesUp.size() == 0){
                    System.out.println("No significantly up-regulated Genes.");
                }
                else {
                    System.out.println("Up-regulated genes: " + genesUp.size());
                    double z_up = NetworkStatistics.calcZStatistic(network, genesUp, iterations);
                }
                Set<Integer> genesDown = GeneIdParser.parseGeneIds(deDir + experiment + "_down.tsv");
                if(genesDown.size() == 0) {
                    System.out.println("No significantly down-regulated Genes.");
                } else {
                    double z_down = NetworkStatistics.calcZStatistic(network, genesDown, iterations);
                }
            }

        }
    }
}

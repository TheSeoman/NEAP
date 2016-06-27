/**
 * Created by seoman on 5/24/16.
 */

public class Main {
    public static void main(String[] args) {
        String experiment = "GDS4395";
        Network network = NetworkParser.readBinaryNetwork("/home/sch/schmidtju/IntellijProjects/NEAP/networks/all_tissues", "/home/sch/schmidtju/IntellijProjects/NEAP/all_genes.txt");
//        ExpressionData cancer = ExpressionParser.parseSoftGz("/home/seoman/Documents/NEAP/Prostate Cancer/" + experiment + "_full.soft.gz", new String[]{"baseline"}, "Gene ID", "counts");
//        ExpressionParser.saveExpressionData(cancer, "/home/seoman/Documents/NEAP/Prostate Cancer/FilteredCountFiles/" + experiment + "_cancer");
//        ExpressionData normal = ExpressionParser.parseSoftGz("/home/seoman/Documents/NEAP/Prostate Cancer/" + experiment + "_full.soft.gz", new String[]{"control"}, "Gene ID", "counts");
//        ExpressionParser.saveExpressionData(normal, "/home/seoman/Documents/NEAP/Prostate Cancer/FilteredCountFiles/" + experiment + "_normal");
//        ExpressionData data = ExpressionParser.parseExpressionData("/home/seoman/Documents/NEAP/Prostate Cancer/FilteredCountFiles/" + experiment + "_cancer", "count");
//        ExpressionStatistics.saveExpressionCorrelations(data, "/home/seoman/Documents/NEAP/Prostate Cancer/Correlations/" + experiment + "_cancer_cor");
        ExpressionData dat = ExpressionParser.parseExpressionData("/home/seoman/Documents/NEAP/Prostate Cancer/FilteredCountFiles/GDS4114_cancer", "counts");
        ExpressionStatistics.saveExpressionCorrelations(dat, "/home/seoman/Documents/NEAP/Prostate Cancer/GDS4114_cancer_corr");
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
}

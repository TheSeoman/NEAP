import java.util.Map;

/**
 * Created by seoman on 5/24/16.
 */
public class Main {
    public static void main(String[] args) {
        ExpressionData dat = SoftParser.parseSoftGz("E:\\OneDrive\\Uni\\Bioinformatik\\Masterpraktikum NEAP\\Prostate Cancer\\GDS4395_full.soft.gz", new String[]{"baseline"}, "Gene ID","counts");
        SoftParser.saveExpressionData(dat, "GDS4395_cancer");
        if(args.length < 2) {
            printUsage();
        } else if(args[0].equals("-generateFasta")){
            try{

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
        System.out.println("motifs.jar -generateFasta <>");
        System.out.println("motifs.jar -calcCorrelationzScore <fimo output file> <output file> <p-value threshhold> <number of motifs> <ignore zeroes (true|false)>");
        System.out.println("motifs.jar -filterForGo <positive GO file> <negative GO file> <correlation zScore file>");
        System.out.println("motifs.jar -countIntoBins <correlation zScore file>");
    }
}

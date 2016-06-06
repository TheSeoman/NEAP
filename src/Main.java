import java.util.Map;

/**
 * Created by seoman on 5/24/16.
 */
public class Main {
    public static void main(String[] args) {
//        MotifList mlist = MotifParser.readMotifs("C:\\Users\\Seoman\\Documents\\NEAP\\fimo.txt", 0.001, 510);
//        MotifStatistics.calculatePCCs(mlist, "C:\\Users\\Seoman\\Documents\\NEAP\\pcc_z.tsv", false, false);
//        MotifStatistics.calculatePCCs(mlist, "/home/seoman/Documents/NEAP/pcc_z_test.tsv", false);
//        System.out.println("done");
//        changeToEntrezId("/home/sch/schmidtju/NEAP/EnsemblEntrezMapping.tsv", "/home/sch/schmidtju/NEAP/EnsemblUpstreamSeqs.fasta", "/home/sch/schmidtju/NEAP/EntrezUpstreamSeqs.fasta");
        GOParser.readPositiveNegativePairs("C:\\Users\\Seoman\\Documents\\NEAP\\GO\\go_positive_set_propagated.csv");
//        MotifParser.filterForGo(GOParser.readPositiveNegativePairs("C:\\Users\\Seoman\\Documents\\NEAP\\GO\\go_positive_set_propagated.csv"),
//                GOParser.readPositiveNegativePairs("C:\\Users\\Seoman\\Documents\\NEAP\\GO\\go_negative_set_propagated.csv"),
//                "C:\\Users\\Seoman\\Documents\\NEAP\\fimo.txt");
    }
}

/**
 * Created by seoman on 5/24/16.
 */
public class Main {
    public static void main(String[] args) {
        MotifList mlist = MotifParser.readMotifs("/home/seoman/Documents/NEAP/fimo_out/fimo.txt", 0.001, 510);
        MotifStatistics.calculatePCCs(mlist, "/home/seoman/Documents/NEAP/pcc_z.tsv", false);
        System.out.println("done");
//        changeToEntrezId("/home/sch/schmidtju/NEAP/EnsemblEntrezMapping.tsv", "/home/sch/schmidtju/NEAP/EnsemblUpstreamSeqs.fasta", "/home/sch/schmidtju/NEAP/EntrezUpstreamSeqs.fasta");
    }
}

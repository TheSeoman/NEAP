import java.util.Map;

/**
 * Created by seoman on 5/24/16.
 */
public class Main {
    public static void main(String[] args) {
        MotifList mlist = MotifParser.readMotifs("/home/seoman/Documents/NEAP/fimo_out/fimo.txt", 0.001, 510);
        MotifStatistics.calculatePCCs(mlist, "/home/seoman/Documents/NEAP/pcc_z_ignoreZero.tsv", true);
//        MotifStatistics.calculatePCCs(mlist, "/home/seoman/Documents/NEAP/pcc_z_test.tsv", false);
//        System.out.println("done");
//        changeToEntrezId("/home/sch/schmidtju/NEAP/EnsemblEntrezMapping.tsv", "/home/sch/schmidtju/NEAP/EnsemblUpstreamSeqs.fasta", "/home/sch/schmidtju/NEAP/EntrezUpstreamSeqs.fasta");
//        Map<String, String> idMap = GtfParser.parseMappingFile("/home/sch/schmidtju/IntellijProjects/NEAP/EnsemblEntrezMapping.tsv",
//                "/home/sch/schmidtju/IntellijProjects/NEAP/EntrezIds", 1, 0);
//        GtfParser.generateFasta(idMap,
//                "/home/proj/biosoft/praktikum/genomes/human/hg19/standard_chr/Homo_sapiens.GRCh37.75.gtf",
//                "/home/sch/schmidtju/IntellijProjects/NEAP/UpstreamSeqsBiomartMap.fasta");
    }
}

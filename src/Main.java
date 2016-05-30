/**
 * Created by seoman on 5/24/16.
 */
public class Main {
    public static void main(String[] args) {
//        MotifList mlist = MotifParser.readMotifs("/home/seoman/Documents/NEAP/fimo_out/fimo.txt", 0.001, 510);
//        MotifStatistics.calculatePCCs(mlist, "/home/seoman/Documents/NEAP/pcc_z.tsv", false);
//        System.out.println("done");
//        changeToEntrezId("/home/sch/schmidtju/NEAP/EnsemblEntrezMapping.tsv", "/home/sch/schmidtju/NEAP/EnsemblUpstreamSeqs.fasta", "/home/sch/schmidtju/NEAP/EntrezUpstreamSeqs.fasta");
        GtfParser.generateFasta("/home/sch/schmidtju/IntellijProjects/NEAP/gene2ensembl",
                "/home/sch/schmidtju/IntellijProjects/NEAP/EntrezIds",
                "/home/proj/biosoft/praktikum/genomes/human/hg19/standard_chr/Homo_sapiens.GRCh37.75.gtf",
                "/home/sch/schmidtju/IntellijProjects/NEAP/selfExtracted.fasta");
//        GtfParser.parseMappingFile("/home/sch/schmidtju/IntellijProjects/NEAP/gene2ensembl",
//                "/home/sch/schmidtju/IntellijProjects/NEAP/EntrezIds");
    }
}

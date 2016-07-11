import java.io.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.Integer;
import java.lang.String;
import java.lang.System;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * Created by schmidtju on 30.05.16.
 */
public class GtfParser {
    public static void generateFastaFromGff3(String idPath, String gff3Path, String fastaPath) {
        Set<String> entrezIds = new HashSet<>();
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(fastaPath));
            BufferedReader br = new BufferedReader(new FileReader(idPath));
            String line;
            while ((line = br.readLine()) != null) {
                entrezIds.add(line);
            }
            br.close();
            InputStream fileStream = new FileInputStream(gff3Path);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
            br = new BufferedReader(decoder);
            String[] split;
            String geneId = "", seq = "", chromosome = "1";
            int start, stop, genesFound = 0;
            boolean strand, region = false;
            Pattern pRegion = Pattern.compile("Dbxref=taxon:9606;chromosome=([^;]+);gbkey=Src;genome=chromosome;mol_type=genomic DNA");
            Pattern pGeneId = Pattern.compile("GeneID:(\\d+)");
            while ((line = br.readLine()) != null) {
                if (!line.startsWith("#")) {
                    split = line.split("\t");
                    if (split[2].equals("region")){
                        Matcher matcher = pRegion.matcher(split[8]);
                        if(matcher.find()){
                            chromosome = matcher.group(1);
                            region = true;
                        } else {
                            region = false;
                        }
                    } else if (region && split[2].equals("gene")) {
                        Matcher matcher = pGeneId.matcher(split[8]);
                        if (matcher.find()) {
                            geneId = matcher.group(1);
                        }
                        if (entrezIds.contains(geneId)) {
                            genesFound++;
                            start = Integer.parseInt(split[3]);
                            stop = Integer.parseInt(split[4]);
                            strand = split[6].equals("+");

                            if (strand)
                                seq = GenomeSequenceExtractor.getSequence(chromosome, start - 1000, start - 1, true);
                            else
                                seq = GenomeSequenceExtractor.getSequence(chromosome, stop + 1, stop + 1000, false);

                            out.write(">" + geneId + "\n");
                            int i = 0;
                            for (; i < seq.length() - 60; i += 60) {
                                out.write(seq.substring(i, i + 60) + "\n");
                            }
                            out.write(seq.substring(i, seq.length()) + "\n");
                        }
                    }
                }
            }
            out.close();
            System.out.println(genesFound);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void generateFastaFromEnsemblGtf(Map<String, String> idMap, String gtfPath, String fastaPath) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(fastaPath));
            String line;

            BufferedReader br = new BufferedReader(new FileReader(gtfPath));
            String[] split;
            String geneId = "", seq = "", chromosome = "1";
            int start, stop;
            boolean strand;
            Pattern pGeneId = Pattern
                    .compile("gene_id \"([^\"]+)\";");
            while ((line = br.readLine()) != null) {
                split = line.split("\t");
                if (split[2].equals("gene")) {
                    Matcher matcher = pGeneId.matcher(split[8]);
                    if (matcher.find()) {
                        geneId = matcher.group(1);
                    }
                    if (idMap.containsKey(geneId)) {
                        chromosome = split[0];
                        start = Integer.parseInt(split[3]);
                        stop = Integer.parseInt(split[4]);
                        strand = split[6].equals("+");

                        if (strand)
                            seq = GenomeSequenceExtractor.getSequence(chromosome, start - 1000, start - 1, true);
                        else
                            seq = GenomeSequenceExtractor.getSequence(chromosome, stop + 1, stop + 1000, false);

                        out.write(">" + idMap.get(geneId) + "\n");
                        int i = 0;
                        for (; i < seq.length() - 60; i += 60) {
                            out.write(seq.substring(i, i + 60) + "\n");
                        }
                        out.write(seq.substring(i, seq.length()) + "\n");

                    }
                }
            }
            out.close();
            br.close();
        } catch (FileNotFoundException fne) {
            System.out.println("File not found.");
        } catch (IOException ioe) {
            System.out.println("IOExecption.");
        }
    }
}

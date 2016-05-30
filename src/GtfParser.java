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

/**
 * Created by schmidtju on 30.05.16.
 */
public class GtfParser {
    public static Map<String, String> parseMappingFile(String mappingPath, String idPath){
        Map<String, String> ensembl2Entrez = new HashMap<>();
        Set<String> entrezIds = new HashSet<>();
        try {
//            Map<String, String> entrez2Ensembl = new HashMap<>();
            String line;
            String[] split;
            BufferedReader br = new BufferedReader(new FileReader(idPath));
            while ((line = br.readLine()) != null) {
                entrezIds.add(line);
            }
            br.close();
            br = new BufferedReader(new FileReader(mappingPath));
            while ((line = br.readLine()) != null) {
                split = line.split("\t");

                if(split.length > 2 && entrezIds.contains(split[1])){
                    ensembl2Entrez.put(split[2], split[1]);
                }
            }
        } catch (IOException e){
            e.printStackTrace();
        }
        return ensembl2Entrez;
    }

    public static void generateFasta(String mappingPath, String idPath, String gtfPath, String fastaPath) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(fastaPath));
            Map<String, String> idMap = parseMappingFile(mappingPath, idPath);
            String line;

            BufferedReader br = new BufferedReader(new FileReader(gtfPath));
            String[] split;
            String geneId = "", seq = "", chromosome = "1";
            int start, stop;
            boolean strand;
            Pattern pAttributes = Pattern
                    .compile("gene_id \"([^\"]+)\";");
            while ((line = br.readLine()) != null) {
                split = line.split("\t");
                if (split[2].equals("gene")) {
                    Matcher matcher = pAttributes.matcher(split[8]);
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

                        out.write(">" + idMap.get(geneId)+ "\n");
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

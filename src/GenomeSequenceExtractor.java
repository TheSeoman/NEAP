/**
 * Created by schmidtju on 30.05.16.
 */
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;import java.lang.String;import java.lang.StringBuilder;


public class GenomeSequenceExtractor {
    public static String getSequence(String chromosome, long start, long end, boolean strand){ //true = +, false = - strand
        String sequence = "";
        try {
            RandomAccessFile raf = new RandomAccessFile(new File("/home/proj/biosoft/praktikum/genomes/human/hg19/standard_chr/SPLITTED/" + String.valueOf(chromosome) + ".fasta"), "r");
            int s = raf.readLine().length();
            long startPos = s + start + start/60;
            long endPos = s + end + end/60 + 1;
            raf.seek(startPos);
            byte[] b = new byte[(int)((endPos - startPos))];
            raf.read(b);
            sequence = new String(b, "UTF-8").replaceAll("\n", "");
            raf.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        if(!strand)
            sequence = reverse(sequence);
        return sequence;
    }

    private static String reverse(String dna) {
        StringBuilder rev = new StringBuilder();
        for (int i = dna.length() - 1; i >= 0; i--) {
            switch (dna.charAt(i)) {
                case 'A':
                    rev.append('T');
                    break;
                case 'T':
                    rev.append('A');
                    break;
                case 'G':
                    rev.append('C');
                    break;
                case 'C':
                    rev.append('G');
                    break;
                case ' ':
                    rev.append(' ');
                    break;
            }
        }
        return rev.toString();
    }
}

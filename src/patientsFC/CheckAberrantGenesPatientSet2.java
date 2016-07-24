package patientsFC;

/**
 * Created by Stefan on 21.07.2016.
 */

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.*;

/*
 * Object network reads in a network from file "path" and uses the file "genes" to map EntrezIDs to their position.
 * Genes are given in the all_genes.txt, it is import so specifically use this file.
 * The getEdge(genA,genB) method allows direct access to the edge weight between genA and GenB.
 */


public class CheckAberrantGenesPatientSet2 {

    public static void main(String[] args) throws IOException {
        String startTime = Calendar.getInstance().getTime().toString().split("\\s+")[3];
        System.out.println("Start at "+startTime);

        String all_genes = "C:/Users/Stefan/Desktop/MaPra/all_genes.txt";
        String path = "C:/Users/Stefan/Desktop/TissuesGIANT/RAF/prostate_gland";

        File[] patientFolder = new File("C:/Users/Stefan/Desktop/BLOCKPHASE/NEAP/Prostate Cancer/patient fcs/PATIENT_SET2/").listFiles();
        File malacardsFile = new File("C:/Users/Stefan/Desktop/BLOCKPHASE/NEAP/Prostate Cancer/Malacards/all_unique_prad.txt");

        double fcThreshold = 1.0;

        CheckAberrantGenesPatientSet2 n = new CheckAberrantGenesPatientSet2(path, all_genes, malacardsFile, patientFolder, fcThreshold);

        System.out.println(fcPatientMap.size());

        String endTime = Calendar.getInstance().getTime().toString().split("\\s+")[3];
        System.out.println("End OWN at "+endTime);
    }

    private int geneSize;
    private ByteBuffer[] input;
    static HashMap<Integer, Integer> genMap = new HashMap<Integer, Integer>();
    private static HashMap<Integer, ArrayList<Integer>> fcPatientMap;
    private static HashSet<Integer> malacardsGenes = new HashSet<Integer>();

    private boolean isInsideMap;

    public CheckAberrantGenesPatientSet2(String path, String genes, File malacardsFile, File[] patientFolder, double fcThreshold) throws IOException {

        readInGenes(genes);
        geneSize = genMap.size();
        input = new ByteBuffer[geneSize];
        readRaf(path);

        createAberrantGenes(patientFolder, fcThreshold);

        readMalacardsGenes(malacardsFile);

        print(fcPatientMap);


    }

    private void print(HashMap<Integer, ArrayList<Integer>> fcMapPatients) {
        int counter = 0;
        for (Integer i : fcMapPatients.keySet()) {
            if (malacardsGenes.contains(i)) {
                int occurrences = Collections.frequency(fcMapPatients.get(i), 1);
                if (occurrences >= 80) {
                    System.out.println("Malacards Gene : " + i);
                    counter++;
                    //System.out.println(i+"("+occurrences+") : "+fcMapPatients.get(i));
                }
            }
        }
        System.out.println("Counter = "+counter);
    }

    private void readMalacardsGenes(File malacardFile) {
        String sLine = null;

        // get genes of interest
        try (BufferedReader bur = new BufferedReader(new FileReader(malacardFile))) {
            while ((sLine = bur.readLine()) != null) {
                malacardsGenes.add(Integer.parseInt(sLine.split("\\t")[1]));
            }

            bur.close();
        } catch (FileNotFoundException eFNF) {
            eFNF.printStackTrace();
        } catch (IOException eIO) {
            eIO.printStackTrace();
        }
    }

    /**
     * get the aberrant genes in every patient
     *
     * @param patientFolder folder to patient sets in fold change
     */
    private void createAberrantGenes(File[] patientFolder, double fcThreshold) {
        fcPatientMap = new HashMap<Integer, ArrayList<Integer>>();

        String sLine = null;

        for (File file : patientFolder) {
            try (BufferedReader bur = new BufferedReader(new FileReader(file))) {
                while ((sLine = bur.readLine()) != null) {
                    int id = Integer.parseInt(sLine.split("\\t")[0]);

                    // check if gene of interest
                    //if (malacardsGenes.contains(id)) {
                    double fcValue = Double.parseDouble(sLine.split("\\t")[1]);

                    ArrayList<Integer> l = fcPatientMap.get(id);

                    // 1 if aberrant
                    // 0 if not
                    int fcAberrant = (fcValue >= fcThreshold || fcValue <= ((-1) * fcThreshold)) ? 1 : 0;

                    if (l == null) {
                        ArrayList<Integer> cur = new ArrayList<>();
                        cur.add(fcAberrant);
                        fcPatientMap.put(id, cur);
                    } else {
                        l.add(fcAberrant);
                        fcPatientMap.put(id, l);
                    }
                    // }
                }

                bur.close();
            } catch (FileNotFoundException eFNF) {
                eFNF.printStackTrace();
            } catch (IOException eIO) {
                eIO.printStackTrace();
            }
        }
    }


	/*
     * readInGenes: read in list of genes(nodes) present in CheckAberrantGenesPatientSet2 and generate hashmap.
	 */

    private static void readInGenes(String genes) {

        System.out.println("Read in genes in network..");
        InputStream is = null;
        try {
            is = new FileInputStream(genes);
            BufferedReader br = new BufferedReader(new InputStreamReader(is));

            String line = br.readLine();
            int count = 0;

            while (line != null) {

                line = line.trim();
                genMap.put(Integer.parseInt(line), count);
                count++;


                line = br.readLine();
            }

            br.close();
        } catch (IOException e) {
            System.err.print("File can not be read.");
            e.printStackTrace();
        }

    }




	/*
     * readRaf: read in as many bytes as necessary and save them in ByteBuffer array
	 * array length: number of nodes in network
	 * buffer lengths: incrementing=> half matrix
	 */

    private void readRaf(String file) throws IOException {
        RandomAccessFile raf = new RandomAccessFile(file, "r");
        FileChannel channel = raf.getChannel();

        try {
            ByteBuffer bb;
            for (int i = 0; i < geneSize; i++) {
                bb = ByteBuffer.allocate((i+1) * 8);
                channel.read(bb);
                bb.flip();
                input[i] = bb;
                bb.clear();

            }
        } finally {

            channel.close();
            raf.close();
        }

    }

	/*
     * getEdge: inputs two entrezIDs and returns the corresponding edge weight
	 */

    public double getEdge(Integer a, Integer b) {
        this.setInsideMap(false);
        a = genMap.get(a);
        b = genMap.get(b);
        if (a != null && b != null) {
            this.setInsideMap(true);
            if (a <= b) {
                return (input[b].getDouble((a) * 8));

            } else {

                return (input[a].getDouble((b) * 8));

            }

        }
        return 0.0;
    }

    public boolean isInsideMap() {
        return isInsideMap;
    }

    public void setInsideMap(boolean isInsideMap) {
        this.isInsideMap = isInsideMap;
    }
}

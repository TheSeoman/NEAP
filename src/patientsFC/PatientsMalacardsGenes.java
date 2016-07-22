package patientsFC;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by Stefan on 22.07.2016.
 */
public class PatientsMalacardsGenes {

    public static void main(String[] args) throws IOException {
        String path = "C:/Users/Stefan/Desktop/BLOCKPHASE/NEAP/Prostate Cancer/";

        File malacardsFile = new File(path+"Malacards/all_unique_prad.txt");

        File[] patientFolder = new File(path+"patient fcs/PATIENT_SET1/").listFiles();

        double threshold = 2.0;

        new PatientsMalacardsGenes(malacardsFile, patientFolder, threshold);
    }

    private HashSet<Integer> malacardsGenes;
    private HashMap<String, HashMap<Integer, Double>> patientMalacardsFC;

    public PatientsMalacardsGenes(File malacardsFile, File[] patientFolder, double threshold) throws IOException {
        BufferedReader bur = openReader(malacardsFile);

        malacardsGenes = new HashSet<Integer>();

        String sLine = null;

        while ((sLine = bur.readLine()) != null) {
            malacardsGenes.add(Integer.parseInt(sLine.split("\\t")[1]));
        }

        bur.close();

        patientMalacardsFC = new HashMap<String, HashMap<Integer, Double>>();

        for (File patientFile : patientFolder) {
            BufferedReader burP = openReader(patientFile);

            String curPatient = patientFile.getName().split("\\.")[0];

            while ((sLine = burP.readLine()) != null) {
                int gene = Integer.parseInt(sLine.split("\\t")[0]);

                // check if malacards gene
                if (malacardsGenes.contains(gene)) {

                    double fcValue = Double.parseDouble(sLine.split("\\t")[1]);

                    // check if aberrant gene
                    if (fcValue >= threshold || fcValue <= (threshold * (-1))) {

                        HashMap<Integer, Double> cur = patientMalacardsFC.get(curPatient);

                        if (cur == null) {
                            HashMap<Integer, Double> m = new HashMap<Integer, Double>();
                            m.put(gene, fcValue);
                            patientMalacardsFC.put(curPatient, m);
                        } else {
                            cur.put(gene, fcValue);
                            patientMalacardsFC.put(curPatient, cur);
                        }
                    }
                }
            }

            burP.close();
        }

        for (String s : patientMalacardsFC.keySet()) {
            HashMap<Integer, Double> m = patientMalacardsFC.get(s);
            double x = 0.0;
            for (Double d : m.values()) {
                x += Math.abs(d);
            }
            System.out.println(s+" : "+patientMalacardsFC.get(s).size()+"\t"+(x/patientMalacardsFC.get(s).size()));
        }
    }

    private BufferedReader openReader(File file) throws IOException {
        try {
            return new BufferedReader(new FileReader(file));
        } catch (FileNotFoundException eFNF) {
            eFNF.printStackTrace();
        }
        return null;
    }

}

package patientsFC;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by Stefan on 22.07.2016.
 */
public class AberrantGenes {

    public static void main(String[] args) throws IOException {
        String path = "C:/Users/Stefan/Desktop/BLOCKPHASE/";
        File[] patientFolder = new File(path+"NEAP/Prostate Cancer/patient fcs/PATIENT_SET1/").listFiles();
        File pradPatientsSet1 = new File(path+"PRAD_patients_set_1.txt");
        File[] tcgaPatients = new File(path+"NEAP/Prostate Cancer/FoldChange/prostate").listFiles();
        File malacardsFile = new File(path+"NEAP/Prostate Cancer/Malacards/all_unique_prad.txt");


        double threshold = 1.0;

        new AberrantGenes(pradPatientsSet1, patientFolder, tcgaPatients, malacardsFile, threshold);
    }

    /**
     *
     * @param pradPatientsFile
     * @param patientFolder
     * @param tcgaFolder
     * @param malacardsFile
     * @param threshold
     * @throws IOException
     */
    public AberrantGenes(File pradPatientsFile, File[] patientFolder, File[] tcgaFolder, File malacardsFile, double threshold) throws IOException {
        receivePradPatients(pradPatientsFile);
        setNumberOfPatients(getNumberOfUniquePatients(tcgaFolder, patientFolder));
        parseTCGAPatientsAndSet(patientFolder, tcgaFolder, threshold);

        HashSet<Integer> malacardsGenes = new HashSet<Integer>();

        BufferedReader bur = openReader(malacardsFile);
        String sLine = null;

        while ((sLine = bur.readLine()) != null) {
            int id = Integer.parseInt(sLine.split("\\t")[1]);
            malacardsGenes.add(id);
        }

        int c = 0;
        for (Integer i : malacardsGenes) {
            if (aberrantGeneMap.containsKey(i)) {
                c++;
            } else {
                System.out.println(i);
            }
        }
        System.out.println(c + "!!");

        bur.close();
    }

    private HashMap<Integer, int[]> aberrantGeneMap;
    private HashMap<String, Integer> patientNumberMap;
    private HashSet<String> pradPatientSet;
    private int NUMBER_OF_PATIENTS;



    /**
     *
     * @param pradPatientsFile
     * @return Set with PRAD patients in pradPatientsFile
     * @throws IOException
     */
    private HashSet<String> receivePradPatients(File pradPatientsFile) throws IOException {
        // get prad patients of patient set1
        pradPatientSet = new HashSet<String>();
        BufferedReader br = openReader(pradPatientsFile);
        String sLine = null;

        while ((sLine = br.readLine()) != null) {
            pradPatientSet.add(sLine.trim());
        }

        br.close();

        return pradPatientSet;
    }

    /**
     * @param patientFolder
     * @param tcgaFolder
     * @param threshold
     * @return Map with key: gene; value: aberrant genes array
     * @throws IOException
     */
    private HashMap<Integer, int[]> parseTCGAPatientsAndSet(File[] patientFolder, File[] tcgaFolder, double threshold) throws IOException {
        aberrantGeneMap = new HashMap<Integer, int[]>();

        // buffered reader line
        String sLine = null;

        // tcga
        for (File file : tcgaFolder) {
            BufferedReader bur = openReader(file);

            String patient = file.getName().split("\\.")[0];

            while ((sLine = bur.readLine()) != null) {
                int id = Integer.parseInt(sLine.split("\\t")[0]);
                double fcValue = Double.parseDouble(sLine.split("\\t")[1]);

                this.fillMap(aberrantGeneMap, patient, id, fcValue, threshold);
            }

            bur.close();
        }

        // patient set1
        for (File file : patientFolder) {
            String patient = file.getName().split("\\.")[0];

            // take only prad patients into account
            if (pradPatientSet.contains(patient)) {

                BufferedReader bur = openReader(file);

                while ((sLine = bur.readLine()) != null) {
                    int id = Integer.parseInt(sLine.split("\\t")[0]);
                    double fcValue = Double.parseDouble(sLine.split("\\t")[1]);

                    this.fillMap(aberrantGeneMap, patient, id, fcValue, threshold);
                }

                bur.close();
            }
        }

        System.out.println(aberrantGeneMap.size());
        return aberrantGeneMap;
    }

    /**
     * @param folderTCGA
     * @param folderSet
     * @return number of all patients in the folders
     * @throws IOException
     */
    private int getNumberOfUniquePatients(File[] folderTCGA, File[] folderSet) throws IOException {
        int count = 0;
        patientNumberMap = new HashMap<String, Integer>();
        for (File f : folderTCGA) {
            String patient = f.getName().split("\\.")[0];
            patientNumberMap.put(patient, count);
            count++;
        }
        for (File f : folderSet) {
            String patient = f.getName().split("\\.")[0];
            if (pradPatientSet.contains(patient)) {
                patientNumberMap.put(patient, count);
                count++;
            }
        }
        return count;
    }

    /**
     * @param map       aberrantGeneMap
     * @param patient   current patient
     * @param id        gene id
     * @param fcValue   fold change for current gene
     * @param threshold fold change threshold
     * @return new map
     */
    private HashMap<Integer, int[]> fillMap(HashMap<Integer, int[]> map, String patient, int id, double fcValue, double threshold) {
        int[] curGene = map.get(id);

        int aberrant = (fcValue >= threshold || fcValue <= (threshold * (-1))) ? 1 : 0;

        if (curGene == null) {
            int[] cur = new int[getNumberOfPatients()];
            cur[0] = aberrant;
            map.put(id, cur);
        } else {
            int patientPosition = patientNumberMap.get(patient);
            curGene[patientPosition] = aberrant;
            map.put(id, curGene);
        }

        return map;
    }

    /**
     * @param file to be opened
     * @return buffered reader for foile
     * @throws IOException
     */
    private BufferedReader openReader(File file) throws IOException {
        try {
            return new BufferedReader(new FileReader(file));
        } catch (FileNotFoundException eFNF) {
            eFNF.printStackTrace();
        }
        return null;
    }

    public int getNumberOfPatients() {
        return NUMBER_OF_PATIENTS;
    }

    public void setNumberOfPatients(int NUMBER_OF_PATIENTS) {
        this.NUMBER_OF_PATIENTS = NUMBER_OF_PATIENTS;
    }

}

package patientsFC;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by Stefan on 21.07.2016.
 */
public class SubnetworkPatients {

    public static void main(String[] args) throws IOException {
        File genesOfInterest = new File("C:/Users/Stefan/Desktop/BLOCKPHASE/NEAP/Prostate Cancer/Malacards/all_unique_prad.txt");
        File[] patientFolder = new File("C:/Users/Stefan/Desktop/PatientSet2/").listFiles();

        new SubnetworkPatients(patientFolder, genesOfInterest, 1.0);
    }

    /**
     * @param patientFolder   patient files with fold changes
     * @param genesOfInterest malacards genes of prostate cancer
     * @throws IOException
     */
    public SubnetworkPatients(File[] patientFolder, File genesOfInterest, double fcThreshold) throws IOException {
        HashSet<Integer> malacardsGenes = new HashSet<Integer>();

        String sLine = null;

        // get genes of interest
        try (BufferedReader bur = new BufferedReader(new FileReader(genesOfInterest))) {
            while ((sLine = bur.readLine()) != null) {
                malacardsGenes.add(Integer.parseInt(sLine.split("\\t")[1]));
            }

            bur.close();
        } catch (FileNotFoundException eFNF) {
            eFNF.printStackTrace();
        } catch (IOException eIO) {
            eIO.printStackTrace();
        }

        HashMap<Integer, ArrayList<Integer>> fcPatientMap = new HashMap<Integer, ArrayList<Integer>>();
        HashMap<Integer, String> patients = new HashMap<Integer, String>();
        int posPatients = 0;

        for (File file : patientFolder) {
            patients.put(posPatients++, file.getName().split("\\.")[0]);
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
//        print(fcPatientMap, malacardsGenes, 0.8);
//        System.out.println(patients.get(0));
        checkPatient(fcPatientMap, patients, 0.8);
    }

    private void checkPatient(HashMap<Integer, ArrayList<Integer>> map, HashMap<Integer, String> patients, double threshold) throws IOException {
        HashMap<String, HashMap<Integer, Integer>> patientProbabilityPRAD_Map = new HashMap<String, HashMap<Integer, Integer>>();
        for (Integer i : map.keySet()) {
            int aberrantGenes = Collections.frequency(map.get(i), 1);

            if (aberrantGenes >= (int) (threshold * 100)) {
                for (int k = 0; k < map.get(i).size(); k++) {
                    String curPatient = patients.get(k);
                    int curGeneAberrantOrNot = map.get(i).get(k);

//                    System.out.println(k+"  "+curPatient+" "+i+" "+curGeneAberrantOrNot);

                    HashMap<Integer, Integer> temp = patientProbabilityPRAD_Map.get(curPatient);

                    if (temp == null) {
                        HashMap<Integer, Integer> curTemp = new HashMap<Integer, Integer>();
                        curTemp.put(i, curGeneAberrantOrNot);
                        patientProbabilityPRAD_Map.put(curPatient, curTemp);
                    } else {
                        temp.put(i, curGeneAberrantOrNot);
                        patientProbabilityPRAD_Map.put(curPatient, temp);
                    }
                }
            }
        }
        for (String s : patientProbabilityPRAD_Map.keySet()) {
            HashMap<Integer, Integer> m = patientProbabilityPRAD_Map.get(s);
            BufferedWriter wr = new BufferedWriter(new FileWriter(new File("C:/Users/Stefan/Desktop/PRAD_PATIENT_SET2/"+s+".txt")));
            System.out.print(s);
            for (Integer i : m.keySet()) {
                System.out.print(" "+i+" "+m.get(i));
                wr.write(i+"\t"+m.get(i));
                wr.newLine();
            }
            System.out.println();
            wr.close();
        }

    }

    private void print(HashMap<Integer, ArrayList<Integer>> fcMapPatients, HashSet<Integer> malacards, double threshold) {
        int counter = 0;
        int c = 0;
        for (Integer i : fcMapPatients.keySet()) {
            int occurrences = Collections.frequency(fcMapPatients.get(i), 1);
            if (occurrences >= (int) (threshold * 100)) {
                System.out.println(c+++" : "+i+"("+occurrences+") : "+fcMapPatients.get(i).size());

                if (malacards.contains(i)) {
                    counter++;
                }
            }
        }
        System.out.println("Prostate Cancer Genes = "+counter);
    }

}

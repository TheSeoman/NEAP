import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

/**
 * Created by Seoman on 22.06.2016.
 */
public class ExpressionParser {
    public static ExpressionData parseSoftGz(String path, String[] conditions, String idColName, String type) {
        HashMap<String, Integer> idMap = new HashMap<String, Integer>();
        int idc = 0;
        List<double[]> values = new ArrayList<>();
        List<Integer> occurrences = new ArrayList<>();
        String line;
        Map<String, Set<String>> samples = new HashMap<String, Set<String>>();
        String subsetId = "";
        String subsetDescription = "";
        String[] subsetSampleId = new String[0];
        String[] tableHead;
        int idCol = 0;
        int[] sampleCols = null;
        boolean head = true;
        for (String condition : conditions) {
            samples.put(condition, new HashSet<String>());
        }
        Set<String> finalSamples = null;

        try {
            InputStream fileStream = new FileInputStream(path);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
            BufferedReader br = new BufferedReader(decoder);
            while ((line = br.readLine()) != null) {
                if (head) {
                    if (line.startsWith("^SUBSET")) {
                        for (String condition : conditions) {
                            if (subsetDescription.indexOf(condition) != -1) {
                                for (String sample : subsetSampleId) {
                                    samples.get(condition).add(sample);
                                }
                            }
                        }
                        subsetId = "";
                        subsetDescription = "";
                        subsetSampleId = new String[0];
                    } else if (line.startsWith("!subset_dataset_id")) {
                        subsetId = line.split(" = ")[1];
                    } else if (line.startsWith("!subset_description")) {
                        subsetDescription = line.split(" = ")[1];
                    } else if (line.startsWith("!subset_sample_id")) {
                        subsetSampleId = line.split(" = ")[1].split(",");
                    } else if (line.startsWith("!dataset_table_begin")) {
                        for (String condition : conditions) {
                            if (subsetDescription.indexOf(condition) != -1) {
                                for (String sample : subsetSampleId) {
                                    samples.get(condition).add(sample);
                                }
                            }
                        }
                        head = false;
                        line = br.readLine();
                        tableHead = line.split("\t");
                        if (samples.size() == 0) {
                            System.out.println("No samples found for keyword: " + conditions);
                            break;
                        } else {
                            finalSamples = samples.get(conditions[0]);
                            for (Set<String> sample : samples.values()) {
                                finalSamples.retainAll(sample);
                            }
                            sampleCols = new int[finalSamples.size()];
                            int c = 0;
                            for (int i = 0; i < tableHead.length; i++) {
                                if (finalSamples.contains(tableHead[i])) {
                                    sampleCols[c] = i;
                                    c++;
                                } else if (tableHead[i].equals(idColName)) {
                                    idCol = i;
                                }
                            }
                        }
                    }
                } else {
                    if (line.startsWith("!dataset_table_end")) {
                        break;
                    } else {
                        String[] split = line.split("\t");
                        if (idCol < split.length && !split[idCol].equals("")) {
                            double[] value = new double[finalSamples.size()];
                            for (int i = 0; i < sampleCols.length; i++) {
                                if (!split[sampleCols[i]].equals("null")) {
                                    value[i] = Double.parseDouble(split[sampleCols[i]]);
                                } else {
                                    value[i] = Double.NaN;
                                }
                            }
                            for (String geneId : split[idCol].split("///")) {
                                if (idMap.containsKey(geneId)) {
                                    int idInt = idMap.get(geneId);
                                    int prevOccurrences = occurrences.get(idInt);
                                    occurrences.set(idInt, prevOccurrences + 1);
                                    double[] oldValues = values.get(idInt);
                                    for (int i = 0; i < value.length; i++) {
                                        oldValues[i] = ((oldValues[i] * prevOccurrences) + value[i]) / (prevOccurrences + 1);
                                    }
                                } else {
                                    occurrences.add(1);
                                    idMap.put(geneId, idc);
                                    idc++;
                                    values.add(value);
                                }
                            }
                        }
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        double[][] valuesArr = new double[values.size()][samples.size()];
        int i = 0;
        for (double[] value : values) {
            valuesArr[i] = value;
            i++;
        }
        return new ExpressionData(type, valuesArr, idMap, finalSamples.toArray(new String[samples.size()]));
    }

    public static void saveExpressionData(ExpressionData ex, String path) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(path));
            out.write("id");
            for (String sample : ex.getSamples()) {
                out.write("\t" + sample);
            }
            out.write("\n");
            for (String id : ex.getIdMap().keySet()) {
                out.write(id);
                for (double value : ex.getValues()[ex.getIdMap().get(id)]) {
                    out.write("\t" + value);
                }
                out.write("\n");
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static ExpressionData parseExpressionData(String path, String type, boolean head) {
        List<double[]> values = new ArrayList<>();
        Map<String, Integer> idMap = new HashMap<>();
        String[] samples = null;
        int c = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            String[] split;
            if (head) {
                line = br.readLine();
                split = line.split("\t");
                samples = Arrays.copyOfRange(split, 1, split.length);
            }
            while ((line = br.readLine()) != null) {
                split = line.split("\t");
                String id = split[0];
                double[] value = new double[split.length - 1];
                for (int i = 1; i < split.length; i++) {
                    value[i - 1] = Double.parseDouble(split[i]);
                }
                values.add(value);
                idMap.put(id, c);
                c++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        double[][] valuesArr = new double[values.size()][samples.length];
        int i = 0;
        for (double[] value : values) {
            valuesArr[i] = value;
            i++;
        }
        if (samples == null || samples.equals("expression")) {
            samples = new String[valuesArr[0].length];
            for (i = 0; i < samples.length; i++) {
                samples[i] = path.substring(path.lastIndexOf("/") + 1);
            }
        }
        return new ExpressionData(type, valuesArr, idMap, samples);
    }

    public static ExpressionData parseTCGAExpressionData(String path, String type, Map<String, String> ensembleMapping) {
        List<double[]> values = new ArrayList<>();
        Map<String, Integer> idMap = new HashMap<>();
        String[] samples = new String[]{path.substring(path.lastIndexOf("/") + 1, path.indexOf("."))};
        int c = 0;
        try {
            InputStream fileStream = new FileInputStream(path);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
            BufferedReader br = new BufferedReader(decoder);
            String line;
            String[] split;
            while ((line = br.readLine()) != null) {
                split = line.split("\t");
                if (split[0].contains(".")) {
                    split[0] = split[0].substring(0, split[0].indexOf("."));
                }
                if (!ensembleMapping.containsKey(split[0])) {
                    continue;
                }
                String id = ensembleMapping.get(split[0]);
                double value = Double.parseDouble(split[1]);
                if (!idMap.containsKey(id)) {
                    idMap.put(id, c);
                    values.add(new double[]{value});
                    c++;
                } else {
                    double prevValue = values.get(idMap.get(id))[0];
                    if (value > prevValue) {
                        values.get(idMap.get(id))[0] = value;
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        double[][] valuesArr = new double[values.size()][samples.length];
        int i = 0;
        for (double[] value : values) {
            valuesArr[i] = value;
            i++;
        }
        return new ExpressionData(type, valuesArr, idMap, samples);
    }

    public static Map<String, String> parseTCGAJSON(String path) {
        Map<String, String> caseToFile = new HashMap<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line, fileName = "", caseId = "";
            while ((line = br.readLine()) != null) {
                if (line.contains("file_name")) {
                    fileName = line.split("\"")[3];
                } else if (line.contains("case_id")) {
                    caseId = line.split("\"")[3];
                    caseToFile.put(caseId, fileName);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return caseToFile;
    }

    public static void saveAverageFoldChange(ExpressionData data, int[] cols1, int[] cols2, String path) {
        double[] fcs = ExpressionStatistics.calcAverageFoldChange(data, cols1, cols2);
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(path));
            for (String id : data.getIdMap().keySet()) {
                double fc = fcs[data.getIdMap().get(id)];
                if (fc != -1) {
                    out.write(id + "\t" + Math.log(fc) + "\n");
                }
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void savePatientsFoldChanges(String patientsDir) {
        File[] files = new File(patientsDir).listFiles();
        Set<String> patientIds = new HashSet<>();
        for (File file : files) {
            if (file.isFile()) {
                patientIds.add(file.getName().substring(0, file.getName().indexOf(".")));
            }
        }
        for(String id : patientIds){
            ExpressionData healtyCounts = parseExpressionData(patientsDir + id + ".normal", "counts", true);
            ExpressionData tumorCounts = parseExpressionData(patientsDir + id + ".tumor", "counts", true);
            ExpressionData combined = mergeExpressionData(healtyCounts, tumorCounts, 0);

            saveAverageFoldChange(combined, new int[]{0}, new int[]{1}, "/home/sch/schmidtju/IntellijProjects/NEAP/Prostate Cancer/patient fcs/" + id + ".fc.tsv");


        }
    }


    public static void mergeTCGACountFiles(String healtyJSON, String tumorJSON, String TCGAdir, String ensembl2entrezMapping) {
        Map<String, String> ensembl2entrez = GeneIdParser.parseMappingFile(ensembl2entrezMapping, 0, 1);
        Map<String, String> healthy = parseTCGAJSON(healtyJSON);
        Map<String, String> tumor = parseTCGAJSON(tumorJSON);
        ExpressionData healthyTotal = null;
        ExpressionData tumorTotal = null;
        for (String case1Id : healthy.keySet()) {
            for (String case2Id : parseTCGAJSON(tumorJSON).keySet()) {
                if (case1Id.equals(case2Id)) {
                    String healtyFile = TCGAdir + healthy.get(case1Id);
                    String tumorFile = TCGAdir + tumor.get(case2Id);
                    ExpressionData healthyData = parseTCGAExpressionData(healtyFile, "counts", ensembl2entrez);
                    ExpressionData tumorData = parseTCGAExpressionData(tumorFile, "counts", ensembl2entrez);
                    ExpressionData combined = mergeExpressionData(healthyData, tumorData, 0);
                    if(healthyTotal == null){
                        healthyTotal = healthyData;
                        tumorTotal = tumorData;
                    } else {
                        healthyTotal = mergeExpressionData(healthyTotal, healthyData, 0);
                        tumorTotal = mergeExpressionData(tumorTotal, tumorData, 0);
                    }

                    ExpressionData total = mergeExpressionData(healthyTotal, tumorTotal, 0);

                    saveExpressionData(total, "/home/sch/schmidtju/IntellijProjects/NEAP/Prostate Cancer/TCGAcases/total.count.tsv");

//                    saveExpressionData(combined, "/home/seoman/Documents/NEAP/Prostate Cancer/TCGAcases/" + case1Id + ".count.tsv");
//                    saveAverageFoldChange(combined, new int[]{0}, new int[]{1}, "/home/seoman/Documents/NEAP/Prostate Cancer/FoldChange/" + case1Id + ".fc.tsv");
                }
            }
        }
    }

    public static ExpressionData mergeExpressionData(ExpressionData data1, ExpressionData data2, double defaultValue) {
        ExpressionData data = null;
        if (!data1.getType().equals(data2.getType())) {
            System.out.println("Can't merge Expression data of different types");
        } else {
            Map<String, Integer> idMap = new HashMap<>();
            int data1len = data1.getSamples().length;
            int data2len = data2.getSamples().length;
            String[] samples = new String[data1len + data2len];
            for (int i = 0; i < data1len; i++) {
                samples[i] = data1.getSamples()[i];
            }
            for (int i = 0; i < data2len; i++) {
                samples[i + data1len] = data2.getSamples()[i];
            }
            int c = 0;
            for (String id1 : data1.getIdMap().keySet()) {
                idMap.put(id1, c);
                c++;
            }
            for (String id2 : data2.getIdMap().keySet()) {
                if (!idMap.containsKey(id2)) {
                    idMap.put(id2, c);
                    c++;
                }
            }
            double[][] values = new double[c + 1][data1len + data2len];
            c = 0;
            for (String id : idMap.keySet()) {
                double[] value = new double[data1len + data2len];
                if (data1.getIdMap().containsKey(id)) {
                    double[] data1value = data1.getValues()[data1.getIdMap().get(id)];
                    for (int i = 0; i < data1len; i++) {
                        values[c][i] = data1value[i];
                    }
                } else {
                    for (int i = 0; i < data1len; i++) {
                        values[c][i] = defaultValue;
                    }
                }
                if (data2.getIdMap().containsKey(id)) {
                    double[] data2value = data2.getValues()[data1.getIdMap().get(id)];
                    for (int i = 0; i < data2len; i++) {
                        values[c][i + data1len] = data2value[i];
                    }
                } else {
                    for (int i = 0; i < data1len; i++) {
                        values[c][i + data1len] = defaultValue;
                    }
                }
                c++;
            }
            data = new ExpressionData("count", values, idMap, samples);
        }
        return data;
    }

}

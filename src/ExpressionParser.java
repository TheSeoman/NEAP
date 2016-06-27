import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.Inflater;

/**
 * Created by Seoman on 22.06.2016.
 */
public class ExpressionParser {
    public static ExpressionData parseSoftGz(String path, String[] conditions, String idColName, String type) {
        HashMap<String, Integer> idMap = new HashMap<String, Integer>();
        int idc = 0;
        List<double[]> values = new ArrayList<>();
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
                                    double[] oldValues = values.get(idMap.get(geneId));
                                    for (int i = 0; i < value.length; i++) {
                                        oldValues[i] += value[i];
                                    }
                                } else {
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
            int length = ex.getValues()[0].length;
            for (String id : ex.getIdMap().keySet()) {
                out.write(id);
                if(id.equals("6548"))
                    System.out.println("");
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

    public static ExpressionData parseExpressionData(String path, String type) {
        List<double[]> values = new ArrayList<>();
        Map<String, Integer> idMap = new HashMap<>();
        String[] samples = null;
        int c = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line = br.readLine();
            String[] split = line.split("\t");
            samples = Arrays.copyOfRange(split, 1, split.length);
            while ((line = br.readLine()) != null) {
                split = line.split("\t");
                String id = split[0];
                double[] value = new double[split.length - 1];
                for(int i = 1; i < split.length; i++){
                    value[i-1] = Double.parseDouble(split[i]);
                }
                values.add(value);
                idMap.put(id, c);
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

}

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Created by seoman on 6/27/16.
 */
public class ExpressionStatistics {

    public static double[][] calcExpressionCorrelations(ExpressionData data){
        double[][] values = data.getValues();
        double[][] pccs = new double[values.length][values.length];
        for(int i = 0; i < values.length; i++){
            pccs[i][i] = 1;
            for(int j = i+1; j < values.length; j++){
                pccs[i][j] = Statistics.calculatePCC(values[i], values[j]);
                pccs[j][i] = pccs[i][j];
            }
            if((((double)i * 100 / values.length)) % 1 == 0)
                System.out.println("Progress calculationg pccs: " + ((i * 100 / values.length)) + "%");
        }
        return pccs;
    }

    public static void saveExpressionCorrelations(ExpressionData ex, String path){
        double[][] pccs = calcExpressionCorrelations(ex);
        System.out.println("Save correlations.");
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(path));

            for (String id1 : ex.getIdMap().keySet()) {
                for(String id2 : ex.getIdMap().keySet()){
                    if(id1.compareTo(id2) > 0){
                        out.write(id1 + "\t" + id2 + pccs[ex.getIdMap().get(id1)][ex.getIdMap().get(id2)]);
                        out.write("\n");
                    }
                }
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}

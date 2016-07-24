package patientsFC;

import java.io.IOException;

/**
 * Created by Stefan on 24.07.2016.
 */
public class McSubnetworkApproach2 {

    public static void main(String[] args) throws IOException {

        String allGenes = "C:/Users/Stefan/Desktop/MaPra/all_genes.txt";
        String pathGIANT = "C:/Users/Stefan/Desktop/TissuesGIANT/RAF/prostate_gland";

        SubNetwork n = new SubNetwork(pathGIANT, allGenes, false);
    }



}

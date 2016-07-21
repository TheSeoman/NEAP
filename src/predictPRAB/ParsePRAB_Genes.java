package predictPRAB;

import java.io.*;
import java.util.ArrayList;

public class ParsePRAB_Genes {
	
	/**
	 * store PRAB related genes manually filtered from malacards
	 * 
	 * @param file
	 *            tab separated malacards gene file with
	 *            <p>
	 *            UniprotID
	 *            </p>
	 *            <p>
	 *            Gene ID
	 *            </p>
	 * @return Integer List with Gene IDs related to PRAB
	 */
	public static ArrayList<Integer> getProstateCancerGenes(File file) {
		ArrayList<Integer> geneList = new ArrayList<Integer>();

		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String sLine = null;

			while ((sLine = br.readLine()) != null) {
				geneList.add(Integer.parseInt(sLine.split("\\t")[1]));
			}

			br.close();
		} catch (FileNotFoundException eFNF) {
			eFNF.printStackTrace();
		} catch (IOException eIO) {
			eIO.printStackTrace();
		}

		return geneList;
	}

}

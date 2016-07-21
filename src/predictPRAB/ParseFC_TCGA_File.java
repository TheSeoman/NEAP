package predictPRAB;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

public class ParseFC_TCGA_File {

	/**
	 * read in tcga fold changes and store in double array
	 * 
	 * @param idList
	 *            Gene IDs from PRAB malacards genes
	 * @param file
	 *            TCGA fold change file
	 * @return double array for pearson correlation
	 */
	public static HashMap<Integer, Double> readTCGA(ArrayList<Integer> idList, File file) {
		HashMap<Integer, Double> tcgaMap = new HashMap<Integer, Double>();

		try (BufferedReader bur = new BufferedReader(new FileReader(file))) {
			String sLine = bur.readLine(); // header

			while ((sLine = bur.readLine()) != null) {
				String[] entries = sLine.split("\t");

				int id = Integer.parseInt(entries[0]);
				double logFC = Double.parseDouble(entries[1]);

				if (idList.contains(id)) {
					tcgaMap.put(id, logFC);
				}
			}

			bur.close();
		} catch (FileNotFoundException eFNF) {
			eFNF.printStackTrace();
		} catch (IOException eIO) {
			eIO.printStackTrace();
		}

		return tcgaMap;
	}

}

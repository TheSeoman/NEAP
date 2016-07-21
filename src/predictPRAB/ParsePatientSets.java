package predictPRAB;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

public class ParsePatientSets {

	/**
	 * calculate pcc for each patient set with fold changes (own)
	 * 
	 * @param folder
	 *            Patient Sets folder
	 * @param idList
	 *            Gene IDs from PRAB malacards genes
	 * @return
	 */
	public static HashMap<Integer, Double> getPCCPatientSets(File folderPatient, ArrayList<Integer> idList) {
		HashMap<Integer, Double> patienMap = new HashMap<Integer, Double>();
		try (BufferedReader bur = new BufferedReader(new FileReader(folderPatient))) {

			String sLine = null;

			while ((sLine = bur.readLine()) != null) {
				String[] entries = sLine.split("\t");

				int id = Integer.parseInt(entries[0]);
				double logFC = Double.parseDouble(entries[1]);

				if (idList.contains(id)) {
					patienMap.put(id, logFC);
				}
			}

			bur.close();
		} catch (FileNotFoundException eFNF) {
			eFNF.printStackTrace();
		} catch (IOException eIO) {
			eIO.printStackTrace();
		}
		return patienMap;
	}
}

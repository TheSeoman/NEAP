package predictPRAB;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashMap;

public class MainAppl {

	public static void main(String[] args) throws IOException {
		if (args.length != 18) {
			System.out.println(MainAppl.SYNTAX_HELP);
			System.out.println();
			System.out.println(EXAMPLE_INVOKE);
			System.out.println();
			System.out.println("Aborting program...");
			return;
		}

		String startTime = Calendar.getInstance().getTime().toString().split("\\s+")[3];
		System.out.println("Start at " + startTime);

		String tissue = MainAppl.getTissue(args);
		String raf = MainAppl.getRAF_Path(args);
		String rafFile = raf + tissue;
		String networkType = MainAppl.getNetworkType(args);

		String genesFile = MainAppl.getAllGenesFile(args);

		GIANT_OW_Network n = new GIANT_OW_Network(rafFile, genesFile);

		System.out.println(rafFile + "  read in successfully...");

		File malacardsFile = new File(MainAppl.getMalacards(args));
		ArrayList<Integer> prabList = ParsePRAB_Genes.getProstateCancerGenes(malacardsFile);

		double threshold = 0.75;
		try {
			threshold = Double.parseDouble(MainAppl.getThreshold(args));

			if (threshold < 0.0 || threshold > 1.0) {
				System.out.println("Threshold must be a floating number between 0 and 1");
				System.out.println("Aborting program...");
				return;
			}
		} catch (NumberFormatException eNF) {
			System.out.println("Threshold must be a floating number between 0 and 1");
			System.out.println("Aborting program...");
			return;
		}
		
		System.out.println("Edge weight threshold = " + threshold);
		// double threshold = 1.0;
		ArrayList<Integer> idList = getIDList(prabList, n, threshold);
		ArrayList<String> tissueList = new ArrayList<String>();

		tissueList.add("prostate");
		tissueList.add("breast");
		tissueList.add("thyroid");
		tissueList.add("lung");

		String tcgaPath = MainAppl.getTCGA_Path(args);
		File[] patientPath = new File(MainAppl.getPatientSetPath(args)).listFiles();
		outputMap = cmpPatientSetTissues(tissueList, tcgaPath, idList, patientPath);

		String outputFile = MainAppl.getOutputPath(args);

		System.out.println(
				"Start writing output to file " + outputFile + "ProstateCancer_" + networkType + ".prediction1");

		writeOutput(outputFile, networkType, tissue);

		String endTime = Calendar.getInstance().getTime().toString().split("\\s+")[3];
		System.out.println("End at " + endTime);
	}

	private static final String END_TCGA_MERGED = "_MERGED_FILE_DE.tsv";

	private static final String SYNTAX_HELP = "Required command line flags:\n"
			+ "-network <GIANT or OW-GIANT network>\n"
			+ "-tissue <Tissue you want to work with, preferrably prostate_gland, must be named after RAF file without \"_all\">\n"
			+ "-raf <RAF path [provided in: /home/proj/biocluster/NEAP_SS16/woerheide/output/  ||  /home/proj/biocluster/NEAP_SS16/muellerho/networks/]>\n"
			+ "-allgenes <All genes file (must have for parsing Networks) [provided in: /home/proj/biocluster/NEAP_SS16/weberste/Solution5/all_genes.txt]>\n"
			+ "-threshold <Define threshold (0.0:1.0) for edges in networks to considere only those genes above or equal to threshold [defualt: 0.75]>\n"
			+ "-tcga <Path to merged TCGA files containing Fold changes [provided in: /home/proj/biocluster/NEAP_SS16/weberste/Solution5/TCGA/]>\n"
			+ "-patients <Path to Patient sets containing Fold changes [provided in: /home/proj/biocluster/NEAP_SS16/weberste/Solution5/PatientFCS/]>\n"
			+ "-malacards <Gene set for malacards genes related to prostate cancer which is a tab separated file with columns: UniprotID and GeneID [provided in: /home/proj/biocluster/NEAP_SS16/weberste/Solution5/MalacardsGenes/prostate_cancer_genes.txt]>\n"
			+ "-o <Output path (no file name required)>";

	private static final String EXAMPLE_INVOKE = "Try using:\n" + "java -jar task5.jar -network giant "
			+ "-tissue prostate_gland -raf /home/proj/biocluster/NEAP_SS16/woerheide/output/ "
			+ "-allgenes /home/proj/biocluster/NEAP_SS16/weberste/Solution5/all_genes.txt " + "-threshold 0.9 "
			+ "-tcga /home/proj/biocluster/NEAP_SS16/weberste/Solution5/TCGA/ "
			+ "-patients /home/proj/biocluster/NEAP_SS16/weberste/Solution5/PatientFCS/ "
			+ "-malacards /home/proj/biocluster/NEAP_SS16/weberste/Solution5/MalacardsGenes/prostate_cancer_genes.txt "
			+ "-o /path/to/outputfile/";

	private static HashMap<String, double[]> outputMap;

	/**
	 * calculate z score
	 * 
	 * @param pccs
	 * @return
	 */
	public static double[] normalize(double[] pccs) {
		double[] z = new double[pccs.length];
		z = Arrays.copyOf(pccs, pccs.length);

		double mean = 0;
		int c = 0;
		for (int i = 0; i < pccs.length; i++) {
			if (Double.isNaN(z[i])) {
				continue;
			}
			mean += z[i];
			c++;
		}

		mean /= c;
		double sdev = 0;
		for (int i = 0; i < pccs.length; i++) {
			if (Double.isNaN(z[i])) {
				continue;
			}
			z[i] -= mean;
			sdev += Math.pow(z[i], 2);

		}
		sdev = Math.sqrt(sdev / (c - 1));
		for (int i = 0; i < pccs.length; i++) {
			if (Double.isNaN(z[i])) {
				continue;
			}
			z[i] /= sdev;
		}
		return z;
	}

	/**
	 * write output
	 * 
	 * @param outputFile
	 * @param tissue
	 */
	private static void writeOutput(String outputFile, String networkType, String tissue) {
		try (BufferedWriter wr = new BufferedWriter(
				new FileWriter(outputFile +  networkType + "_ProstateCancer_PRAB.prediction1"))) {
			for (String s : outputMap.keySet()) {
				double[] pccArray = outputMap.get(s);
				double[] curZ_Score = normalize(outputMap.get(s));

				wr.write(s + "");

				double pZscore = curZ_Score[0];
				double bZscore = curZ_Score[1];
				double tZscore = curZ_Score[2];
				double lZscore = curZ_Score[3];

				double maxOther = Math.max(tZscore, Math.max(lZscore, bZscore));

				// prostate cancer prediction limitation
				if (pZscore - 1.0 > bZscore && pZscore - 1.0 > tZscore && pZscore - 1.0 > lZscore
						&& pccArray[0] > 0.0) {
					wr.write("\t" + tissue);
					wr.write("\t" + Math.round((pZscore - maxOther - 1.0) * 10000) / 10000.0);
					wr.write("\tProstate_Cancer");
				} else {
					wr.write("\t" + tissue);
					wr.write("\t0.000");
					wr.write("\tNO");
				}
				wr.newLine();
			}

			wr.close();
		} catch (FileNotFoundException eFNF) {
			eFNF.printStackTrace();
		} catch (IOException eIO) {
			eIO.printStackTrace();
		}
	}

	/**
	 * get prostate cancer malacards genes
	 * 
	 * @param prabList
	 * @param n
	 * @param threshold
	 * @return
	 */
	private static ArrayList<Integer> getIDList(ArrayList<Integer> prabList, GIANT_OW_Network n, double threshold) {
		ArrayList<Integer> idList = new ArrayList<Integer>();
		for (Integer i : prabList) {
			// for (Integer i : Network.genMap.keySet()) {
			for (Integer j : GIANT_OW_Network.genMap.keySet()) {
				double score = n.getEdge(i, j);
				if (score >= threshold) {
					if (!idList.contains(i))
						idList.add(i);

					if (!idList.contains(j))
						idList.add(j);
				}
			}
		}
		return idList;
	}

	/**
	 * compare patient sets with tissue sets (only genes from idList
	 * [malacards])
	 * 
	 * @param tissueList
	 * @param tcgaPath
	 * @param idList
	 * @param patientPath
	 * @return
	 */
	private static HashMap<String, double[]> cmpPatientSetTissues(ArrayList<String> tissueList, String tcgaPath,
			ArrayList<Integer> idList, File[] patientPath) {
		outputMap = new HashMap<String, double[]>();
		int tissuePosition = 0;

		for (String curTissue : tissueList) {
			File curTCGAFile = new File(tcgaPath + curTissue + END_TCGA_MERGED);

			HashMap<Integer, Double> tcgaMap = ParseFC_TCGA_File.readTCGA(idList, curTCGAFile);

			for (File patientFile : patientPath) {
				HashMap<Integer, Double> patientMap = ParsePatientSets.getPCCPatientSets(patientFile, idList);

				int c = 0;
				for (Integer i : patientMap.keySet()) {
					if (tcgaMap.containsKey(i)) {
						c++;
					}
				}

				double[] xsPatients = new double[c];
				int posPatients = 0;
				double[] ysTCGA = new double[c];
				int posTCGA = 0;

				for (Integer i : patientMap.keySet()) {
					if (tcgaMap.containsKey(i)) {
						xsPatients[posPatients] = patientMap.get(i);
						posPatients++;
						ysTCGA[posTCGA] = tcgaMap.get(i);
						posTCGA++;
					}
				}

				double pcc = CalcPearsonCorrelation.pcc(xsPatients, ysTCGA);

				String curPatientSet = patientFile.getName().toString().split("\\.")[0];

				// patient set ____ tissue ___ pcc
				if (!outputMap.containsKey(curPatientSet)) {
					double[] pccArray = new double[4];
					pccArray[tissuePosition] = pcc;
					outputMap.put(curPatientSet, pccArray);
				} else {
					double[] pccArray = outputMap.get(curPatientSet);
					pccArray[tissuePosition] = pcc;
					outputMap.put(curPatientSet, pccArray);
				}
			}
			System.out.println(curTissue + " currently compared");
			tissuePosition++;
		}

		return outputMap;
	}

	/**
	 * loop through command line arguments
	 * 
	 * @param args
	 *            all command line arguments given
	 * @param cmpString
	 *            compare string for each argument
	 * @return argument for given compare string
	 * @throws NullPointerException
	 */
	private static String goThroughArgs(String[] args, String cmpString) throws NullPointerException {
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals(cmpString)) {
				try {
					return args[i + 1];
				} catch (ArrayIndexOutOfBoundsException e) {
					return null;
				}
			}
		}
		return null;
	}

	/*
	 * parse command line arguments
	 */

	private static String getNetworkType(String[] args) throws NullPointerException {
		try {
			return MainAppl.goThroughArgs(args, "-network");
		} catch (Exception e) {
			return null;
		}
	}

	private static String getTissue(String[] args) throws NullPointerException {
		try {
			return MainAppl.goThroughArgs(args, "-tissue");
		} catch (Exception e) {
			return null;
		}
	}

	private static String getRAF_Path(String[] args) throws NullPointerException {
		try {
			return MainAppl.goThroughArgs(args, "-raf");
		} catch (Exception e) {
			return null;
		}
	}

	private static String getAllGenesFile(String[] args) throws NullPointerException {
		try {
			return MainAppl.goThroughArgs(args, "-allgenes");
		} catch (Exception e) {
			return null;
		}
	}

	private static String getThreshold(String[] args) throws NullPointerException {
		try {
			return MainAppl.goThroughArgs(args, "-threshold");
		} catch (Exception e) {
			return null;
		}
	}

	private static String getTCGA_Path(String[] args) throws NullPointerException {
		try {
			return MainAppl.goThroughArgs(args, "-tcga");
		} catch (Exception e) {
			return null;
		}
	}

	private static String getPatientSetPath(String[] args) throws NullPointerException {
		try {
			return MainAppl.goThroughArgs(args, "-patients");
		} catch (Exception e) {
			return null;
		}
	}

	private static String getMalacards(String[] args) throws NullPointerException {
		try {
			return MainAppl.goThroughArgs(args, "-malacards");
		} catch (Exception e) {
			return null;
		}
	}

	private static String getOutputPath(String[] args) throws NullPointerException {
		try {
			return MainAppl.goThroughArgs(args, "-o");
		} catch (Exception e) {
			return null;
		}
	}

}

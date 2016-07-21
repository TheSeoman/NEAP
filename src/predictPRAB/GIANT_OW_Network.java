package predictPRAB;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.HashMap;

public class GIANT_OW_Network {

	private boolean isInsideMap = false;

	private int geneSize;
	private ByteBuffer[] input;
	static HashMap<Integer, Integer> genMap = new HashMap<Integer, Integer>();
	
	public GIANT_OW_Network(String path, String genes) throws IOException {

		readInGenes(genes);
		geneSize = genMap.size();
		input = new ByteBuffer[geneSize];
		readRaf(path);
	}
	
	/*
	 * readInGenes: read in list of genes(nodes) present in Network and generate
	 * hashmap.
	 */

	private static void readInGenes(String genes) {

		System.out.println("Read in genes in network..");
		InputStream is = null;
		try {
			is = new FileInputStream(genes);
			BufferedReader br = new BufferedReader(new InputStreamReader(is));

			String line = br.readLine();
			int count = 0;

			while (line != null) {

				line = line.trim();
				genMap.put(Integer.parseInt(line), count);
				count++;

				line = br.readLine();
			}

			br.close();
		} catch (IOException e) {
			System.err.print("File can not be read.");
			e.printStackTrace();
		}

	}

	/*
	 * readRaf: read in as many bytes as necessary and save them in ByteBuffer
	 * array array length: number of nodes in network buffer lengths:
	 * incrementing=> half matrix
	 */

	private void readRaf(String file) throws IOException {
		RandomAccessFile raf = new RandomAccessFile(file, "r");
		FileChannel channel = raf.getChannel();

		try {
			ByteBuffer bb;
			for (int i = 0; i < geneSize; i++) {
				bb = ByteBuffer.allocate((i + 1) * 8);
				channel.read(bb);
				bb.flip();
				input[i] = bb;
				bb.clear();

			}
		} finally {

			channel.close();
			raf.close();
		}

	}

	/*
	 * getEdge: inputs two entrezIDs and returns the corresponding edge weight
	 */

	public double getEdge(Integer a, Integer b) {
		this.setInsideMap(false);
		a = genMap.get(a);
		b = genMap.get(b);
		if (a != null && b != null) {
			this.setInsideMap(true);
			if (a <= b) {
				return (input[b].getDouble((a) * 8));

			} else {

				return (input[a].getDouble((b) * 8));

			}

		}
		System.out.println("NAAAA");
		return 0.0;
	}

	public boolean isInsideMap() {
		return isInsideMap;
	}

	public void setInsideMap(boolean isInsideMap) {
		this.isInsideMap = isInsideMap;
	}

	
}

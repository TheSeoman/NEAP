import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.HashMap;

/*
 * Object network reads in a network from file "path" and uses the file "genes" to map EntrezIDs to their position.
 * Genes are given in the all_genes.txt, it is import so specifically use this file.
 * The getEdge(genA,genB) method allows direct access to the edge weight between genA and GenB. 
 */
      

public class Network {

	private int geneSize;
	private ByteBuffer[] input;
	public static HashMap<Integer, Integer> genMap = new HashMap<Integer, Integer>();

	public Network(String path, String genes) throws IOException {

		readInGenes(genes);
		geneSize=genMap.size();
		input= new ByteBuffer[geneSize];
		readRaf(path);
	}

	
	/*
	 * readInGenes: read in list of genes(nodes) present in CheckAberrantGenesPatientSet2 and generate hashmap.
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

				line=line.trim();
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
	 * readRaf: read in as many bytes as necessary and save them in ByteBuffer array 
	 * array length: number of nodes in network
	 * buffer lengths: incrementing=> half matrix 
	 */
	
	private void readRaf(String file) throws IOException {
		RandomAccessFile raf = new RandomAccessFile(file, "r");
		FileChannel channel = raf.getChannel();

		try{
			ByteBuffer bb;
			for(int i=0;i<geneSize;i++){
				bb= ByteBuffer.allocate((i+1)*8);
				channel.read(bb);
				bb.flip();
				input[i]=bb;
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
		a=genMap.get(a);
		b=genMap.get(b);

		if (a <= b){
			return(input[b].getDouble((a)*8));

		}else{

			return(input[a].getDouble((b)*8));

		}
	}

	public double getConnectivityScore(Integer u, Integer v){
		double puv = getEdge(v, u);
		double dv = 0;
		int vq = 0;
		for(int t : genMap.values()){
			double pvt = getEdge(v, t);
			if(pvt > 0.1) {
				vq++;
				dv += pvt;
			}
		}
		dv /= vq;
		return puv/dv;
	}
}

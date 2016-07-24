package patientsFC;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.HashMap;

public class TxtToRAF {

    private static final int TIMES_ALLOCATE = 8;
    private File file;
    private String output;
    private HashMap<Integer, Integer> genes;
    private ByteBuffer[] network;
    private int numberGeneInteractions;
    private double sumScore;

    public TxtToRAF(File file, String output, HashMap<Integer, Integer> genes) throws IOException {
        this.file = file;
        this.output = output;
        this.genes = genes;
        createBB();
        readIn();
        writeRaf();
    }

    private void createBB() {
        network = new ByteBuffer[genes.size()];
        System.out.println("SIZE = "+genes.size()+" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        System.out.println("Create ByteBuffer");
        for (int i = 0; i < genes.size(); i++) {
            ByteBuffer bb = ByteBuffer.allocate((i+1) * TIMES_ALLOCATE);
            network[i] = bb;
        }
        System.out.println("BB created");
    }

    private void writeRaf() throws IOException {
        RandomAccessFile raf = new RandomAccessFile(output, "rw");

        FileChannel channel = raf.getChannel();

        System.out.println("RAF start...");

        try {
            for (int i = 0; i < genes.size(); i++) {
                network[i].clear();
                channel.write(network[i]);
            }
        } finally {
            if (channel != null)
                channel.close();
            if (raf != null)
                raf.close();
        }

        System.out.println("RAF successfully ended...");
    }

    private void addEdge(Integer a, Integer b, Double weight) {
        System.out.println("a = "+a+"\t"+"b = "+b+"\t"+" weight = "+weight);
        if (a <= b) {
            network[b].position((a) * TIMES_ALLOCATE);
            network[b].putDouble(weight);
        } else {
            network[a].position((b) * TIMES_ALLOCATE);
            network[a].putDouble(weight);
        }
    }

    private void readIn() {
        try (BufferedReader bur = new BufferedReader(
                new FileReader(file))) {

            String sLine = null;

            System.out.println(file.getName()+" currently read in...");

            while ((sLine = bur.readLine()) != null) {
                String[] entries = sLine.split("\\t");
                System.out.println("entries0 = "+entries[0]+"\tentries1 ="+entries[1]+"\tentries2 = "+entries[2]);
                System.out.println(Integer.parseInt(entries[0]));
                addEdge(genes.get(Integer.parseInt(entries[0])), genes.get(Integer.parseInt(entries[1])),
                        (Double.parseDouble(entries[2])));
            }

            bur.close();
        } catch (FileNotFoundException eFNF) {
            eFNF.printStackTrace();
        } catch (IOException eIO) {
            eIO.printStackTrace();
        }
    }

    public int getNumberGeneInteractions() {
        return numberGeneInteractions;
    }

    public void setNumberGeneInteractions(int numberGeneInteractions) {
        this.numberGeneInteractions = numberGeneInteractions;
    }

    public double getSumScore() {
        return sumScore;
    }

    public void setSumScore(double sumScore) {
        this.sumScore = sumScore;
    }

}

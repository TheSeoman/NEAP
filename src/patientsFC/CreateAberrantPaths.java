package patientsFC;

import java.util.*;

/**
 * Created by Stefan on 23.07.2016.
 */
public class CreateAberrantPaths {

    private static HashSet<LinkedList<Integer>> allPaths;
    private static LinkedList<Integer> visited;

    /**
     * @param aberrant                aberrant genes map
     * @param neighbors               neighbors map
     * @param minimumSizeSubnetwork   determines minimum number of nodes on path
     * @param percentageAberrantGenes determines percentage of aberrant genes of patients in gene
     * @param tcgaCases               number of tcga cases to compare with
     * @param threshold               threshold for aberrant patient genes
     * @return all paths with at least minimumSizeSubnetwork nodes and at least percentageAberrantGenes% aberrant genes
     */
    public static HashMap<LinkedList<Integer>, int[]> calcUnion(HashMap<Integer, int[]> aberrant, HashMap<Integer, HashSet<Integer>> neighbors, int minimumSizeSubnetwork, double percentageAberrantGenes, int tcgaCases, double threshold) {
        allPaths = new HashSet<LinkedList<Integer>>();
        HashMap<LinkedList<Integer>, Boolean> o = new HashMap<LinkedList<Integer>, Boolean>();

        for (int x = minimumSizeSubnetwork; x < neighbors.size(); x++) {
            for (Integer i : neighbors.keySet()) {
                for (Integer j : neighbors.keySet()) {
                    visited = new LinkedList<Integer>();
                    findAllPossiblePaths(neighbors, aberrant, i, j, x, threshold);
                }
            }
            // unique paths
            for (LinkedList<Integer> l : allPaths) {
                o.put(l, true);
            }
        }

        HashMap<LinkedList<Integer>, int[]> output = new HashMap<LinkedList<Integer>, int[]>();

        for (LinkedList<Integer> l : o.keySet()) {
            int[] union = new int[tcgaCases];
            for (int i = 0; i < union.length; i++) {
                union[i] = 1;
            }

            int countAberrantGenes = 0;

            for (Integer i : l) {
                int[] a = aberrant.get(i);
                countAberrantGenes = 0;
                int[] temp = new int[tcgaCases];
                for (int j = 0; j < a.length; j++) {
                    temp[j] = (a[j] == 1 && union[j] == 1) ? 1 : 0;
                    if (temp[j] == 1) {
                        countAberrantGenes++;
                    }
                }
                System.arraycopy(temp, 0, union, 0, union.length);
            }
            if (((double) countAberrantGenes / (double) tcgaCases) >= percentageAberrantGenes) {
                output.put(l, union);
            }
        }

        for (LinkedList<Integer> l : output.keySet()) {
            int[] a = output.get(l);

            System.out.print(l+":");
            int c = 0;
            for (int i = 0; i < a.length; i++) {
                if (a[i] == 1) c++;
            }
            System.out.print("("+c+":"+tcgaCases+")  ");
            System.out.println((100 * ((double) c / (double) tcgaCases))+"% of aberrant patient genes...");
        }

        return output;
    }

    /**
     * find all possible paths from a to b with exactly sizeOfMcSub (Integer) nodes
     *
     * @param m           neighbors map
     * @param a           node source
     * @param b           node destination
     * @param sizeOfMcSub number of nodes in path
     */
    @SuppressWarnings("Since15")
    private static void findAllPossiblePaths(HashMap<Integer, HashSet<Integer>> m, HashMap<Integer, int[]> aberrantMap, int a, int b, int sizeOfMcSub, double threshold) {
        int[] curAberrant = aberrantMap.get(a);
        int counter = 0;
        for (int i = 0; i < curAberrant.length; i++) {
            if (curAberrant[i] == 1)
                counter++;
        }
        if ((double) counter / (double) curAberrant.length < threshold) {
            return;
        }
        if (!visited.contains(a))
            visited.add(a);
        for (Integer i : m.get(a)) {
            curAberrant = aberrantMap.get(a);
            counter = 0;
            for (int k = 0; k < curAberrant.length; k++) {
                if (curAberrant[k] == 1)
                    counter++;
            }
            if ((double) counter / (double) curAberrant.length < threshold) {
                return;
            }
            if (visited.contains(i)) {
                continue;
            }
            if (i == b) {
                visited.add(i);
                if (visited.size() == sizeOfMcSub) {
                    visited.sort(new Comparator<Integer>() {
                        @Override
                        public int compare(Integer x, Integer y) {
                            return x < y ? -1 : x == y ? 0 : 1;
                        }
                    });
                    visited.addAll(0, visited);
                    visited.removeLast();
                    if (!allPaths.contains(visited)) {
                        allPaths.add(visited);
                    }
                }
                visited.removeLast();
                break;
            }
        }

        for (Integer i : m.get(a)) {
            if (visited.contains(i) || i == b) {
                continue;
            }
            visited.addLast(i);
            findAllPossiblePaths(m, aberrantMap, i, b, sizeOfMcSub, threshold);
            visited.removeLast();
        }
    }

}

import java.util.Set;

/**
 * Created by seoman on 7/24/16.
 */
public class AbberrantNeighbor {
    Network net;
    Set<Integer> genesOfInterest;

    public AbberrantNeighbor(Network net, Set<Integer> genesOfInterest) {
        this.net = net;
        this.genesOfInterest = genesOfInterest;
    }

    public double testAssociation(Set<Integer> abberrantGenes) {
        double totalAssociation = 0;
        int c = 0;
        for (int id1 : abberrantGenes) {
            if (net.genMap.containsKey(id1)) {
                double geneAssociation = 0.0;
                for (int id2 : genesOfInterest) {
                    if (net.genMap.containsKey(id2)) {
                        double edge = net.getEdge(id1, id2);
                        if(edge > 0.1)
                            geneAssociation = geneAssociation + (1 - geneAssociation) * Math.pow(edge, 2);
                    }
                }
                totalAssociation += geneAssociation;
                c++;
            }
        }
        return totalAssociation / c;
    }
}

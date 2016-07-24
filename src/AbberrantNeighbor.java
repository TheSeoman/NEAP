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
        for (int id1 : abberrantGenes) {
            double geneAssociation = 0.0;
            for (int id2 : genesOfInterest) {
                double edge = net.getEdge(id1, id2);
                geneAssociation = geneAssociation + (1 - geneAssociation) * edge;
            }
        }
        return totalAssociation / abberrantGenes.size();
    }
}

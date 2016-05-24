import java.util.Map;

/**
 * Created by schmidtju on 09.05.16.
 */
public class Network {
    private byte[][] net;
    private Map<Integer, Integer> idMap;

    public Network(byte[][] net, Map<Integer, Integer> idMap){
        this.net = net;
        this.idMap = idMap;
    }

    public byte[][] getNet() {
        return net;
    }

    public Map<Integer, Integer> getIdMap() {
        return idMap;
    }
}

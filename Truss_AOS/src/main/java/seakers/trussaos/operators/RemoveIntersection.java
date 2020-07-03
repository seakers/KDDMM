package seakers.trussaos.operators;

import seakers.aos.operator.CheckParents;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;
import java.util.ArrayList;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;

/**
 * Repair operator that improves design feasibility. Two operating modes are available depending on whether design
 * stability is also to be considered
 *
 * @author roshan94
 */

public class RemoveIntersection implements Variation, CheckParents {

    /**
     * Boolean value determining whether to consider stability while employing operator
     */
    private final boolean KeepStable;

    //private final TrussRepeatableArchitecture architecture;

    //private final int[][] FullConnectivityArray;

    private final double[][] NodalConnectivityArray = {{}};

    /**
     * Constructor for RemoveIntersection class
     * @param KeepStable
     */
    public RemoveIntersection(boolean KeepStable){
        this.KeepStable = KeepStable;
    }

    @Override
    public int getArity() {return 1;}

    @Override
    public Solution[] evolve(Solution[] sols) {
        TrussRepeatableArchitecture architecture = (TrussRepeatableArchitecture) sols[0];
        TrussRepeatableArchitecture architrectureCopy = (TrussRepeatableArchitecture) architecture.copy();

        if(KeepStable) {

        }
        else {

        }

    }

    private int[] findIntersectingTrusses (double[][] NodeLocations) {


    }

    private int findOrientation (double[] point1, double[] point2, double[] point3){
        double slopeNumerator1 = ((point2[1] - point1[1])*(point3[0] - point2[0]));
        double slopeNumerator2 = ((point3[1] - point2[1])*(point2[0] - point1[0]));
        return Double.compare(slopeNumerator1, slopeNumerator2);
    }

    private boolean determineIntersection(double[][] line1, double[][] line2){

    }

    @Override
    public boolean check(Solution[] sols) {

    }
}






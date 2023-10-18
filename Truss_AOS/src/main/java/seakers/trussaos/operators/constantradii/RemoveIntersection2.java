package seakers.trussaos.operators.constantradii;

import com.mathworks.engine.MatlabEngine;
//import seakers.aos.operator.CheckParents;
import org.moeaframework.core.PRNG;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;
//import java.util.ArrayList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;

//import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;

/**
 * Repair operator for the constant radii problem that improves design feasibility. This version is a greedier operator than the previous
 * version "RemoveIntersection"
 *
 * @author roshan94
 */

public class RemoveIntersection2 implements Variation {

    private final boolean arteryProblem;
    private static MatlabEngine engine;
    private final double[][] nodalConnectivityArray;
    private final double sidenum;
    private final double sel;
    private final int numHeurObjectives;
    private final int numHeurConstraints;

    public RemoveIntersection2(boolean arteryProblem, MatlabEngine eng, double[][] nodalConnArray, double sidenum, double sel, int numHeuristicObjectives, int numHeuristicConstraints) {
        this.arteryProblem = arteryProblem;
        engine = eng;
        this.nodalConnectivityArray = nodalConnArray;
        this.sidenum = sidenum;
        this.sel = sel;
        this.numHeurObjectives = numHeuristicObjectives;
        this.numHeurConstraints = numHeuristicConstraints;
    }

    @Override
    public int getArity() {
        return 1;
    }

    @Override
    public Solution[] evolve(Solution[] sols) {
        TrussRepeatableArchitecture architecture = new TrussRepeatableArchitecture(sols[0],sidenum,numHeurObjectives,numHeurConstraints);
        double[][] connectivityArray = architecture.getConnectivityArrayFromSolution(sols[0]);
        TrussRepeatableArchitecture newArchitecture;
        double[][] newConnectivityArray;

        Object output = null;
        try {
            output = engine.feval("intersectLogger",nodalConnectivityArray,connectivityArray);
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
        double[][] intersectionMatrix = (double[][])output;
        if (intersectionMatrix == null) {
            newConnectivityArray = connectivityArray.clone();
        } else {
            double[][] intersectingMembers = new double[intersectionMatrix.length][2];
            double[] numberOfIntersections = new double[intersectionMatrix.length];

            for (int i = 0; i < intersectingMembers.length; i++) {
                intersectingMembers[i][0] = intersectionMatrix[i][0];
                intersectingMembers[i][1] = intersectionMatrix[i][1];
                numberOfIntersections[i] = intersectionMatrix[i][2];
            }

            // Store indices of intersecting members for the corresponding number of intersections
            HashMap<Double, ArrayList<Integer>> indices = new HashMap<>();
            for (int i = 0; i < numberOfIntersections.length; i++) {
                indices.computeIfAbsent(numberOfIntersections[i], c -> new ArrayList<>()).add(i);
            }

            // Find member with maximum number of intersections
            ArrayList<Integer> maxIntersectionIndices = indices.get(Arrays.stream(numberOfIntersections).max().getAsDouble());
            double memberIndexToDelete = maxIntersectionIndices.get(PRNG.nextInt(maxIntersectionIndices.size()));
            double[] memberToDelete = new double[2];
            System.arraycopy(intersectingMembers[(int) memberIndexToDelete], 0, memberToDelete, 0, 2);

            // Remove member from original design
            newConnectivityArray = new double[connectivityArray.length-1][2];
            int next = 0;
            for (int i = 0; i < connectivityArray.length; i++) {
                if (connectivityArray[i][0] == (int) memberToDelete[0]) {
                    if (connectivityArray[i][1] == (int) memberToDelete[1]) {
                        continue;
                    }
                }
                System.arraycopy(connectivityArray[i], 0, newConnectivityArray[next], 0, 2);
                next += 1;
            }
        }
        newArchitecture = architecture.getArchitectureFromConnectivityArray(newConnectivityArray, arteryProblem);

        return new Solution[]{newArchitecture};
    }
}

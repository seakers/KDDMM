package seakers.trussaos.operators;

import com.mathworks.engine.MatlabEngine;
//import seakers.aos.operator.CheckParents;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;
//import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.ExecutionException;

//import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;

/**
 * Repair operator that improves design feasibility. Two operating modes are available depending on whether design
 * stability is also to be considered
 *
 * @author roshan94
 */

public class RemoveIntersection implements Variation {

    /**
     * Boolean value determining whether to consider stability while employing operator
     */
    private final boolean KeepStable;

    private static MatlabEngine engine;

    //private final TrussRepeatableArchitecture architecture;

    //private final int[][] FullConnectivityArray;

    private final double[][] nodalConnectivityArray;

    /**
     * Constructor for RemoveIntersection class
     * @param KeepStable
     * @param nodalConnArray
     * @param eng
     */
    public RemoveIntersection(boolean KeepStable, MatlabEngine eng, double[][] nodalConnArray){
        this.KeepStable = KeepStable;
        engine = eng;
        this.nodalConnectivityArray = nodalConnArray;
    }

    @Override
    public int getArity() {return 1;}

    @Override
    public Solution[] evolve(Solution[] sols) {
        TrussRepeatableArchitecture architecture = new TrussRepeatableArchitecture(sols[0]);
        //TrussRepeatableArchitecture architectureCopy = new TrussRepeatableArchitecture(architecture.deepCopy());
        int[][] connectivityArray = architecture.getConnectivityArrayFromSolution(sols[0]);
        try {
            int[][] intersectingTrusses = findIntersectingTrusses(connectivityArray);
        } catch (ExecutionException | InterruptedException e) {
            e.printStackTrace();
        }
        TrussRepeatableArchitecture newArchitecture;
        int[][] newConnectivityArray = new int[connectivityArray.length-1][2];
        if(KeepStable) {
            int[] oldNumberOfConnections = getNumberOfConnections(connectivityArray);
            boolean stabilitySameOrBetter = false;
            int numberOfAttempts = 0;
            while (!stabilitySameOrBetter && numberOfAttempts < 5) {
                try {
                    newConnectivityArray = removeIntersection(connectivityArray);
                } catch (ExecutionException | InterruptedException e) {
                    e.printStackTrace();
                }
                int[] newNumberOfConnections = getNumberOfConnections(newConnectivityArray);
                boolean[] connectionsImproved = compareNumberOfConnections(oldNumberOfConnections,newNumberOfConnections);
                for (boolean b : connectionsImproved) {
                    stabilitySameOrBetter = b;
                }
                numberOfAttempts += 1;
            }
            newArchitecture = architecture.getArchitectureFromConnectivityArray(newConnectivityArray);
        }
        else {
            try {
                newConnectivityArray = removeIntersection(connectivityArray);
            } catch (ExecutionException | InterruptedException e) {
                e.printStackTrace();
            }
            newArchitecture = architecture.getArchitectureFromConnectivityArray(newConnectivityArray);
        }
        return new Solution[]{newArchitecture};
    }

    private boolean[] compareNumberOfConnections (int[] oldConnections, int[] newConnections) {
        int[] oldMidConnections = new int[4];
        int[] oldCornerConnections = new int[4];
        int[] newMidConnections = new int[4];
        int[] newCornerConnections = new int[4];
        System.arraycopy(oldConnections,5,oldMidConnections,0,4);
        System.arraycopy(oldConnections,1,oldCornerConnections,0,4);
        System.arraycopy(newConnections,5,newMidConnections,0,4);
        System.arraycopy(newConnections,1,oldCornerConnections,0,4);
        boolean cornerConnectionsImproved = false;
        boolean midPointConnectionsImproved = false;
        for (int i : oldCornerConnections) {
            for (int j : newCornerConnections) {
                if (j>=i && j>=2) {
                    cornerConnectionsImproved = true;
                }
                else {
                    cornerConnectionsImproved = false;
                }
            }
        }
        for (int i : oldMidConnections) {
            for (int j : newMidConnections) {
                if (j>=i && j>=3) {
                    midPointConnectionsImproved = true;
                }
                else {
                    midPointConnectionsImproved = false;
                }
            }
        }
        return new boolean[] {cornerConnectionsImproved,midPointConnectionsImproved};
    }

    private int[] getNumberOfConnections (int[][] designConnectivityArray) {
        int[] numberOfMidPointConnections = new int[4]; // Nodes = [2, 4, 6, 8]
        int[] numberOfCornerConnections = new int[4]; // Nodes = [1, 3, 7, 9]
        int[] numberOfCenterConnections = new int[1]; // Nodes = 5
        for (int i = 0; i < designConnectivityArray.length; i++) {
            for (int j = 0; j < 2; j++) {
                if (designConnectivityArray[i][j] == 1) {
                    numberOfCornerConnections[0] += 1;
                }
                else if (designConnectivityArray[i][j] == 2) {
                    numberOfMidPointConnections[0] += 1;
                }
                else if (designConnectivityArray[i][j] == 3) {
                    numberOfCornerConnections[1] += 1;
                }
                else if (designConnectivityArray[i][j] == 4) {
                    numberOfMidPointConnections[1] += 1;
                }
                else if (designConnectivityArray[i][j] == 5) {
                    numberOfCenterConnections[0] += 1;
                }
                else if (designConnectivityArray[i][j] == 6) {
                    numberOfMidPointConnections[2] += 1;
                }
                else if (designConnectivityArray[i][j] == 7) {
                    numberOfCornerConnections[2] += 1;
                }
                else if (designConnectivityArray[i][j] == 8) {
                    numberOfMidPointConnections[3] += 1;
                }
                else if (designConnectivityArray[i][j] == 9) {
                    numberOfCornerConnections[3] += 1;
                }
            }
        }
        int[] numberOfConnections = new int[numberOfCenterConnections.length + numberOfCornerConnections.length + numberOfMidPointConnections.length];
        System.arraycopy(numberOfCenterConnections,0,numberOfConnections,0,numberOfCenterConnections.length);
        System.arraycopy(numberOfCornerConnections,0,numberOfConnections,numberOfCenterConnections.length,numberOfCornerConnections.length);
        System.arraycopy(numberOfMidPointConnections,0,numberOfConnections,numberOfCornerConnections.length,numberOfMidPointConnections.length);
        return numberOfConnections;
    }

    private int[][] removeIntersection(int[][] designConnectivityArray) throws ExecutionException, InterruptedException {
        int[][] trussIntersections = findIntersectingTrusses(designConnectivityArray);
        int[][] newConnectivityArray;
        if (trussIntersections.length == 0) {
            newConnectivityArray = designConnectivityArray.clone();
        }
        else {
            Random ran = new Random();
            int intersectionChoice = ran.nextInt(trussIntersections.length);
            int[] intersectingTrussIndexPair = trussIntersections[intersectionChoice];
            int trussIndexChoiceToDelete = ran.nextInt(intersectingTrussIndexPair.length);
            int trussIndexToDelete = intersectingTrussIndexPair[trussIndexChoiceToDelete];
            newConnectivityArray = new int[designConnectivityArray.length-1][2];
            int next = 0;
            for (int i = 0; i < designConnectivityArray.length; i++) {
                if (i == trussIndexToDelete) {
                    continue;
                }
                System.arraycopy(designConnectivityArray[i], 0, newConnectivityArray[next], 0, 2);
                next += 1;
            }
        }
        return newConnectivityArray;
    }

    private int[][] findIntersectingTrusses (int[][] designConnectivityArray) throws ExecutionException, InterruptedException {
        int numberOfTrusses = designConnectivityArray.length;
        int[][] intersectingTrusses = new int[numberOfTrusses*(numberOfTrusses-1)/2][2];
        double[][] currentLineNodePositions = new double[2][nodalConnectivityArray[0].length];
        double[][] nextLineNodePositions = new double[2][nodalConnectivityArray[0].length];
        int numberOfIntersectingTrussPairs = 0;
        int[] currentTrussPair = new int[2];
        int[] nextTrussPair = new int[2];
        for (int i = 0; i < designConnectivityArray.length - 1; i++) {
            for (int j = i+1; j < designConnectivityArray.length; j++) {
                for (int k = 0; k < 2; k++) {
                    currentTrussPair[k] = designConnectivityArray[i][k];
                    nextTrussPair[k] = designConnectivityArray[j][k];
                    for (int l = 0; l < nodalConnectivityArray[0].length; l++) {
                        currentLineNodePositions[k][l] = nodalConnectivityArray[currentTrussPair[k]-1][l];
                        nextLineNodePositions[k][l] = nodalConnectivityArray[nextTrussPair[k]-1][l];
                    }
                }
                boolean linesIntersect = determineIntersection(currentLineNodePositions,nextLineNodePositions);
                if (linesIntersect) {
                    intersectingTrusses[numberOfIntersectingTrussPairs][0] = i;
                    intersectingTrusses[numberOfIntersectingTrussPairs][1] = j;
                    numberOfIntersectingTrussPairs += 1;
                }
            }
        }
        int[][] trueIntersectingTrusses = new int[numberOfIntersectingTrussPairs][2];
        for (int i = 0; i < numberOfIntersectingTrussPairs; i++) {
            for (int j = 0; j < 2; j++) {
                trueIntersectingTrusses[i][j] = intersectingTrusses[i][j];
            }
        }
        return trueIntersectingTrusses;
    }

    private boolean determineIntersection(double[][] line1, double[][] line2) throws ExecutionException, InterruptedException {
        // line 1 = [p1_x, p1_y; p2_x, p2_y], line 2 = [q1_x, q1_y; q2_x, q2_y]
        double[] p1 = new double[line1[0].length];
        double[] p2 = new double[line1[0].length];
        double[] q1 = new double[line1[0].length];
        double[] q2 = new double[line1[0].length];
        for (int i = 0; i < line1[0].length; i++) {
            p1[i] = line1[0][i];
            p2[i] = line1[1][i];
            q1[i] = line2[0][i];
            q2[i] = line2[1][i];
        }
        //boolean intersects;
        Object obj = null;
        obj = engine.feval("findLineSegIntersection",p1,p2,q1,q2);
        return (boolean)obj;
    }
}






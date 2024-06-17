package seakers.trussaos.operators.integerradii;

import com.mathworks.engine.MatlabEngine;
//import seakers.aos.operator.CheckParents;
import seakers.trussaos.architecture.IntegerRepeatableArchitecture;
//import java.util.ArrayList;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;

import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;

public class RemoveIntersectionIntegerRadii implements Variation{

    private final boolean KeepStable;
    private static MatlabEngine engine;
    private final double[][] nodalConnectivityArray;
    private final double sidenum;
    private final double sel;
    private final double[] radii;
    private final double[] youngsModulii;
    private final int numHeurObjectives;
    private final int numHeurConstraints;

    public RemoveIntersectionIntegerRadii(boolean KeepStable, MatlabEngine eng, double[][] nodalConnArray, double sidenum, double sel, double[] radii, double[] youngsModulii, int numHeurObjectives, int numHeurConstraints){
        this.KeepStable = KeepStable;
        engine = eng;
        this.nodalConnectivityArray = nodalConnArray;
        this.sidenum = sidenum;
        this.sel = sel;
        this.radii = radii;
        this.youngsModulii = youngsModulii;
        this.numHeurObjectives = numHeurObjectives;
        this.numHeurConstraints = numHeurConstraints;
    }

    @Override
    public int getArity() {
        return 1;
    }

    @Override
    public Solution[] evolve(Solution[] solutions) {
        IntegerRepeatableArchitecture architecture = new IntegerRepeatableArchitecture(solutions[0], (int) sidenum, numHeurObjectives, numHeurConstraints, radii, youngsModulii);
        double[][] connectivityArray = architecture.getFullCAFromSolutionIntegerRadii();
        double[][] completeConnectivityArray = architecture.getCompleteConnectivityArrayFromSidenum();
        double[] radiusArray = architecture.getFullRadiusArray();
        int youngsModulusChoice = architecture.getYoungsModulusChoice();
        double[] newRadiusArray = new double[radiusArray.length];

        IntegerRepeatableArchitecture newArchitecture;
        double[][] newConnectivityArray = new double[connectivityArray.length-1][2];
        if(KeepStable) {
            //int[] oldNumberOfConnections = getNumberOfConnections(connectivityArray);
            double[] oldNumberOfConnections = new double[(int) (sidenum*sidenum)];
            try {
                oldNumberOfConnections = getNumberOfConnectionsRepeatable(connectivityArray);
            } catch (ExecutionException | InterruptedException e) {
                e.printStackTrace();
            }
            boolean stabilitySameOrBetter = false;
            int numberOfAttempts = 0;
            double[][] trussIntersections = null;
            double deletionTrussIndex = 0;
            while (!stabilitySameOrBetter && numberOfAttempts < 5) {
                try {
                    trussIntersections = findIntersectingTrusses(connectivityArray);
                } catch (ExecutionException | InterruptedException e) {
                    e.printStackTrace();
                }
                assert trussIntersections != null;
                if (trussIntersections.length == 0) {
                    newRadiusArray = radiusArray.clone();
                    newConnectivityArray = connectivityArray.clone();
                } else {
                    try {
                        deletionTrussIndex = getTrussIndexToDelete(connectivityArray, completeConnectivityArray, trussIntersections);
                    } catch (InterruptedException | ExecutionException e) {
                        e.printStackTrace();
                    }
                    newRadiusArray = removeRadiusMemberFromArray((int) deletionTrussIndex, radiusArray);
                    try {
                        newConnectivityArray = getNewConnectivityArray(connectivityArray,deletionTrussIndex);
                    } catch (ExecutionException | InterruptedException e) {
                        e.printStackTrace();
                    }

                    //int[] newNumberOfConnections = getNumberOfConnections(newConnectivityArray);
                    double[] newNumberOfConnections = new double[(int) (sidenum*sidenum)];
                    try {
                        newNumberOfConnections = getNumberOfConnectionsRepeatable(newConnectivityArray);
                    } catch (ExecutionException | InterruptedException e) {
                        e.printStackTrace();
                    }
                    boolean[] connectionsImproved = compareNumberOfConnections(oldNumberOfConnections,newNumberOfConnections);
                    for (boolean b : connectionsImproved) {
                        stabilitySameOrBetter = b;
                    }
                    numberOfAttempts += 1;
                }

            }
            newArchitecture = architecture.getIntegerArchitectureFromRadiusArray(newRadiusArray, youngsModulusChoice);
        }
        else {
            double[][] trussIntersections = null;
            double deletionTrussIndex = 0;
            try {
                trussIntersections = findIntersectingTrusses(connectivityArray);
            } catch (ExecutionException | InterruptedException e) {
                e.printStackTrace();
            }
            assert trussIntersections != null;
            if (trussIntersections.length == 0) {
                newRadiusArray = radiusArray.clone();
                //newConnectivityArray = connectivityArray.clone();
            } else {
                try {
                    deletionTrussIndex = getTrussIndexToDelete(connectivityArray, completeConnectivityArray, trussIntersections);
                } catch (InterruptedException | ExecutionException e) {
                    e.printStackTrace();
                }
                newRadiusArray = removeRadiusMemberFromArray((int) deletionTrussIndex, radiusArray);
                //try {
                //newConnectivityArray = getNewConnectivityArray(connectivityArray, deletionTrussIndex);
                //} catch (ExecutionException | InterruptedException e) {
                //e.printStackTrace();
                //}
            }
            newArchitecture = architecture.getIntegerArchitectureFromRadiusArray(newRadiusArray, youngsModulusChoice);

        }
        return new Solution[]{newArchitecture};
    }

    private boolean[] compareNumberOfConnections (double[] oldConnections, double[] newConnections) {
        // only for 3x3 node grid case
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
                if (j>=i && j>=3) {
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

    private double[] getNumberOfConnectionsRepeatable (double[][] designConnectivityArray) throws ExecutionException, InterruptedException {
        Object[] outputs = null;
        outputs = engine.feval("connectivityCounter",sidenum,designConnectivityArray,nodalConnectivityArray,sel);
        return (double[])outputs[0];
    }

    private double getTrussIndexToDelete (double[][] designConnectivityArray, double[][] fullConnectivityArray, double[][] intersectingTrusses) throws ExecutionException, InterruptedException {
        //double[][] trussIntersections = findIntersectingTrusses(designConnectivityArray);
        int intersectionChoice = PRNG.nextInt(intersectingTrusses.length);
        double[] intersectingTrussIndexPair = intersectingTrusses[intersectionChoice];
        int trussIndexChoiceToDelete = PRNG.nextInt(intersectingTrussIndexPair.length);
        double trussIndexToRemove = intersectingTrussIndexPair[trussIndexChoiceToDelete];
        double[] trussPairToRemove = designConnectivityArray[(int) trussIndexToRemove];
        double memberPosition = 0;
        for (int i = 0; i < fullConnectivityArray.length; i++) {
            if (fullConnectivityArray[i][0] == trussPairToRemove[0]) {
                if (fullConnectivityArray[i][1] > trussPairToRemove[1])
                    break;
            }
            if (fullConnectivityArray[i][0] > trussPairToRemove[0])
                break;
            memberPosition = i+1;
        }
        memberPosition = Math.min(memberPosition, (fullConnectivityArray.length-1));
        return memberPosition-1;
    }

    private double[][] getNewConnectivityArray(double[][] designConnectivityArray, double trussIndexToRemove) throws ExecutionException, InterruptedException {
        double[][] newConnectivityArray;
        newConnectivityArray = new double[designConnectivityArray.length-1][2];
        int next = 0;
        for (int i = 0; i < designConnectivityArray.length; i++) {
            if (i == trussIndexToRemove) {
                continue;
            }
            System.arraycopy(designConnectivityArray[i], 0, newConnectivityArray[next], 0, 2);
            next += 1;
        }
        return newConnectivityArray;
    }

    private double[][] findIntersectingTrusses (double[][] designConnectivityArray) throws ExecutionException, InterruptedException {
        int numberOfTrusses = designConnectivityArray.length;
        //int[][] intersectingTrusses = new int[numberOfTrusses*(numberOfTrusses-1)/2][2];
        ArrayList<ArrayList<Double>> intersectingTrusses = new ArrayList<>();
        double[][] currentLineNodePositions = new double[2][nodalConnectivityArray[0].length];
        double[][] nextLineNodePositions = new double[2][nodalConnectivityArray[0].length];
        //int numberOfIntersectingTrussPairs = 0;
        double[] currentTrussPair = new double[2];
        double[] nextTrussPair = new double[2];
        for (int i = 0; i < designConnectivityArray.length - 1; i++) {
            for (int j = i+1; j < designConnectivityArray.length; j++) {
                for (int k = 0; k < 2; k++) {
                    currentTrussPair[k] = designConnectivityArray[i][k];
                    nextTrussPair[k] = designConnectivityArray[j][k];
                    for (int l = 0; l < nodalConnectivityArray[0].length; l++) {
                        currentLineNodePositions[k][l] = nodalConnectivityArray[(int) currentTrussPair[k]-1][l];
                        nextLineNodePositions[k][l] = nodalConnectivityArray[(int) nextTrussPair[k]-1][l];
                    }
                }
                boolean linesIntersect = determineIntersection(currentLineNodePositions,nextLineNodePositions);
                if (linesIntersect) {
                    ArrayList<Double> intersectingPair = new ArrayList<>();
                    intersectingPair.add((double) i);
                    intersectingPair.add((double) j);
                    intersectingTrusses.add(intersectingPair);
                    //numberOfIntersectingTrussPairs += 1;
                }
            }
        }
        return convertArrayListTo2DArray(intersectingTrusses);
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

    private double[][] convertArrayListTo2DArray (ArrayList<ArrayList<Double>> arrayList2D) {
        double[][] array2D = new double[arrayList2D.size()][2];
        for (int i = 0; i < arrayList2D.size(); i++) {
            ArrayList<Double> row = arrayList2D.get(i);
            double[] rowArray = row.stream().mapToDouble(j->j).toArray();
            array2D[i] = rowArray;
        }
        return array2D;
    }

    private double[] removeRadiusMemberFromArray (int memberToRemoveIndex, double[] oldRadiusArray) {
        memberToRemoveIndex = Math.min(memberToRemoveIndex, (oldRadiusArray.length - 1));
        oldRadiusArray[memberToRemoveIndex] = 0;
        return oldRadiusArray;
    }
}

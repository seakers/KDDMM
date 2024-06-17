package seakers.trussaos.architecture;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.RealVariable;
import seakers.trussaos.problems.VariableRadiusTrussProblem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Design class for the real encoded design vector (variable radii) in the 2D NxN nodal grid case
 *
 * @author roshan94
 */

public class VariableRadiiRepeatableArchitecture extends Solution{
    private static final long serialVersionUID = -3246138689443538032L;
    private int numberOfTrusses;
    private double[] radiusArray;
    private double[] radiusLowerBounds;
    private double[] radiusUpperBounds;
    private double sidenum;
    private int[][] completeConnectivityArray;
    private int[][] designConnectivityArray;
    private static double lowestRadiusFactor;

    public VariableRadiiRepeatableArchitecture(Solution solution, double sidenum, double[] radiusLowerBounds, double[] radiusUpperBounds, double smallestRadiusFactor) {
        super(solution);
        this.sidenum = sidenum;
        lowestRadiusFactor = smallestRadiusFactor;

        // Generate Complete Connectivity Array from sidenum
        completeConnectivityArray = getCompleteConnectivityArrayFromSidenum();

        // Extract the radius array from the solution
        radiusArray = getCompleteRadiusArrayFromSolution(solution);
        this.radiusLowerBounds = radiusLowerBounds;
        this.radiusUpperBounds = radiusUpperBounds;

        // Obtain design connectivity array
        designConnectivityArray = getFullConnectivityArrayVariableRadii();
        numberOfTrusses = getTrussCountVariableRadii(radiusArray);
    }

    private double[] getRadiusArray(Solution soln) {
        int numVars = soln.getNumberOfVariables();
        double[] radii = new double[numVars];
        for (int index = 0; index < numVars; index++){
            double decision =  EncodingUtils.getReal(soln.getVariable(index));
            radii[index] = decision;
        }
        return radii;
    }

    public double[] getCompleteRadiusArrayFromSolution(Solution soln){
        return getCompleteRadiusArrayFromSolutionNxN(soln);
    }

    public int[][] getCompleteConnectivityArrayFromSidenum(){
        int memberCount = 0;
        int totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sidenum*sidenum))/(CombinatoricsUtils.factorial((int) ((sidenum*sidenum) - 2)) * CombinatoricsUtils.factorial(2)));
        int[][] completeConnectivityArray = new int[totalNumberOfMembers][2];
        for (int i = 0; i < ((sidenum*sidenum)-1); i++) {
            for (int j = i+1; j < (sidenum*sidenum); j++) {
                completeConnectivityArray[memberCount][0] = i+1;
                completeConnectivityArray[memberCount][1] = j+1;
                memberCount += 1;
            }
        }
        return completeConnectivityArray;
    }

    private int[][] getFullConnectivityArrayVariableRadii() {
        //int TrussCount = getTrussCountVariableRadii(radiusArray);
        int[][] ConnArray = new int[radiusArray.length][2];
        int currentTrussCount = 0;
        for (int index = 0; index < radiusArray.length; index++){
            double decision = radiusArray[index];
            if (decision > lowestRadiusFactor*radiusUpperBounds[index]){
                ConnArray[currentTrussCount] = completeConnectivityArray[index];
                currentTrussCount += 1;
            }
        }
        int[][] designConnArray;
        designConnArray = sliceFullConnectivityArrayVariableRadii(ConnArray, currentTrussCount);
        return designConnArray;
    }

    public int[][] getFullCAFromSolutionVariableRadii(){
        return designConnectivityArray;
    }

    public double[] getFullRadiusArray(){
        return radiusArray;
    }

    private int getTrussCountVariableRadii(double[] radiusArray){
        int numVars = radiusArray.length;
        int TrussCount = 0;
        for (int index = 0; index < numVars; index++){
            double radiusValue =  radiusArray[index];
            if (radiusValue > lowestRadiusFactor*radiusUpperBounds[index]){
                TrussCount += 1;
            }
        }
        return TrussCount;
    }

    private int[][] sliceFullConnectivityArrayVariableRadii (int[][] fullConnectivityArray, int trussCount) {
        int[][] connectivityArray = new int[trussCount][2];
        for (int i = 0; i < trussCount; i++) {
            for (int j = 0; j < 2; j++) {
                connectivityArray[i][j] = fullConnectivityArray[i][j];
            }
        }
        return connectivityArray;
    }

    public int getTrussCountFromSolutionVariableRadii () {
        return numberOfTrusses;
    }

    public VariableRadiiRepeatableArchitecture getVariableRadiiArchitectureFromRadiusArray (double[] fullRadiusArray) {
        int numVariables = findNumberOfVariablesFromSidenum();
        double[] designRepeatableRadiusArray = getRepeatableRadiusArrayFromCompleteArray(fullRadiusArray);
        Solution architecture = new Solution(numVariables,2);
        for (int i = 0; i < designRepeatableRadiusArray.length; i++) {
            RealVariable var = new RealVariable(radiusLowerBounds[i],radiusUpperBounds[i]);
            EncodingUtils.setReal(var,designRepeatableRadiusArray[i]);
            architecture.setVariable(i,var);
        }
        return new VariableRadiiRepeatableArchitecture(architecture,sidenum,radiusLowerBounds,radiusUpperBounds,lowestRadiusFactor);
    }

    private boolean[] identifyRepeatedEdgeMembers () {
        int[][] completeConnArray = getCompleteConnectivityArrayFromSidenum();
        int[] topNodes = getTopEdgeNodes();
        int memberCount = 0;
        boolean[] isReapeatedEdgeMember = new boolean[completeConnArray.length];
        for (int i = 0; i < ((sidenum*sidenum)-1); i++) {
            boolean rightEdge = false;
            boolean topEdge = false;
            int node = i+1;
            for (int j = node+1; j < ((sidenum*sidenum)+1); j++) {
                if (node > ((sidenum*sidenum)-sidenum)) { // identifying right edge members
                    if (j > ((sidenum*sidenum)-sidenum)) {
                        isReapeatedEdgeMember[memberCount] = true;
                        rightEdge = true;
                    }
                }
                if (IntStream.of(topNodes).anyMatch(x -> x == node)) { // identifying top edge members
                    int finalJ = j;
                    if (IntStream.of(topNodes).anyMatch(x -> x == finalJ)) {
                        isReapeatedEdgeMember[memberCount] = true;
                        topEdge = true;
                    }
                }
                if (!rightEdge && !topEdge) {
                    isReapeatedEdgeMember[memberCount] = false;
                }
                memberCount += 1;
            }
        }
        return isReapeatedEdgeMember;
    }

    public double[] getRepeatableRadiusArrayFromCompleteArray (double[] completeRadiusArray) {
        ArrayList<Double> repeatableRadiusArrayList = new ArrayList<>();
        boolean[] isRepeated = identifyRepeatedEdgeMembers();
        for (int i = 0; i < completeRadiusArray.length; i++) {
            if (!isRepeated[i]) {
                repeatableRadiusArrayList.add(completeRadiusArray[i]);
            }
        }
        return repeatableRadiusArrayList.stream().mapToDouble(x->x).toArray();
    }

    private int[][] getRepeatableConnectivityArrayFromSidenum() {
        int memberCount = 0;
        int numberOfVariables = findNumberOfVariablesFromSidenum();
        int[] topNodes = getTopEdgeNodes();
        int[][] repeatableConnectivityArray = new int[numberOfVariables][2];
        for (int i = 0; i < ((sidenum*sidenum)-1); i++) {
            int node = i+1;
            for (int j = node+1; j < ((sidenum*sidenum)+1); j++) {
                if (node > ((sidenum*sidenum)-sidenum)) { // identifying right edge members
                    if (j > ((sidenum * sidenum) - sidenum)) {
                        continue;
                    }
                }
                if (IntStream.of(topNodes).anyMatch(x -> x == node)) { // identifying top edge members
                    int finalJ = j;
                    if (IntStream.of(topNodes).anyMatch(x -> x == finalJ)) {
                        continue;
                    }
                }
                repeatableConnectivityArray[memberCount][0] = i+1;
                repeatableConnectivityArray[memberCount][1] = j;
                memberCount += 1;
            }
        }
        return repeatableConnectivityArray;
    }

    public double[] getCompleteRadiusArrayFromSolutionNxN (Solution solution) {
        double[] repeatableRadiusArray = getRadiusArray(solution);
        //int totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sidenum*sidenum))/(CombinatoricsUtils.factorial((int) ((sidenum*sidenum) - 2)) * CombinatoricsUtils.factorial(2)));
        int[] topNodes = getTopEdgeNodes();
        int[][] repeatableConnArray = getRepeatableConnectivityArrayFromSidenum();
        ArrayList<Double> completeRadiusArrayList = new ArrayList<>();
        int memberIndex = 0;
        for (int i = 0; i < ((sidenum*sidenum)-1); i++) {
            int node = i+1;
            for (int j = node+1; j < ((sidenum*sidenum)+1); j++) {
                boolean rightEdge = false;
                boolean topEdge = false;
                if (node > ((sidenum*sidenum)-sidenum)) { // identifying right edge members
                    if (j > ((sidenum*sidenum)-sidenum)) {
                        int[] repeatedMember = new int[]{(int) (node - ((sidenum-1) * sidenum)), (int) (j - ((sidenum-1) * sidenum))}; // corresponding left edge member
                        int repeatedMemberIndex = getMemberIndex(repeatableConnArray, repeatedMember);
                        completeRadiusArrayList.add(repeatableRadiusArray[repeatedMemberIndex]);
                        rightEdge = true;
                    }
                }
                if (IntStream.of(topNodes).anyMatch(x -> x == node)) { // identifying top edge members
                    int finalJ = j;
                    if (IntStream.of(topNodes).anyMatch(x -> x == finalJ)) {
                        int[] repeatedMember = new int[]{(int) (node - (sidenum-1)), (int) (j - (sidenum-1))}; // corresponding bottom edge member
                        int repeatedMemberIndex = getMemberIndex(repeatableConnArray, repeatedMember);
                        completeRadiusArrayList.add(repeatableRadiusArray[repeatedMemberIndex]);
                        topEdge = true;
                    }
                }
                if (!rightEdge && !topEdge) {
                    completeRadiusArrayList.add(repeatableRadiusArray[memberIndex]);
                    memberIndex += 1;
                }
            }
        }
        return completeRadiusArrayList.stream().mapToDouble(x->x).toArray();
    }

    private int findNumberOfVariablesFromSidenum () {
        int totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sidenum*sidenum))/(CombinatoricsUtils.factorial((int) ((sidenum*sidenum) - 2)) * CombinatoricsUtils.factorial(2)));
        int numberOfRepeatableMembers = (int) (2 * (CombinatoricsUtils.factorial((int) sidenum)/(CombinatoricsUtils.factorial((int) (sidenum - 2)) * CombinatoricsUtils.factorial(2))));
        return (totalNumberOfMembers - numberOfRepeatableMembers);
    }

    private int getMemberIndex (int[][] designConnArray, int[] member) {
        int index = 0;
        for (int i = 0; i < designConnArray.length; i++) {
            if (designConnArray[i][0] == member[0]) {
                if (designConnArray[i][1] == member[1]) {
                    index = i;
                    break;
                }
            }
        }
        return index;
    }

    private int[] getTopEdgeNodes () {
        ArrayList<Integer> topNodes = new ArrayList<Integer>();
        boolean reachedRightEdge = false;
        int node = (int) sidenum;
        while (!reachedRightEdge) {
            if (node > (sidenum*sidenum)) {
                reachedRightEdge = true;
            } else {
                topNodes.add(node);
                node += sidenum;
            }
        }
        return topNodes.stream().mapToInt(i->i).toArray();
    }

}

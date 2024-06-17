package seakers.trussaos.architecture;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.RealVariable;

import java.util.ArrayList;
import java.util.stream.IntStream;


/**
 * Design class for the integer encoded design vector (variable radii) in the 2D NxN nodal grid case
 *
 * @author roshan94
 */

public class IntegerRepeatableArchitecture extends Solution{
    private static final long serialVersionUID = -3246138689443542872L;
    private int numberOfTrusses;
    private double[] radiusArray;
    private double[] radiiChoices;
    private double[] YoungsModulii;
    private double sidenum;
    private double[][] completeConnectivityArray;
    private double[][] designConnectivityArray;
    private final int numHeurObjectives;
    private final int numHeurConstraints;
    private final int YoungsModulusChoice;

    public IntegerRepeatableArchitecture(Solution solution, double sidenum, int numHeuristicObjectives, int numHeuristicConstraints, double[] radiiChoices, double[] YoungsModulii) {
        super(solution);
        this.sidenum = sidenum;
        this.numHeurConstraints = numHeuristicConstraints;
        this.numHeurObjectives = numHeuristicObjectives;

        // Generate Complete Connectivity Array from sidenum
        completeConnectivityArray = getCompleteConnectivityArrayFromSidenum();

        // Extract the radius array from the solution
        this.radiiChoices = radiiChoices;
        this.YoungsModulii = YoungsModulii;
        radiusArray = getCompleteIntegerRadiusArrayFromSolution(solution);

        // Obtain design connectivity array
        designConnectivityArray = getFullConnectivityArrayIntegerRadii();
        numberOfTrusses = getTrussCountIntegerRadii(radiusArray);

        int numTrussVariables = findNumberOfTrussVariablesFromSidenum();
        YoungsModulusChoice = EncodingUtils.getInt(solution.getVariable(numTrussVariables));
    }

    public int getYoungsModulusChoice() {
        return YoungsModulusChoice;
    }

    private double[] getRadiusArray(Solution soln) {
        int numVars = soln.getNumberOfVariables();
        double[] radii = new double[numVars-1];
        for (int index = 0; index < numVars-1; index++){  // all except the last decision variable (which determines the material)
            double decision =  radiiChoices[(int) EncodingUtils.getReal(soln.getVariable(index))];
            radii[index] = decision;
        }
        return radii;
    }

    public double[] getCompleteIntegerRadiusArrayFromSolution(Solution soln){
        return getCompleteIntegerRadiusArrayFromSolutionNxN(soln);
    }

    public double[][] getCompleteConnectivityArrayFromSidenum(){
        int memberCount = 0;
        int totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sidenum*sidenum))/(CombinatoricsUtils.factorial((int) ((sidenum*sidenum) - 2)) * CombinatoricsUtils.factorial(2)));
        double[][] completeConnectivityArray = new double[totalNumberOfMembers][2];
        for (int i = 0; i < ((sidenum*sidenum)-1); i++) {
            for (int j = i+1; j < (sidenum*sidenum); j++) {
                completeConnectivityArray[memberCount][0] = i+1;
                completeConnectivityArray[memberCount][1] = j+1;
                memberCount += 1;
            }
        }
        return completeConnectivityArray;
    }

    public double[][] getFullConnectivityArrayIntegerRadii() {
        //int TrussCount = getTrussCountVariableRadii(radiusArray);
        double[][] ConnArray = new double[radiusArray.length][2];
        int currentTrussCount = 0;
        for (int index = 0; index < radiusArray.length; index++){
            double decision = radiusArray[index];
            if (decision != 0D){
                ConnArray[currentTrussCount] = completeConnectivityArray[index];
                currentTrussCount += 1;
            }
        }
        double[][] designConnArray;
        designConnArray = sliceFullConnectivityArrayVariableRadii(ConnArray, currentTrussCount);
        return designConnArray;
    }

    public double[][] getFullCAFromSolutionIntegerRadii(){
        return designConnectivityArray;
    }

    public double[] getFullRadiusArray(){
        return radiusArray;
    }

    private int getTrussCountIntegerRadii(double[] radiusArray){
        int numVars = radiusArray.length;
        int TrussCount = 0;
        for (int index = 0; index < numVars; index++){
            double radiusValue =  radiusArray[index];
            if (radiusValue != 0D){
                TrussCount += 1;
            }
        }
        return TrussCount;
    }

    private double[][] sliceFullConnectivityArrayVariableRadii (double[][] fullConnectivityArray, int trussCount) {
        double[][] connectivityArray = new double[trussCount][2];
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

    public IntegerRepeatableArchitecture getIntegerArchitectureFromRadiusArray (double[] fullRadiusChoiceArray, double choiceOfYoungsModulus) {
        int numTrussVariables = findNumberOfTrussVariablesFromSidenum();
        double[] designRepeatableRadiusChoiceArray = getRepeatableIntegerRadiusArrayFromCompleteArray(fullRadiusChoiceArray);
        Solution architecture = new Solution(numTrussVariables+1,2);
        for (int i = 0; i < designRepeatableRadiusChoiceArray.length; i++) {
            RealVariable var = new RealVariable(designRepeatableRadiusChoiceArray[i],0,radiusArray.length);
            architecture.setVariable(i,var);
        }
        RealVariable modulusVar = new RealVariable(choiceOfYoungsModulus,0,YoungsModulii.length);
        architecture.setVariable(numTrussVariables+1,modulusVar);
        return new IntegerRepeatableArchitecture(architecture,sidenum,numHeurObjectives,numHeurConstraints,radiiChoices,YoungsModulii);
    }

    private boolean[] identifyRepeatedEdgeMembers () {
        double[][] completeConnArray = getCompleteConnectivityArrayFromSidenum();
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

    public double[] getRepeatableIntegerRadiusArrayFromCompleteArray (double[] completeRadiusArray) {
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
        int numberOfVariables = findNumberOfTrussVariablesFromSidenum();
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

    public double[] getCompleteIntegerRadiusArrayFromSolutionNxN (Solution solution) { // ACTUAL RADII VALUES (NOT INDICES)!
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

    private int findNumberOfTrussVariablesFromSidenum () {
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

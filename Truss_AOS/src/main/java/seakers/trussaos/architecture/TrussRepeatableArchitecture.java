package seakers.trussaos.architecture;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variable;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.EncodingUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Design class for the binary design vector in the 2D NxN nodal grid case
 *
 * @author roshan94
 */

public class TrussRepeatableArchitecture extends Solution{
    private static final long serialVersionUID = -3246610489443538032L;

    private int NumberOfTrusses;

    private boolean[] CompleteBooleanDesign;

    //private int[][] FullConnectivityArray = {{1,2}, {1,3}, {1,4}, {1,5}, {1,6}, {1,7}, {1,8}, {1,9},
                                            //{2,3}, {2,4}, {2,5}, {2,6}, {2,7}, {2,8}, {2,9},
                                            //{3,4}, {3,5}, {3,6}, {3,7}, {3,8}, {3,9},
                                            //{4,5}, {4,6}, {4,7}, {4,8}, {4,9},
                                            //{5,6}, {5,7}, {5,8}, {5,9},
                                            //{6,7}, {6,8}, {6,9},
                                            //{7,8}, {7,9},
                                            //{8,9}};

    private double[][] CompleteConnectivityArray;

    private double[][] DesignConnectivityArray;

    private final double sidenum;

    private final int numHeurObjectives;

    private final int numHeurConstraints;

    public TrussRepeatableArchitecture(Solution solution, double sidenum, int numHeuristicObjectives, int numHeuristicConstraints) {
        super(solution);
        this.sidenum = sidenum;
        this.numHeurConstraints = numHeuristicConstraints;
        this.numHeurObjectives = numHeuristicObjectives;

        //// Extract the design as a boolean array
        boolean[] BooleanDesign = getBooleanDesignArray(solution);
        this.CompleteConnectivityArray = getCompleteConnectivityArrayFromSidenum();
        //design = EncodingUtils.getBinary(solution.getVariable());

        //// Obtain corresponding Connectivity Array
        DesignConnectivityArray= ConvertToFullConnectivityArray(BooleanDesign);

        //// Compute number of trusses in the design
        CompleteBooleanDesign = getCompleteBooleanDesignFromRepeatableDesign(BooleanDesign);
        NumberOfTrusses = getTrussCount(CompleteBooleanDesign);

    }

    private double[][] ConvertToFullConnectivityArray(boolean[] CurrentDesign){
        boolean[] FullDesign = getCompleteBooleanDesignFromRepeatableDesign(CurrentDesign);
        int TrussCount = 0;
        double[][] ConnArray = new double[FullDesign.length][2];
        //design = EncodingUtils.getBinary(solution.getVariable());
        for (int index = 0; index < FullDesign.length; index++){
            boolean decision = FullDesign[index];
            if (decision){
                ConnArray[TrussCount] = CompleteConnectivityArray[index];
                TrussCount += 1;
            }
        }
        double[][] designConnArray;
        designConnArray = sliceFullConnectivityArray(ConnArray, TrussCount);
        return designConnArray;
    }

    private int getTrussCount(boolean[] CompleteDesign){
        int numVars = CompleteDesign.length;
        int TrussCount = 0;
        for (int index = 0; index < numVars; index++){
            boolean decision =  CompleteDesign[index];
            if (decision){
                TrussCount += 1;
            }
        }
        return TrussCount;
    }

    private boolean[] getBooleanDesignArray (Solution soln){
        int numVars = soln.getNumberOfVariables();
        boolean[] design = new boolean[numVars];
        for (int index = 0; index < numVars; index++){
            boolean decision =  EncodingUtils.getBoolean(soln.getVariable(index));
            design[index] = decision;
        }
        return design;
    }

    public boolean[] getBooleanDesignFromSolution (Solution solution) {
        return getCompleteBooleanArrayFromSolutionNxN(solution);
    }

    public double[][] getConnectivityArrayFromSolution (Solution solution) {
        //// Extract the design as a boolean array
        boolean[] BooleanDesign = getBooleanDesignArray(solution);

        //// Obtain corresponding Connectivity Array
        return ConvertToFullConnectivityArray(BooleanDesign);
    }

    private double[][] sliceFullConnectivityArray (double[][] fullConnectivityArray, int trussCount) {
        double[][] connectivityArray = new double[trussCount][2];
        for (int i = 0; i < trussCount; i++) {
            for (int j = 0; j < 2; j++) {
                connectivityArray[i][j] = fullConnectivityArray[i][j];
            }
        }
        return connectivityArray;
    }

    public int getTrussCountFromSolution () {
        return NumberOfTrusses;
    }

    public TrussRepeatableArchitecture getArchitectureFromConnectivityArray (double[][] connArray, boolean arteryProblem) {
        boolean[] designFullBooleanArray = new boolean[CompleteConnectivityArray.length];
        boolean contains = false;
        int[] designFirstNodes = new int[connArray.length];
        int[] designSecondNodes = new int[connArray.length];
        for (int i = 0; i < connArray.length; i++) {
            designFirstNodes[i] = (int) connArray[i][0];
            designSecondNodes[i] = (int) connArray[i][1];
        }
        for (int i = 0; i < CompleteConnectivityArray.length; i++) {
            int firstNode = (int) CompleteConnectivityArray[i][0];
            int secondNode = (int) CompleteConnectivityArray[i][1];
            for (int j = 0; j < connArray.length; j++) {
                if (designFirstNodes[j] == firstNode) {
                    if (designSecondNodes[j] == secondNode) {
                        contains = true;
                        break;
                    }
                }
            }
            designFullBooleanArray[i] = contains;
            contains = false;
        }
        boolean[] designRepeatableBooleanArray = getRepeatableBooleanArrayFromCompleteArray(designFullBooleanArray);
        Solution architecture;
        if (arteryProblem) {
            architecture = new Solution(findNumberOfVariablesFromSidenum(),2+numHeurObjectives, 2+numHeurConstraints);
        } else {
            architecture = new Solution(findNumberOfVariablesFromSidenum(),2+numHeurObjectives, 3+numHeurConstraints);
        }
        for (int i = 0; i < designRepeatableBooleanArray.length; i++) {
            BinaryVariable var = new BinaryVariable(1);
            EncodingUtils.setBoolean(var,designRepeatableBooleanArray[i]);
            architecture.setVariable(i,var);
        }
        return new TrussRepeatableArchitecture(architecture,sidenum,numHeurObjectives,numHeurConstraints);
    }

    private double[][] getCompleteConnectivityArrayFromSidenum(){
        int memberCount = 0;
        //int[] nodesArray = IntStream.range(1,sidenum*sidenum).toArray();
        int totalNumberOfMembers;
        if (sidenum >= 5) {
            int sidenumSquared = (int) (sidenum*sidenum);
            totalNumberOfMembers =  sidenumSquared * (sidenumSquared - 1)/2;
        }
        else {
            totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sidenum*sidenum))/(CombinatoricsUtils.factorial((int) ((sidenum*sidenum) - 2)) * CombinatoricsUtils.factorial(2)));
        }
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

    private int findNumberOfVariablesFromSidenum () {
        int totalNumberOfMembers;
        if (sidenum >= 5) {
            int sidenumSquared = (int) (sidenum*sidenum);
            totalNumberOfMembers =  sidenumSquared * (sidenumSquared - 1)/2;
        }
        else {
            totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sidenum*sidenum))/(CombinatoricsUtils.factorial((int) ((sidenum*sidenum) - 2)) * CombinatoricsUtils.factorial(2)));
        }
        int numberOfRepeatableMembers = (int) (2 * (CombinatoricsUtils.factorial((int) sidenum)/(CombinatoricsUtils.factorial((int) (sidenum - 2)) * CombinatoricsUtils.factorial(2))));
        return (totalNumberOfMembers - numberOfRepeatableMembers);
    }

    public boolean[] getRepeatableBooleanArrayFromCompleteArray (boolean[] completeBooleanArray) {
        ArrayList<Boolean> repeatableBooleanArray = new ArrayList<>();
        boolean[] isRepeated = identifyRepeatedEdgeMembers();
        for (int i = 0; i < completeBooleanArray.length; i++) {
            if (!isRepeated[i]) {
                repeatableBooleanArray.add(completeBooleanArray[i]);
            }
        }
        return convertBooleanArrayListToArray(repeatableBooleanArray);
    }

    private boolean[] identifyRepeatedEdgeMembers () {
        double[][] completeConnArray = getCompleteConnectivityArrayFromSidenum();
        int[] topNodes = getTopEdgeNodes();
        int memberCount = 0;
        boolean[] isReapeatedEdgeMember = new boolean[completeConnArray.length];
        for (int i = 0; i < ((sidenum*sidenum)-1); i++) {
            boolean rightEdge = false;
            boolean topEgde = false;
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
                        topEgde = true;
                    }
                }
                if (!rightEdge && !topEgde) {
                    isReapeatedEdgeMember[memberCount] = false;
                }
                memberCount += 1;
            }
        }
        return isReapeatedEdgeMember;
    }

    private double[][] getRepeatedConnectivityArrayFromSidenum () {
        int memberCount = 0;
        int numberOfVariables = findNumberOfVariablesFromSidenum();
        int[] topNodes = getTopEdgeNodes();
        double[][] repeatableConnectivityArray = new double[numberOfVariables][2];
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

    private boolean[] getCompleteBooleanDesignFromRepeatableDesign (boolean[] repeatableArray) {
        //int totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sidenum*sidenum))/(CombinatoricsUtils.factorial((int) ((sidenum*sidenum) - 2)) * CombinatoricsUtils.factorial(2)));
        int[] topNodes = getTopEdgeNodes();
        double[][] repeatableConnArray = getRepeatedConnectivityArrayFromSidenum();
        ArrayList<Boolean> completeBooleanArrayList = new ArrayList<>();
        int memberIndex = 0;
        for (int i = 0; i < ((sidenum*sidenum)-1); i++) {
            int node = i+1;
            for (int j = node+1; j < ((sidenum*sidenum)+1); j++) {
                boolean rightEdge = false;
                boolean topEdge = false;
                if (node > ((sidenum*sidenum)-sidenum)) { // identifying right edge members
                    if (j > ((sidenum*sidenum)-sidenum)) {
                        rightEdge = true;
                        int[] repeatedMember = new int[]{(int) (node - ((sidenum-1) * sidenum)), (int) (j - ((sidenum-1) * sidenum))}; // corresponding left edge member
                        int repeatedMemberIndex = getMemberIndex(repeatableConnArray, repeatedMember);
                        completeBooleanArrayList.add(repeatableArray[repeatedMemberIndex]);
                    }
                }
                if (IntStream.of(topNodes).anyMatch(x -> x == node)) { // identifying top edge members
                    int finalJ = j;
                    if (IntStream.of(topNodes).anyMatch(x -> x == finalJ)) {
                        topEdge = true;
                        int[] repeatedMember = new int[]{(int) (node - (sidenum-1)), (int) (j - (sidenum-1))}; // corresponding bottom edge member
                        int repeatedMemberIndex = getMemberIndex(repeatableConnArray, repeatedMember);
                        completeBooleanArrayList.add(repeatableArray[repeatedMemberIndex]);
                    }
                }
                if (!rightEdge && !topEdge) {
                    completeBooleanArrayList.add(repeatableArray[memberIndex]);
                    memberIndex += 1;
                }
            }
        }
        return convertBooleanArrayListToArray(completeBooleanArrayList);
    }

    public boolean[] getCompleteBooleanArrayFromSolutionNxN (Solution solution) {
        boolean[] repeatableBooleanArray = getBooleanDesignArray(solution);
        return getCompleteBooleanDesignFromRepeatableDesign(repeatableBooleanArray);
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

    private int getMemberIndex (double[][] connectivityArray, int[] member) {
        int index = 0;
        for (int i = 0; i < connectivityArray.length; i++) {
            if (connectivityArray[i][0] == member[0]) {
                if (connectivityArray[i][1] == member[1]) {
                    index = i;
                    break;
                }
            }
        }
        return index;
    }

    private boolean[] convertBooleanArrayListToArray (ArrayList<Boolean> booleanArrayList) {
        boolean[] booleanArray = new boolean[booleanArrayList.size()];
        for (int i = 0; i < booleanArrayList.size(); i++) {
            booleanArray[i] = booleanArrayList.get(i);
        }
        return booleanArray;
    }



}

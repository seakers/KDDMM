package seakers.trussaos.initialization;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
//import org.moeaframework.core.Variable;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.EncodingUtils;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;
import seakers.trussaos.operators.constantradii.ImproveOrientation;

//import java.lang.reflect.Array;
import java.util.ArrayList;
//import java.util.HashSet;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * Generates a biased random initialized population for the Truss GA problem wherein equal number of designs with 6, 12.
 * 18 and 24 members in the 32 dimensional design vector are generated. In other words, 4 groups of designs are randomly
 * generated.
 *
 * @author roshan94
 */

public class BiasedInitialization implements Initialization {
    private final Problem problem;
    private final int populationSize;
    private final boolean lowMemberBiasing;
    private final boolean partialCollapsibilityEnforced;
    private final boolean nodalPropertiesEnforced;
    private final boolean orientationEnforced;
    private final boolean intersectionEnforced;
    private final double[][] globalNodePositions;
    private final double targetStiffnessRatio;
    private final double sideNodeNumber;
    private final int numHeurObjectives;
    private final int numHeurConstraints;
    private static double orientaionAcceptanceMargin = 10;

    public BiasedInitialization(Problem problem, int populationSize, boolean lowMemberBiasing, boolean partialCollapsibilityEnforced, boolean nodalPropertiesEnforced, boolean orientationEnforced, boolean intersectionEnforced, double[][] globalNodePositions, double targetStiffnessRatio, int sideNodeNumber, int numHeuristicObjectives, int numHeuristicConstraints) {
        this.problem = problem;
        this.populationSize = populationSize;
        this.lowMemberBiasing = lowMemberBiasing;
        this.partialCollapsibilityEnforced = partialCollapsibilityEnforced;
        this.nodalPropertiesEnforced = nodalPropertiesEnforced;
        this.orientationEnforced = orientationEnforced;
        this.intersectionEnforced = intersectionEnforced; // Currently only used in conjunction with some cases (eg. with orientation)
        this.globalNodePositions = globalNodePositions;
        this.targetStiffnessRatio = targetStiffnessRatio;
        this.numHeurObjectives = numHeuristicObjectives;
        this.numHeurConstraints = numHeuristicConstraints;
        this.sideNodeNumber = sideNodeNumber;
    }

    @Override
    public Solution[] initialize() {
        int solutionsPerGroup = (int) Math.round(populationSize/4);

        double normalProb = 0.5;
        if (lowMemberBiasing) {
            normalProb = 0.0893;
        }

        Solution[] initialPopulation = new Solution[this.populationSize];
        int solutionCount = 0;

        double[] numberOfMembers = new double[4]; // four clusters of solutions to be initialized
        if (nodalPropertiesEnforced && !intersectionEnforced) {
            int numberOfVariables = findNumberOfVariablesFromSidenum();
            numberOfMembers = new double[]{numberOfVariables/2D, numberOfVariables*7D/12D, numberOfVariables*2D/3D, numberOfVariables*3D/4D};
        }
        //if (intersectionEnforced && !nodalPropertiesEnforced) {
            //numberOfMembers = new double[]{6, 9, 12, 15};
        //}
        if (intersectionEnforced && nodalPropertiesEnforced) {
            int numberOfVariables = findNumberOfVariablesFromSidenum();
            numberOfMembers = new double[]{numberOfVariables/4D, numberOfVariables*7D/24D, numberOfVariables*1D/3D, numberOfVariables*3D/8D};
        }

        Random rand = new Random();
        boolean val;

        if (!partialCollapsibilityEnforced && !nodalPropertiesEnforced && !orientationEnforced && !intersectionEnforced) {
            for (int i = 0;i < 4;i++) {
                for (int j = 0;j < solutionsPerGroup;j++) {
                    Solution solution = this.problem.newSolution();
                    for (int k = 0;k < solution.getNumberOfVariables();++k) {
                        val = rand.nextFloat() < normalProb;
                        BinaryVariable var = new BinaryVariable(1);
                        EncodingUtils.setBoolean(var,val);
                        solution.setVariable(k,var);
                    }
                    initialPopulation[solutionCount] = solution;
                    solutionCount++;
                }
            }
        }

        if (nodalPropertiesEnforced && !partialCollapsibilityEnforced && !orientationEnforced) {
            int numberOfVariables = findNumberOfVariablesFromSidenum();
            for (int i = 0;i < 4;i++) {
                double TrueProb = numberOfMembers[i]/numberOfVariables;
                for (int j = 0;j < solutionsPerGroup;j++) {
                    Solution solution = this.problem.newSolution();
                    for (int k = 0;k < solution.getNumberOfVariables();++k) {
                        val = rand.nextFloat() < TrueProb;
                        BinaryVariable var = new BinaryVariable(1);
                        EncodingUtils.setBoolean(var,val);
                        solution.setVariable(k,var);
                    }
                    initialPopulation[solutionCount] = solution;
                    solutionCount++;
                }
            }
        }

        if (partialCollapsibilityEnforced && !nodalPropertiesEnforced && !orientationEnforced) {
            double diagProb = 0.75;
            boolean[] booleanDiagonals = identifyDiagonalMembers();
            boolean[] booleanDiagonalsRepeatable = getRepeatableBooleanArrayFromCompleteArray(booleanDiagonals);
            for (int i = 0; i < populationSize; i++) {
                Solution solution = this.problem.newSolution();
                for (int k = 0; k < solution.getNumberOfVariables(); ++k) {
                    if (booleanDiagonalsRepeatable[k]) {
                        val = rand.nextFloat() < diagProb;
                    } else {
                        //val = rand.nextFloat() < 0.5;
                        val = rand.nextFloat() < normalProb;
                    }
                    BinaryVariable var = new BinaryVariable(1);
                    EncodingUtils.setBoolean(var,val);
                    solution.setVariable(k,var);
                }
                initialPopulation[solutionCount] = solution;
                solutionCount++;
            }
        }

        if (partialCollapsibilityEnforced && nodalPropertiesEnforced && !orientationEnforced) {
            double diagProb = 0.75;
            boolean[] booleanDiagonals = identifyDiagonalMembers();
            boolean[] booleanDiagonalsRepeatable = getRepeatableBooleanArrayFromCompleteArray(booleanDiagonals);
            int numberOfVariables = findNumberOfVariablesFromSidenum();
            for (int i = 0;i < 4;i++) {
                double TrueProb = numberOfMembers[i]/numberOfVariables;
                for (int j = 0;j < solutionsPerGroup;j++) {
                    Solution solution = this.problem.newSolution();
                    for (int k = 0;k < solution.getNumberOfVariables();++k) {
                        if (booleanDiagonalsRepeatable[k]) {
                            val = rand.nextFloat() < diagProb;
                        } else {
                            val = rand.nextFloat() < TrueProb;
                        }
                        BinaryVariable var = new BinaryVariable(1);
                        EncodingUtils.setBoolean(var,val);
                        solution.setVariable(k,var);
                    }
                    initialPopulation[solutionCount] = solution;
                    solutionCount++;
                }
            }
        }

        ImproveOrientation improveOrientation;
        double targetOrientation;
        if (!partialCollapsibilityEnforced && !nodalPropertiesEnforced && orientationEnforced && !intersectionEnforced) {
            improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
            targetOrientation = improveOrientation.findTargetOrientation(targetStiffnessRatio);
            int m = 0;
            while (m < populationSize) {
                Solution solution = this.problem.newSolution();
                for (int n = 0; n < solution.getNumberOfVariables(); ++n) {
                    //val = rand.nextFloat() < 0.5;
                    val = rand.nextFloat() < normalProb;
                    BinaryVariable var = new BinaryVariable(1);
                    EncodingUtils.setBoolean(var, val);
                    solution.setVariable(n, var);
                }
                TrussRepeatableArchitecture trussArch = new TrussRepeatableArchitecture(solution, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
                double[][] connectivityArray = trussArch.getConnectivityArrayFromSolution(solution);

                double[] architectureOrientations = improveOrientation.findMemberOrientations(connectivityArray);

                double totalOrientation = 0;
                for (int i = 0; i < architectureOrientations.length; i++) {
                    totalOrientation += architectureOrientations[i];
                }
                double meanOrientation = totalOrientation / architectureOrientations.length;
                if (Math.abs(meanOrientation - targetOrientation) < orientaionAcceptanceMargin) {
                    initialPopulation[m] = solution;
                    m++;
                }
            }
        } else if (partialCollapsibilityEnforced && nodalPropertiesEnforced && orientationEnforced) {
            double diagProb = 0.75;
            int numberOfVariables = findNumberOfVariablesFromSidenum();
            boolean[] booleanDiagonals = identifyDiagonalMembers();
            boolean[] booleanDiagonalsRepeatable = getRepeatableBooleanArrayFromCompleteArray(booleanDiagonals);
            improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
            targetOrientation = improveOrientation.findTargetOrientation(targetStiffnessRatio);
            for (int i = 0;i < 4;i++) {
                double TrueProb = numberOfMembers[i]/numberOfVariables;
                int j = 0;
                while (j < solutionsPerGroup) {
                    Solution solution = this.problem.newSolution();
                    for (int k = 0;k < solution.getNumberOfVariables();++k) {
                        if (booleanDiagonalsRepeatable[k]) {
                            val = rand.nextFloat() < diagProb;
                        } else {
                            val = rand.nextFloat() < TrueProb;
                        }
                        BinaryVariable var = new BinaryVariable(1);
                        EncodingUtils.setBoolean(var,val);
                        solution.setVariable(k,var);
                    }
                    TrussRepeatableArchitecture trussArch = new TrussRepeatableArchitecture(solution, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
                    double[][] connectivityArray = trussArch.getConnectivityArrayFromSolution(solution);

                    double[] architectureOrientations = improveOrientation.findMemberOrientations(connectivityArray);

                    double totalOrientation = 0;
                    for (int s = 0; s < architectureOrientations.length; s++) {
                        totalOrientation += architectureOrientations[s];
                    }
                    double meanOrientation = totalOrientation / architectureOrientations.length;
                    if (Math.abs(meanOrientation - targetOrientation) < orientaionAcceptanceMargin) {
                        initialPopulation[j] = solution;
                        j++;
                    }
                }
            }
        } else if (!partialCollapsibilityEnforced && nodalPropertiesEnforced && orientationEnforced) {
            int numberOfVariables = findNumberOfVariablesFromSidenum();
            improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
            targetOrientation = improveOrientation.findTargetOrientation(targetStiffnessRatio);
            for (int i = 0;i < 4;i++) {
                double TrueProb = numberOfMembers[i]/numberOfVariables;
                int j = 0;
                while (j < solutionsPerGroup) {
                    Solution solution = this.problem.newSolution();
                    for (int k = 0;k < solution.getNumberOfVariables();++k) {
                        val = rand.nextFloat() < TrueProb;
                        BinaryVariable var = new BinaryVariable(1);
                        EncodingUtils.setBoolean(var,val);
                        solution.setVariable(k,var);
                    }
                    TrussRepeatableArchitecture trussArch = new TrussRepeatableArchitecture(solution, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
                    double[][] connectivityArray = trussArch.getConnectivityArrayFromSolution(solution);

                    double[] architectureOrientations = improveOrientation.findMemberOrientations(connectivityArray);

                    double totalOrientation = 0;
                    for (int s = 0; s < architectureOrientations.length; s++) {
                        totalOrientation += architectureOrientations[s];
                    }
                    double meanOrientation = totalOrientation / architectureOrientations.length;
                    if (Math.abs(meanOrientation - targetOrientation) < orientaionAcceptanceMargin) {
                        initialPopulation[j] = solution;
                        j++;
                    }
                }
            }
        } else if (partialCollapsibilityEnforced && !nodalPropertiesEnforced && orientationEnforced && !intersectionEnforced) {
            double diagProb = 0.75;
            boolean[] booleanDiagonals = identifyDiagonalMembers();
            boolean[] booleanDiagonalsRepeatable = getRepeatableBooleanArrayFromCompleteArray(booleanDiagonals);
            improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
            targetOrientation = improveOrientation.findTargetOrientation(targetStiffnessRatio);
            int m = 0;
            while (m < populationSize) {
                Solution solution = this.problem.newSolution();
                for (int n = 0; n < solution.getNumberOfVariables(); ++n) {
                    if (booleanDiagonalsRepeatable[n]) {
                        val = rand.nextFloat() < diagProb;
                    } else {
                        //val = rand.nextFloat() < 0.5;
                        val = rand.nextFloat() < normalProb;
                    }
                    BinaryVariable var = new BinaryVariable(1);
                    EncodingUtils.setBoolean(var, val);
                    solution.setVariable(n, var);
                }
                TrussRepeatableArchitecture trussArch = new TrussRepeatableArchitecture(solution, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
                double[][] connectivityArray = trussArch.getConnectivityArrayFromSolution(solution);

                double[] architectureOrientations = improveOrientation.findMemberOrientations(connectivityArray);

                double totalOrientation = 0;
                for (int i = 0; i < architectureOrientations.length; i++) {
                    totalOrientation += architectureOrientations[i];
                }
                double meanOrientation = totalOrientation / architectureOrientations.length;
                if (Math.abs(meanOrientation - targetOrientation) < orientaionAcceptanceMargin) {
                    initialPopulation[m] = solution;
                    m++;
                }
            }
        }
        if (!partialCollapsibilityEnforced && !nodalPropertiesEnforced && orientationEnforced && intersectionEnforced) {
            int numberOfVariables = findNumberOfVariablesFromSidenum();
            improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
            targetOrientation = improveOrientation.findTargetOrientation(targetStiffnessRatio);
            for (int i = 0;i < 4;i++) {
                double TrueProb = numberOfMembers[i]/numberOfVariables;
                int j = 0;
                while (j < solutionsPerGroup) {
                    Solution solution = this.problem.newSolution();
                    for (int k = 0;k < solution.getNumberOfVariables();++k) {
                        val = rand.nextFloat() < TrueProb;
                        BinaryVariable var = new BinaryVariable(1);
                        EncodingUtils.setBoolean(var,val);
                        solution.setVariable(k,var);
                    }
                    TrussRepeatableArchitecture trussArch = new TrussRepeatableArchitecture(solution, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
                    double[][] connectivityArray = trussArch.getConnectivityArrayFromSolution(solution);

                    double[] architectureOrientations = improveOrientation.findMemberOrientations(connectivityArray);

                    double totalOrientation = 0;
                    for (int s = 0; s < architectureOrientations.length; s++) {
                        totalOrientation += architectureOrientations[s];
                    }
                    double meanOrientation = totalOrientation / architectureOrientations.length;
                    if (Math.abs(meanOrientation - targetOrientation) < orientaionAcceptanceMargin) {
                        initialPopulation[j] = solution;
                        j++;
                    }
                }
            }
        }
        return initialPopulation;
    }

    private boolean[] identifyDiagonalMembers () {
        ArrayList<Boolean> isDiagonalMembers = new ArrayList<Boolean>();

        for (int i = 0; i < ((sideNodeNumber * sideNodeNumber) - 1); i++) {
            int nodeNumber = i+1;
            int closestTopNode = findClosestTopNode(nodeNumber);
            for (int j = nodeNumber+1; j < ((sideNodeNumber * sideNodeNumber)+1); j++) {
                if ((j <= closestTopNode) || ((j - nodeNumber)%sideNodeNumber == 0)) {
                    isDiagonalMembers.add(false);
                } else {
                    isDiagonalMembers.add(true);
                }
            }
        }
        // Create function that converts this to repeatable form and use that
        return convertBooleanArrayListToBooleanArray(isDiagonalMembers);
    }

    private int findClosestTopNode (int node) {
        boolean topReached = false;
        int currentNode = node;
        while (!topReached) {
            if (currentNode%sideNodeNumber == 0) {
                topReached = true;
            } else {
                currentNode += 1;
            }
        }
        return currentNode;
    }

    private boolean[] convertBooleanArrayListToBooleanArray (ArrayList<Boolean> booleanArrayList) {
        boolean[] booleanArray = new boolean[booleanArrayList.size()];
        for (int i = 0; i < booleanArrayList.size(); i++) {
            booleanArray[i] = booleanArrayList.get(i);
        }
        return booleanArray;
    }

    private int findNumberOfVariablesFromSidenum () {
        int totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sideNodeNumber*sideNodeNumber))/(CombinatoricsUtils.factorial((int) ((sideNodeNumber*sideNodeNumber) - 2)) * CombinatoricsUtils.factorial(2)));
        int numberOfRepeatableMembers = (int) (2 * (CombinatoricsUtils.factorial((int) sideNodeNumber)/(CombinatoricsUtils.factorial((int) (sideNodeNumber - 2)) * CombinatoricsUtils.factorial(2))));
        return (totalNumberOfMembers - numberOfRepeatableMembers);
    }

    private boolean[] getRepeatableBooleanArrayFromCompleteArray (boolean[] completeBooleanArray) {
        ArrayList<Boolean> repeatableBooleanArray = new ArrayList<>();
        boolean[] isRepeated = identifyRepeatedEdgeMembers();
        for (int i = 0; i < completeBooleanArray.length; i++) {
            if (!isRepeated[i]) {
                repeatableBooleanArray.add(completeBooleanArray[i]);
            }
        }
        return convertBooleanArrayListToBooleanArray(repeatableBooleanArray);
    }

    private boolean[] identifyRepeatedEdgeMembers () {
        int[][] completeConnArray = getCompleteConnectivityArrayFromSidenum();
        int[] topNodes = getTopEdgeNodes();
        int memberCount = 0;
        boolean[] isReapeatedEdgeMember = new boolean[completeConnArray.length];
        for (int i = 0; i < ((sideNodeNumber*sideNodeNumber)-1); i++) {
            boolean rightEdge = false;
            boolean topEgde = false;
            int node = i+1;
            for (int j = node+1; j < ((sideNodeNumber*sideNodeNumber)+1); j++) {
                if (node > ((sideNodeNumber*sideNodeNumber)-sideNodeNumber)) { // identifying right edge members
                    if (j > ((sideNodeNumber*sideNodeNumber)-sideNodeNumber)) {
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

    private int[][] getCompleteConnectivityArrayFromSidenum(){
        int memberCount = 0;
        //int[] nodesArray = IntStream.range(1,sidenum*sidenum).toArray();
        int totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sideNodeNumber*sideNodeNumber))/(CombinatoricsUtils.factorial((int) ((sideNodeNumber*sideNodeNumber) - 2)) * CombinatoricsUtils.factorial(2)));
        int[][] completeConnectivityArray = new int[totalNumberOfMembers][2];
        for (int i = 0; i < ((sideNodeNumber*sideNodeNumber)-1); i++) {
            for (int j = i+1; j < (sideNodeNumber*sideNodeNumber); j++) {
                completeConnectivityArray[memberCount][0] = i+1;
                completeConnectivityArray[memberCount][1] = j+1;
                memberCount += 1;
            }
        }
        return completeConnectivityArray;
    }

    private int[] getTopEdgeNodes () {
        ArrayList<Integer> topNodes = new ArrayList<Integer>();
        boolean reachedRightEdge = false;
        int node = (int) sideNodeNumber;
        while (!reachedRightEdge) {
            topNodes.add(node);
            if (node > (sideNodeNumber*sideNodeNumber)) {
                reachedRightEdge = true;
            } else {
                node += sideNodeNumber;
            }
        }
        return topNodes.stream().mapToInt(i->i).toArray();
    }
}

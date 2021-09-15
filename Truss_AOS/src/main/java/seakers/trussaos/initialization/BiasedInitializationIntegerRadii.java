package seakers.trussaos.initialization;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
//import org.moeaframework.core.Variable;
import org.moeaframework.core.variable.RealVariable;
import seakers.trussaos.architecture.IntegerRepeatableArchitecture;
import seakers.trussaos.operators.integerradii.ImproveOrientationIntegerRadii;

//import java.lang.reflect.Array;
import java.util.ArrayList;
//import java.util.HashSet;
//import java.util.Random;
import java.util.stream.IntStream;

/**
 * Generates a biased random initialized population for the Truss GA integer radii problem based on the heuristic(s) to be enforces
 *
 * @author roshan94
 */

public class BiasedInitializationIntegerRadii implements Initialization{

    private final Problem problem;
    private final int populationSize;
    private final boolean partialCollapsibilityConstrained;
    private final boolean nodalPropertiesConstrained;
    private final boolean orientationConstrained;
    private final boolean intersectionEnforced;
    private final double[][] globalNodePositions;
    private final double targetStiffnessRatio;
    private final double sideNodeNumber;
    private final double[] radii;
    private final double[] youngsModulii;
    private final int numHeurObjectives;
    private final int numHeurConstraints;
    private static double orientaionAcceptanceMargin = 10;

    public BiasedInitializationIntegerRadii(Problem problem, int populationSize, double[] radii, double[] YoungsModulii, boolean partialCollapsibilityConstrained, boolean nodalPropertiesConstrained, boolean orientationConstrained, boolean intersectionEnforced, double[][] globalNodePositions, double targetStiffnessRatio, int sideNodeNumber, int numHeuristicObjectives, int numHeuristicConstraints) {
        this.problem = problem;
        this.populationSize = populationSize;
        this.partialCollapsibilityConstrained = partialCollapsibilityConstrained;
        this.nodalPropertiesConstrained = nodalPropertiesConstrained;
        this.orientationConstrained = orientationConstrained;
        this.intersectionEnforced = intersectionEnforced; // Currently only used in conjunction with some cases (eg. with orientation)
        this.globalNodePositions = globalNodePositions;
        this.targetStiffnessRatio = targetStiffnessRatio;
        this.sideNodeNumber = sideNodeNumber;
        this.radii = radii;
        this.youngsModulii = YoungsModulii;
        this.numHeurObjectives = numHeuristicObjectives;
        this.numHeurConstraints = numHeuristicConstraints;
    }

    @Override
    public Solution[] initialize() {
        synchronized (PRNG.getRandom()) {
            int solutionsPerGroup = (int) Math.round(populationSize/4);

            Solution[] initialPopulation = new Solution[this.populationSize];
            int solutionCount = 0;

            double[] numberOfMembers = new double[4];
            if (nodalPropertiesConstrained && !intersectionEnforced) {
                int numberOfTrussVariables = findNumberOfTrussVariablesFromSidenum();
                numberOfMembers = new double[]{numberOfTrussVariables/2D, numberOfTrussVariables*7D/12D, numberOfTrussVariables*2D/3D, numberOfTrussVariables*3D/4D};
            }
            if (intersectionEnforced && !nodalPropertiesConstrained) {
                int numberOfTrussVariables = findNumberOfTrussVariablesFromSidenum();
                numberOfMembers = new double[]{numberOfTrussVariables/6D, numberOfTrussVariables/4D, numberOfTrussVariables/3D, numberOfTrussVariables/2D};
            }
            if (intersectionEnforced && nodalPropertiesConstrained) {
                int numberOfTrussVariables = findNumberOfTrussVariablesFromSidenum();
                numberOfMembers = new double[]{numberOfTrussVariables/4D, numberOfTrussVariables*7D/24D, numberOfTrussVariables*1D/3D, numberOfTrussVariables*3D/8D};
            }
            boolean val;

            if ((nodalPropertiesConstrained || intersectionEnforced) && !partialCollapsibilityConstrained && !orientationConstrained) {
                int numberOfTrussVariables = findNumberOfTrussVariablesFromSidenum();
                for (int i = 0;i < 4;i++) {
                    double TrueProb = numberOfMembers[i]/numberOfTrussVariables;
                    for (int j = 0;j < solutionsPerGroup;j++) {
                        Solution solution = this.problem.newSolution();
                        for (int k = 0;k <= numberOfTrussVariables-1;++k) {
                            val = PRNG.nextFloat() < TrueProb;
                            RealVariable newVar = new RealVariable(0,radii.length);
                            if (val) {
                                newVar.setValue(PRNG.nextInt(1,radii.length-1));
                            } else {
                                newVar.setValue(0);
                            }
                            solution.setVariable(k, newVar);
                        }
                        RealVariable newEVar = new RealVariable(0,youngsModulii.length);
                        newEVar.setValue(PRNG.nextInt(youngsModulii.length));
                        solution.setVariable(numberOfTrussVariables, newEVar);
                        initialPopulation[solutionCount] = solution;
                        solutionCount++;
                    }
                }
            }

            if (partialCollapsibilityConstrained && !nodalPropertiesConstrained && !orientationConstrained) {
                int numberOfTrussVariables = findNumberOfTrussVariablesFromSidenum();
                double diagProb = 0.75;
                boolean[] booleanDiagonals = identifyDiagonalMembers();
                boolean[] booleanDiagonalsRepeatable = getRepeatableBooleanArrayFromCompleteArray(booleanDiagonals);
                for (int i = 0; i < populationSize; i++) {
                    Solution solution = this.problem.newSolution();
                    for (int k = 0; k <= numberOfTrussVariables-1; ++k) {
                        RealVariable newVar = new RealVariable(0,radii.length);
                        if (booleanDiagonalsRepeatable[k]) {
                            val = PRNG.nextFloat() < diagProb;
                            if (val) {
                                newVar.setValue(PRNG.nextInt(1,radii.length-1));
                            } else {
                                newVar.setValue(0);
                            }
                        } else {
                            val = PRNG.nextFloat() < 0.5;
                            if (val) {
                                newVar.setValue(PRNG.nextInt(1,radii.length-1));
                            } else {
                                newVar.setValue(0);
                            }
                        }
                        solution.setVariable(k, newVar);
                    }
                    RealVariable newEVar = new RealVariable(0,youngsModulii.length);
                    newEVar.setValue(PRNG.nextInt(youngsModulii.length));
                    solution.setVariable(numberOfTrussVariables, newEVar);
                    initialPopulation[solutionCount] = solution;
                    solutionCount++;
                }
            }

            if (partialCollapsibilityConstrained && nodalPropertiesConstrained && !orientationConstrained) {
                double diagProb = 0.75;
                boolean[] booleanDiagonals = identifyDiagonalMembers();
                boolean[] booleanDiagonalsRepeatable = getRepeatableBooleanArrayFromCompleteArray(booleanDiagonals);
                int numberOfTrussVariables = findNumberOfTrussVariablesFromSidenum();
                for (int i = 0;i < 4;i++) {
                    double TrueProb = numberOfMembers[i]/numberOfTrussVariables;
                    for (int j = 0;j < solutionsPerGroup;j++) {
                        Solution solution = this.problem.newSolution();
                        for (int k = 0;k <= numberOfTrussVariables-1;++k) {
                            RealVariable newVar = new RealVariable(0,radii.length);
                            if (booleanDiagonalsRepeatable[k]) {
                                val = PRNG.nextFloat() < diagProb;
                                if (val) {
                                    newVar.setValue(PRNG.nextInt(1,radii.length-1));
                                } else {
                                    newVar.setValue(0);
                                }
                            } else {
                                val = PRNG.nextFloat() < TrueProb;
                                if (val) {
                                    newVar.setValue(PRNG.nextInt(1,radii.length-1));
                                } else {
                                    newVar.setValue(0);
                                }
                            }
                            solution.setVariable(k, newVar);
                        }
                        RealVariable newEVar = new RealVariable(0,youngsModulii.length);
                        newEVar.setValue(PRNG.nextInt(youngsModulii.length));
                        solution.setVariable(numberOfTrussVariables, newEVar);
                        initialPopulation[solutionCount] = solution;
                        solutionCount++;
                    }
                }
            }

            ImproveOrientationIntegerRadii improveOrientation;
            double targetOrientation;
            if (!partialCollapsibilityConstrained && !nodalPropertiesConstrained && orientationConstrained && !intersectionEnforced) {
                int numberOfTrussVariables = findNumberOfTrussVariablesFromSidenum();
                improveOrientation = new ImproveOrientationIntegerRadii(globalNodePositions, radii, youngsModulii, targetStiffnessRatio, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
                targetOrientation = improveOrientation.findTargetOrientation(targetStiffnessRatio);
                int m = 0;
                while (m < populationSize) {
                    Solution solution = this.problem.newSolution();
                    for (int n = 0; n <= numberOfTrussVariables-1; ++n) {
                        val = PRNG.nextFloat() < 0.5;
                        RealVariable newVar = new RealVariable(0,radii.length);
                        if (val) {
                            newVar.setValue(PRNG.nextInt(1,radii.length-1));
                        } else {
                            newVar.setValue(0);
                        }
                        solution.setVariable(n, newVar);
                    }
                    RealVariable newEVar = new RealVariable(0,youngsModulii.length);
                    newEVar.setValue(PRNG.nextInt(youngsModulii.length));
                    solution.setVariable(numberOfTrussVariables, newEVar);
                    IntegerRepeatableArchitecture trussArch = new IntegerRepeatableArchitecture(solution, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints, radii, youngsModulii);
                    double[][] connectivityArray = trussArch.getFullCAFromSolutionIntegerRadii();

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
            } else if (partialCollapsibilityConstrained && nodalPropertiesConstrained && orientationConstrained) {
                double diagProb = 0.75;
                int numberOfTrussVariables = findNumberOfTrussVariablesFromSidenum();
                boolean[] booleanDiagonals = identifyDiagonalMembers();
                boolean[] booleanDiagonalsRepeatable = getRepeatableBooleanArrayFromCompleteArray(booleanDiagonals);
                improveOrientation = new ImproveOrientationIntegerRadii(globalNodePositions, radii, youngsModulii, targetStiffnessRatio, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
                targetOrientation = improveOrientation.findTargetOrientation(targetStiffnessRatio);
                for (int i = 0;i < 4;i++) {
                    double TrueProb = numberOfMembers[i]/numberOfTrussVariables;
                    int j = 0;
                    while (j < solutionsPerGroup) {
                        Solution solution = this.problem.newSolution();
                        for (int k = 0;k <= numberOfTrussVariables-1;++k) {
                            RealVariable newVar = new RealVariable(0,radii.length);
                            if (booleanDiagonalsRepeatable[k]) {
                                val = PRNG.nextFloat() < diagProb;
                                if (val) {
                                    newVar.setValue(PRNG.nextInt(1,radii.length-1));
                                } else {
                                    newVar.setValue(0);
                                }
                            } else {
                                val = PRNG.nextFloat() < TrueProb;
                                if (val) {
                                    newVar.setValue(PRNG.nextInt(1,radii.length-1));
                                } else {
                                    newVar.setValue(0);
                                }
                            }
                            solution.setVariable(k, newVar);
                        }
                        RealVariable newEVar = new RealVariable(0,youngsModulii.length);
                        newEVar.setValue(PRNG.nextInt(youngsModulii.length));
                        solution.setVariable(numberOfTrussVariables, newEVar);
                        IntegerRepeatableArchitecture trussArch = new IntegerRepeatableArchitecture(solution, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints, radii, youngsModulii);;
                        double[][] connectivityArray = trussArch.getFullCAFromSolutionIntegerRadii();

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
            } else if (!partialCollapsibilityConstrained && nodalPropertiesConstrained && orientationConstrained) {
                int numberOfTrussVariables = findNumberOfTrussVariablesFromSidenum();
                improveOrientation = new ImproveOrientationIntegerRadii(globalNodePositions, radii, youngsModulii, targetStiffnessRatio, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
                targetOrientation = improveOrientation.findTargetOrientation(targetStiffnessRatio);
                for (int i = 0;i < 4;i++) {
                    double TrueProb = numberOfMembers[i]/numberOfTrussVariables;
                    int j = 0;
                    while (j < solutionsPerGroup) {
                        Solution solution = this.problem.newSolution();
                        for (int k = 0;k <= numberOfTrussVariables-1;++k) {
                            val = PRNG.nextFloat() < TrueProb;
                            RealVariable newVar = new RealVariable(0,radii.length);
                            if (val) {
                                newVar.setValue(PRNG.nextInt(1,radii.length-1));
                            } else {
                                newVar.setValue(0);
                            }
                            solution.setVariable(k, newVar);
                        }
                        RealVariable newEVar = new RealVariable(0,youngsModulii.length);
                        newEVar.setValue(PRNG.nextInt(youngsModulii.length));
                        solution.setVariable(numberOfTrussVariables, newEVar);
                        IntegerRepeatableArchitecture trussArch = new IntegerRepeatableArchitecture(solution, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints, radii, youngsModulii);;
                        double[][] connectivityArray = trussArch.getFullCAFromSolutionIntegerRadii();

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
            } else if (partialCollapsibilityConstrained && !nodalPropertiesConstrained && orientationConstrained && !intersectionEnforced) {
                int numberOfTrussVariables = findNumberOfTrussVariablesFromSidenum();
                double diagProb = 0.75;
                boolean[] booleanDiagonals = identifyDiagonalMembers();
                boolean[] booleanDiagonalsRepeatable = getRepeatableBooleanArrayFromCompleteArray(booleanDiagonals);
                improveOrientation = new ImproveOrientationIntegerRadii(globalNodePositions, radii, youngsModulii, targetStiffnessRatio, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
                targetOrientation = improveOrientation.findTargetOrientation(targetStiffnessRatio);
                int m = 0;
                while (m < populationSize) {
                    Solution solution = this.problem.newSolution();
                    for (int n = 0; n <= numberOfTrussVariables-1; ++n) {
                        RealVariable newVar = new RealVariable(0,radii.length);
                        if (booleanDiagonalsRepeatable[n]) {
                            val = PRNG.nextFloat() < diagProb;
                            if (val) {
                                newVar.setValue(PRNG.nextInt(1,radii.length-1));
                            } else {
                                newVar.setValue(0);
                            }
                        } else {
                            val = PRNG.nextFloat() < 0.5;
                            if (val) {
                                newVar.setValue(PRNG.nextDouble(1,radii.length-1));
                            } else {
                                newVar.setValue(0);
                            }
                        }
                        solution.setVariable(n, newVar);
                    }
                    RealVariable newEVar = new RealVariable(0,youngsModulii.length);
                    newEVar.setValue(PRNG.nextInt(youngsModulii.length));
                    solution.setVariable(numberOfTrussVariables, newEVar);
                    IntegerRepeatableArchitecture trussArch = new IntegerRepeatableArchitecture(solution, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints, radii, youngsModulii);;
                    double[][] connectivityArray = trussArch.getFullCAFromSolutionIntegerRadii();

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
            if (!partialCollapsibilityConstrained && !nodalPropertiesConstrained && orientationConstrained && intersectionEnforced) {
                int numberOfTrussVariables = findNumberOfTrussVariablesFromSidenum();
                improveOrientation = new ImproveOrientationIntegerRadii(globalNodePositions, radii, youngsModulii, targetStiffnessRatio, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints);
                targetOrientation = improveOrientation.findTargetOrientation(targetStiffnessRatio);
                for (int i = 0;i < 4;i++) {
                    double TrueProb = numberOfMembers[i]/numberOfTrussVariables;
                    int j = 0;
                    while (j < solutionsPerGroup) {
                        Solution solution = this.problem.newSolution();
                        for (int k = 0;k <= numberOfTrussVariables-1;++k) {
                            val = PRNG.nextFloat() < TrueProb;
                            RealVariable newVar = new RealVariable(0,radii.length);
                            if (val) {
                                newVar.setValue(PRNG.nextInt(1,radii.length-1));
                            } else {
                                newVar.setValue(0);
                            }
                            solution.setVariable(k, newVar);
                        }
                        RealVariable newEVar = new RealVariable(0,youngsModulii.length);
                        newEVar.setValue(PRNG.nextInt(youngsModulii.length));
                        solution.setVariable(numberOfTrussVariables, newEVar);
                        IntegerRepeatableArchitecture trussArch = new IntegerRepeatableArchitecture(solution, (int) sideNodeNumber, numHeurObjectives, numHeurConstraints, radii, youngsModulii);;
                        double[][] connectivityArray = trussArch.getFullCAFromSolutionIntegerRadii();

                        double[] architectureOrientations = improveOrientation.findMemberOrientations(connectivityArray);

                        double totalOrientation = 0;
                        for (double architectureOrientation : architectureOrientations) {
                            totalOrientation += architectureOrientation;
                        }
                        double meanOrientation = totalOrientation / architectureOrientations.length;
                        if (Math.abs(meanOrientation - targetOrientation) < orientaionAcceptanceMargin) {
                            initialPopulation[j] = solution;
                            j++;
                        }
                    }
                }
            }
            // return (Solution[]) initialPopulation.toArray();
            return initialPopulation;
        }
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

    private int findNumberOfTrussVariablesFromSidenum () {
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
            boolean topEdge = false;
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

    private int[][] getCompleteConnectivityArrayFromSidenum(){
        int memberCount = 0;
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
            if (node > (sideNodeNumber*sideNodeNumber)) {
                reachedRightEdge = true;
            } else {
                topNodes.add(node);
                node += sideNodeNumber;
            }
        }
        return topNodes.stream().mapToInt(i->i).toArray();
    }


}


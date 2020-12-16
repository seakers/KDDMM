package seakers.trussaos.initialization;

import org.moeaframework.core.Initialization;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
//import org.moeaframework.core.Variable;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.EncodingUtils;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;
import seakers.trussaos.operators.ImproveOrientation;

//import java.lang.reflect.Array;
import java.util.ArrayList;
//import java.util.HashSet;
import java.util.Random;

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
    private final boolean feasibilityConstrained;
    private final boolean stabilityConstrained;
    private final boolean orientationConstrained;
    private final double[][] globalNodePositions;
    private final double targetStiffnessRatio;
    private final int sideNodeNumber;
    private static double acceptanceMargin = 10;

    public BiasedInitialization(Problem problem, int populationSize, boolean feasibilityConstrained, boolean stabilityConstrained, boolean orientationConstrained, double[][] globalNodePositions, double targetStiffnessRatio, int sideNodeNumber) {
        this.problem = problem;
        this.populationSize = populationSize;
        this.feasibilityConstrained = feasibilityConstrained;
        this.stabilityConstrained = stabilityConstrained;
        this.orientationConstrained = orientationConstrained;
        this.globalNodePositions = globalNodePositions;
        this.targetStiffnessRatio = targetStiffnessRatio;
        this.sideNodeNumber = sideNodeNumber;
    }

    @Override
    public Solution[] initialize() {
        int solutionsPerGroup = (int) Math.round(populationSize/4);

        Solution[] initialPopulation = new Solution[this.populationSize];
        int solutionCount = 0;

        // ArrayList<Solution> initialPopulation = new ArrayList<Solution>();

        double[] numberOfMembers = new double[4];
        if (feasibilityConstrained && !stabilityConstrained && !orientationConstrained) {
            numberOfMembers = new double[]{6, 9, 12, 15};
        }
        else if (stabilityConstrained && !feasibilityConstrained && !orientationConstrained) {
            numberOfMembers = new double[]{18, 21, 24, 27};
        }
        else if (feasibilityConstrained && stabilityConstrained && !orientationConstrained) {
            numberOfMembers = new double[]{6, 12, 18, 24};
        }

        Random rand = new Random();
        boolean val;

        if (feasibilityConstrained || stabilityConstrained && !orientationConstrained) {
            for (int i = 0;i < 4;i++) {
                double TrueProb = numberOfMembers[i]/32;
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

                    // initialPopulation.add(solution);
                }
            }
        }

        ImproveOrientation improveOrientation;
        double targetOrientation;
        if (!feasibilityConstrained && !stabilityConstrained) {
            improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, sideNodeNumber);
            targetOrientation = improveOrientation.findTargetOrientation(targetStiffnessRatio);
            int m = 0;
            while (m < populationSize) {
                Solution solution = this.problem.newSolution();
                for (int n = 0; n < solution.getNumberOfVariables(); ++n) {
                    val = rand.nextFloat() < 0.5;
                    BinaryVariable var = new BinaryVariable(1);
                    EncodingUtils.setBoolean(var, val);
                    solution.setVariable(n, var);
                }
                TrussRepeatableArchitecture trussArch = new TrussRepeatableArchitecture(solution);
                int[][] connectivityArray = trussArch.getConnectivityArrayFromSolution(solution);

                double[] architectureOrientations = improveOrientation.findMemberOrientations(connectivityArray);

                double totalOrientation = 0;
                for (int i = 0; i < architectureOrientations.length; i++) {
                    totalOrientation += architectureOrientations[i];
                }
                double meanOrientation = totalOrientation / architectureOrientations.length;
                if (Math.abs(meanOrientation - targetOrientation) < acceptanceMargin) {
                    initialPopulation[m] = solution;
                    m++;
                }
            }
        } else if (orientationConstrained) {
            improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, sideNodeNumber);
            targetOrientation = improveOrientation.findTargetOrientation(targetStiffnessRatio);
            for (int i = 0;i < 4;i++) {
                double TrueProb = numberOfMembers[i]/32;
                int j = 0;
                while (j < solutionsPerGroup) {
                    Solution solution = this.problem.newSolution();
                    for (int k = 0;k < solution.getNumberOfVariables();++k) {
                        val = rand.nextFloat() < TrueProb;
                        BinaryVariable var = new BinaryVariable(1);
                        EncodingUtils.setBoolean(var,val);
                        solution.setVariable(k,var);
                    }
                    TrussRepeatableArchitecture trussArch = new TrussRepeatableArchitecture(solution);
                    int[][] connectivityArray = trussArch.getConnectivityArrayFromSolution(solution);

                    double[] architectureOrientations = improveOrientation.findMemberOrientations(connectivityArray);

                    double totalOrientation = 0;
                    for (int s = 0; s < architectureOrientations.length; s++) {
                        totalOrientation += architectureOrientations[s];
                    }
                    double meanOrientation = totalOrientation / architectureOrientations.length;
                    if (Math.abs(meanOrientation - targetOrientation) < acceptanceMargin) {
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

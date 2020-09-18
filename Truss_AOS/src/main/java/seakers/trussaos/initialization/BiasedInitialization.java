package seakers.trussaos.initialization;

import org.moeaframework.core.Initialization;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
//import org.moeaframework.core.Variable;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.EncodingUtils;

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

    public BiasedInitialization(Problem problem, int populationSize) {
        this.problem = problem;
        this.populationSize = populationSize;
    }

    @Override
    public Solution[] initialize() {
        int solutionsPerGroup = (int) Math.round(populationSize/4);

        Solution[] initialPopulation = new Solution[this.populationSize];
        int solutionCount = 0;

        // ArrayList<Solution> initialPopulation = new ArrayList<Solution>();

        double[] numberOfMembers = {6,9,12,15};
        Random rand = new Random();
        boolean val;

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

        // return (Solution[]) initialPopulation.toArray();
        return initialPopulation;
    }
}

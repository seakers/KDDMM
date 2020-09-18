package seakers.trussaos.constraints;

import java.io.Serializable;
import java.util.*;

import org.moeaframework.core.EpsilonBoxDominanceArchive;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.comparator.DominanceComparator;

/**
 * This compound operator applies a series of operators with a specified
 * probability in a random order.
 *
 * The computeProbability method is taken from the SetConsistency class from the EOSS-dev repo
 *
 * @author nozomihitomi, modified by @roshan94
 */
public class KnowledgeStochasticRanking implements DominanceComparator,
        Serializable {

    private static final long serialVersionUID = 3653864833347649396L;

    /**
     * The number of constraints to apply
     */
    private final int numberConstraints;

    /**
     * the probabilities with which to apply the constraint (string property of
     * architecture)
     */
    private final HashMap<String, Double> probabilities;

    /**
     *
     * Current archive
     */
    private EpsilonBoxDominanceArchive archive;

    public KnowledgeStochasticRanking(int numberOperators, EpsilonBoxDominanceArchive currentArchive) {
        super();
        this.numberConstraints = numberOperators;
        this.probabilities = new HashMap<>(numberOperators);
        this.archive = currentArchive;
    }

    public KnowledgeStochasticRanking(int numberOperators, Collection<String> constraints, EpsilonBoxDominanceArchive currentArchive) {
        this.numberConstraints = numberOperators;
        this.archive = currentArchive;
        this.probabilities = new HashMap<>(numberOperators);

        for (String str : constraints) {
            this.probabilities.put(str, 1.0);
        }
    }

    /**
     * The constraint name that is a property of the architecture
     *
     * @param constraint name that is a property of the architecture
     */
    public void appendConstraint(String constraint) {
        probabilities.put(constraint, 1.0);
    }

    /**
     * Updates the probability of applying the constraint
     *
     * @param constraint name that is a property of the architecture
     * @param probability probability of applying the constraint
     */
    public void updateProbability(String constraint, double probability) {
        probabilities.replace(constraint, probability);
    }

    public Collection<String> getConstraints() {
        return probabilities.keySet();
    }

    @Override
    public int compare(Solution solution1, Solution solution2) {
        double constraint1 = 0;
        double constraint2 = 0;
        int numApplied = 0;
        ArrayList<String> constraints = new ArrayList<>(probabilities.keySet());
        Collections.shuffle(constraints);
        for (String str : constraints) {
            if (numApplied >= numberConstraints) {
                break;
            }
            double strProbability = computeProbability(str);
            if (PRNG.nextDouble() < strProbability) {
                constraint1 += (double) solution1.getAttribute(str);
                constraint2 += (double) solution2.getAttribute(str);
            }
            updateProbability(str, strProbability);
        }

        if (constraint1 < constraint2) {
            return -1;
        } else if (constraint1 > constraint2) {
            return 1;
        } else {
            return 0;
        }
    }

    public double computeProbability(String constraintString) {
        double constraintProbability = 0;
        int consistentCount = 0;
        for (Solution s : archive) {
            if((double)s.getAttribute(constraintString) == 0){
                consistentCount++;
            }
        }
        constraintProbability = Math.max((double)consistentCount / (double)archive.getNumberOfDominatingImprovements(), 0.03);

        return constraintProbability;
    }
}

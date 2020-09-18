package seakers.trussaos.constraints;

import java.io.Serializable;
import java.util.Collection;
import org.moeaframework.core.Solution;
import org.moeaframework.core.comparator.DominanceComparator;

/**
 * This checks to see that at least one constraint is satisfied.
 *
 * @author nozomihitomi
 */
public class DisjunctiveNormalForm implements DominanceComparator,
        Serializable {

    private static final long serialVersionUID = 7884934515386899671L;

    private final Collection<String> constraints;

    public DisjunctiveNormalForm(Collection<String> constraints) {
        this.constraints = constraints;
    }

    @Override
    public int compare(Solution solution1, Solution solution2) {
        double constraint1 = 0;
        double constraint2 = 0;
        for (String str : constraints) {
            double violation = (double) solution1.getAttribute(str);
            if (violation == 0) {
                constraint1 = 0;
                break;
            } else {
                constraint1 += violation;
            }
        }
        for (String str : constraints) {
            double violation = (double) solution2.getAttribute(str);
            if (violation == 0) {
                constraint2 = 0;
                break;
            } else {
                constraint2 += violation;
            }
        }

        if (constraint1 < constraint2) {
            return -1;
        } else if (constraint1 > constraint2) {
            return 1;
        } else {
            return 0;
        }
    }
}
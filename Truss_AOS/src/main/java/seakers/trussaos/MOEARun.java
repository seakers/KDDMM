package seakers.trussaos;

import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.CompoundVariation;
import org.moeaframework.core.operator.OnePointCrossover;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.util.TypedProperties;
import seakers.aos.aos.AOSMOEA;
import seakers.aos.creditassignment.offspringparent.OffspringParentDomination;
import seakers.aos.operator.AOSVariation;
import seakers.aos.operator.AOSVariationOP;
import com.mathworks.engine.*;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.core.operator.RandomInitialization;
import org.moeaframework.core.operator.TournamentSelection;
//import org.moeaframework.core.operator.real.PM;
//import org.moeaframework.core.operator.real.SBX;
import seakers.aos.operatorselectors.OperatorSelector;
import seakers.aos.operatorselectors.ProbabilityMatching;
import seakers.trussaos.initialization.BiasedInitialization;
import seakers.trussaos.operators.AddMember;
import seakers.trussaos.operators.RemoveIntersection;
import seakers.trussaos.constraints.DisjunctiveNormalForm;
import seakers.trussaos.constraints.KnowledgeStochasticRanking;

/**
 * Executable class for the eMOEA run with/without AOS for the Truss Optimization problem.
 *
 * @author roshan94
 */

public class MOEARun {

    /**
     * pool of resources
     */
    private static ExecutorService pool;

    /**
     * Executor completion services helps remove completed tasks
     */
    private static CompletionService<Algorithm> ecs;

    /**
     * Matlab Engine for function evaluation
     */
    private static MatlabEngine engine;

    public static void main(String[] args) throws InterruptedException, ExecutionException, EngineException {

        // Define problem parameters
        //String csvPath = "C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\src\\main\\java\\seakers\\trussaos";
        String csvPath = System.getProperty("user.dir");
        double targetStiffnessRatio = 1;
        boolean useFibreStiffness = false;

        boolean biasedInitialization = true;

        /**
         * Mode = 1 for conventional e-MOEA run
         *      = 2 for AOS MOEA run
         *      = 3 for soft constraints
         */
        int mode = 1;

        /**
         * soft_con = 1 for DNF
         *          = 2 for ACH
         */
        int soft_con = 2;

        int numCPU = 1;
        int numRuns = 1;
        pool = Executors.newFixedThreadPool(numCPU);
        ecs = new ExecutorCompletionService<>(pool);
        engine = MatlabEngine.startMatlab();

        //String userPathOutput = "";
        //userPathOutput = engine.feval("userpath",csvPath);

        //create the desired algorithm
        //parameters and operators for search
        double[] epsilonDouble = new double[]{0.01, 0.01};
        TypedProperties properties = new TypedProperties();
        //search paramaters set here
        int popSize = 100;
        int maxEvals = 10000;
        properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);
        double mutationProbability = 1. / 32.;
        properties.setDouble("mutationProbability", mutationProbability);
        Variation singlecross;
        Variation bitFlip;
        Initialization initialization;

        double crossoverProbability = 1.0;

        // For AOS MOEA Run
        boolean useStability = false;
        boolean useFeasibility = false;

        for (int i = 0; i < numRuns; i++) {

            // Create a new problem class
            TrussAOSProblem trussProblem = new TrussAOSProblem(csvPath,useFibreStiffness,targetStiffnessRatio,engine);

            if (biasedInitialization){
                initialization = new BiasedInitialization(trussProblem, popSize);
            }
            else {
                initialization = new RandomInitialization(trussProblem, popSize);
            }

            // Initialize population structure for algorithm
            Population population = new Population();

            EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
            ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
            TournamentSelection selection = new TournamentSelection(2, comp);

            switch(mode) {
                case 1: // Conventional Epsilon MOEA Run
                    properties.setDouble("crossoverProbability", crossoverProbability);
                    singlecross = new OnePointCrossover(crossoverProbability);
                    bitFlip = new BitFlip(mutationProbability);
                    CompoundVariation var = new CompoundVariation(singlecross, bitFlip);

                    Algorithm eMOEA = new EpsilonMOEA(trussProblem, population, archive, selection, var, initialization);
                    ecs.submit(new EvolutionarySearch(eMOEA, properties, csvPath + File.separator + "result", "emoea" + String.valueOf(i) + "_biasedinit", engine, useFibreStiffness, targetStiffnessRatio));
                    break;

                case 2: // AOS MOEA Run

                    //setup for saving results
                    properties.setBoolean("saveQuality", true);
                    properties.setBoolean("saveCredits", true);
                    properties.setBoolean("saveSelection", true);

                    // TESTING
                    //double onePointCrossoverProbability = 1.0;
                    //double twoPointCrossoverProbability = 0.8;
                    //double halfUniformCrossoverProbability = 0.75;
                    //double uniformCrossoverProbability = 0.6;
                    //ArrayList<Variation> operators = new ArrayList<>();
                    //add domain-independent heuristics
                    //operators.add(new CompoundVariation(new OnePointCrossover(onePointCrossoverProbability), new BitFlip(mutationProbability)));
                    //operators.add(new CompoundVariation(new TwoPointCrossover(twoPointCrossoverProbability), new BitFlip(mutationProbability)));
                    //operators.add(new CompoundVariation(new HUX(halfUniformCrossoverProbability), new BitFlip(mutationProbability)));
                    //operators.add(new CompoundVariation(new UniformCrossover(uniformCrossoverProbability), new BitFlip(mutationProbability)));

                    // IMPLEMENTATION WITH ACTUAL REPAIR OPERATORS
                    double[][] globalNodePositions = trussProblem.getNodalConnectivityArray();
                    ArrayList<Variation> operators = new ArrayList<>();
                    Variation addMember = new AddMember(useFeasibility, engine, globalNodePositions);
                    Variation removeIntersection = new RemoveIntersection(useStability, engine, globalNodePositions);
                    operators.add(new CompoundVariation(new OnePointCrossover(crossoverProbability), new BitFlip(mutationProbability)));
                    operators.add(addMember);
                    operators.add(removeIntersection);

                    // DisjunctiveNormalForm dnf = new DisjunctiveNormalForm(constraints);
                    // EpsilonKnoweldgeConstraintComparator epskcc = new EpsilonKnoweldgeConstraintComparator(epsilonDouble, dnf)

                    properties.setDouble("pmin", 0.1);

                    //create operator selector
                    //OperatorSelector operatorSelector = new AdaptivePursuit(operators, 0.8, 0.8, 0.1);
                    OperatorSelector operatorSelector = new ProbabilityMatching(operators, 0.6, 0.1);

                    //create credit assignment
                    //SetImprovementDominance creditAssignment = new SetImprovementDominance(archive, 1, 0);
                    OffspringParentDomination creditAssignment = new OffspringParentDomination(1.0,0.5,0.0);

                    //create AOS
                    //AOSVariation aosStrategy = new AOSVariationSI(operatorSelector, creditAssignment, popSize);
                    AOSVariation aosStrategy = new AOSVariationOP(operatorSelector,creditAssignment,popSize);

                    EpsilonMOEA emoea = new EpsilonMOEA(trussProblem, population, archive,
                            selection, aosStrategy, initialization, comp);
                    AOSMOEA aos = new AOSMOEA(emoea, aosStrategy, true);

                    aos.setName("constraint_adaptive");
                    ecs.submit(new EvolutionarySearch(aos, properties, csvPath + File.separator + "result", aos.getName() + String.valueOf(i), engine, useFibreStiffness, targetStiffnessRatio));
                    break;

                case 3: // Soft Constraints
                    if (soft_con == 1){
                        globalNodePositions = trussProblem.getNodalConnectivityArray();
                        //operators = new ArrayList<>();
                        addMember = new AddMember(useFeasibility, engine, globalNodePositions);
                        removeIntersection = new RemoveIntersection(useStability, engine, globalNodePositions);
                        //operators.add(new CompoundVariation(new OnePointCrossover(crossoverProbability), new BitFlip(mutationProbability)));
                        //operators.add(addMember);
                        //operators.add(removeIntersection);

                        properties.setDouble("crossoverProbability", crossoverProbability);
                        singlecross = new OnePointCrossover(crossoverProbability);
                        bitFlip = new BitFlip(mutationProbability);
                        var = new CompoundVariation(singlecross, bitFlip);

                        HashMap<Variation, String> constraintOperatorMap = new HashMap<>();
                        constraintOperatorMap.put(addMember, "StabilityViolation");
                        constraintOperatorMap.put(removeIntersection, "FeasibilityViolation");

                        HashSet<String> constraints = new HashSet<>(constraintOperatorMap.values());

                        DisjunctiveNormalForm dnf = new DisjunctiveNormalForm(constraints);
                        // EpsilonKnoweldgeConstraintComparator epskcc = new EpsilonKnoweldgeConstraintComparator(epsilonDouble, dnf);

                        selection = new TournamentSelection(2, dnf);

                        emoea = new EpsilonMOEA(trussProblem, population, archive, selection, var, initialization, dnf);
                        ecs.submit(new EvolutionarySearch(emoea, properties, csvPath + File.separator + "result", "emoea_dnf" + String.valueOf(i), engine, useFibreStiffness, targetStiffnessRatio));
                        break;
                    }
                    else if (soft_con == 2) {
                        properties.setDouble("crossoverProbability", crossoverProbability);
                        singlecross = new OnePointCrossover(crossoverProbability);
                        bitFlip = new BitFlip(mutationProbability);
                        var = new CompoundVariation(singlecross, bitFlip);

                        globalNodePositions = trussProblem.getNodalConnectivityArray();
                        addMember = new AddMember(useFeasibility, engine, globalNodePositions);
                        removeIntersection = new RemoveIntersection(useStability, engine, globalNodePositions);

                        HashMap<Variation, String> constraintOperatorMap = new HashMap<>();
                        constraintOperatorMap.put(addMember, "StabilityViolation");
                        constraintOperatorMap.put(removeIntersection, "FeasibilityViolation");

                        HashSet<String> constraints = new HashSet<>(constraintOperatorMap.values());

                        KnowledgeStochasticRanking ksr = new KnowledgeStochasticRanking(constraintOperatorMap.size(), constraintOperatorMap.values(), archive);

                        selection = new TournamentSelection(2, ksr);

                        emoea = new EpsilonMOEA(trussProblem, population, archive, selection, var, initialization, ksr);

                        ecs.submit(new EvolutionarySearch(emoea, properties, csvPath + File.separator + "result", "emoea_ach" + String.valueOf(i), engine, useFibreStiffness, targetStiffnessRatio));
                        break;
                    }

                default:
                    throw new IllegalArgumentException(String.format("%d is an invalid option", mode));

            }
        }

        for (int i = 0; i < numRuns; i++) {
            try {
                Algorithm alg = ecs.take().get();
            } catch (InterruptedException | ExecutionException ex) {
                Logger.getLogger(MOEARun.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        pool.shutdown();
        engine.close();

    }
}

package seakers.trussaos;

import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.CompoundVariation;
import org.moeaframework.core.operator.OnePointCrossover;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.core.operator.TwoPointCrossover;
import org.moeaframework.core.operator.binary.HUX;
import org.moeaframework.core.operator.UniformCrossover;
import org.moeaframework.util.TypedProperties;
import seakers.aos.creditassignment.setimprovement.SetImprovementDominance;
import seakers.aos.history.AOSHistoryIO;
import seakers.aos.aos.AOSMOEA;
import seakers.aos.creditassignment.offspringparent.OffspringParentDomination;
import seakers.aos.operator.AOSVariation;
import seakers.aos.operator.AOSVariationOP;
import seakers.aos.operator.AOSVariationSI;
import seakers.aos.operatorselectors.AdaptivePursuit;
import seakers.aos.operatorselectors.ProbabilityMatching;
import com.mathworks.engine.*;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.moeaframework.Instrumenter;
import org.moeaframework.algorithm.EpsilonMOEA;
import org.moeaframework.analysis.collector.InstrumentedAlgorithm;
import org.moeaframework.core.operator.RandomInitialization;
import org.moeaframework.core.operator.TournamentSelection;
//import org.moeaframework.core.operator.real.PM;
//import org.moeaframework.core.operator.real.SBX;
import seakers.aos.operatorselectors.OperatorSelector;
import seakers.trussaos.operators.AddTruss;
import seakers.trussaos.operators.RemoveIntersection;

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
        boolean useFibreStiffness = true;

        /**
         * Mode = 1 for conventional e-MOEA run
         *      = 2 for AOS MOEA run
         */
        int mode = 2;

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

        for (int i = 0; i < numRuns; i++) {

            // Create a new problem class
            TrussAOSProblem trussProblem = new TrussAOSProblem(csvPath,useFibreStiffness,targetStiffnessRatio,engine);

            initialization = new RandomInitialization(trussProblem, popSize);

            // Initialize population structure for algorithm
            Population population = new Population();
            //KnowledgeStochasticRanking ksr = new KnowledgeStochasticRanking(constraintOperatorMap.size(), constraintOperatorMap.values());
            //DisjunctiveNormalForm dnf = new DisjunctiveNormalForm(constraints);
            //EpsilonKnoweldgeConstraintComparator epskcc = new EpsilonKnoweldgeConstraintComparator(epsilonDouble, dnf);
            EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);
            ChainedComparator comp = new ChainedComparator(new ParetoObjectiveComparator());
            TournamentSelection selection = new TournamentSelection(2, comp);

            // For AOS MOEA Run
            boolean useStability = false;
            boolean useFeasibility = false;

            switch(mode) {
                case 1: // Conventional Epsilon MOEA Run
                    double crossoverProbability = 1.0;
                    properties.setDouble("crossoverProbability", crossoverProbability);
                    singlecross = new OnePointCrossover(crossoverProbability);
                    bitFlip = new BitFlip(mutationProbability);
                    CompoundVariation var = new CompoundVariation(singlecross, bitFlip);

                    Algorithm eMOEA = new EpsilonMOEA(trussProblem, population, archive, selection, var, initialization);
                    ecs.submit(new EvolutionarySearch(eMOEA, properties, csvPath + File.separator + "result", "emoea" + String.valueOf(i), engine, useFibreStiffness, targetStiffnessRatio));
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
                    double onePointCrossoverProbability = 1.0;
                    double[][] globalNodePositions = trussProblem.getNodalConnectivityArray();
                    ArrayList<Variation> operators = new ArrayList<>();
                    Variation addTruss = new AddTruss(useFeasibility, engine, globalNodePositions);
                    Variation removeIntersection = new RemoveIntersection(useStability, engine, globalNodePositions);
                    operators.add(new CompoundVariation(new OnePointCrossover(onePointCrossoverProbability), new BitFlip(mutationProbability)));
                    operators.add(addTruss);
                    operators.add(removeIntersection);

                    //HashMap<Variation, String> constraintOperatorMap = new HashMap<>();
                    //constraintOperatorMap.put(addTruss, "stabilityViolationSum");
                    //constraintOperatorMap.put(removeIntersectipn, "FeasibilityViolationSum");

                    //HashSet<String> constraints = new HashSet<>(constraintOperatorMap.values());

                    //DisjunctiveNormalForm dnf = new DisjunctiveNormalForm(constraints);
                    //EpsilonKnoweldgeConstraintComparator epskcc = new EpsilonKnoweldgeConstraintComparator(epsilonDouble, dnf)

                    properties.setDouble("pmin", 0.03);
                    //create operator selector
                    OperatorSelector operatorSelector = new AdaptivePursuit(operators, 0.8, 0.8, 0.03);
                    //create credit assignment
                    SetImprovementDominance creditAssignment = new SetImprovementDominance(archive, 1, 0);
                    //create AOS
                    AOSVariation aosStrategy = new AOSVariationSI(operatorSelector, creditAssignment, popSize);
                    EpsilonMOEA emoea = new EpsilonMOEA(trussProblem, population, archive,
                            selection, aosStrategy, initialization, comp);
                    AOSMOEA aos = new AOSMOEA(emoea, aosStrategy, true);

                    aos.setName("constraint_adaptive");
                    ecs.submit(new EvolutionarySearch(aos, properties, csvPath + File.separator + "result", aos.getName() + String.valueOf(i), engine, useFibreStiffness, targetStiffnessRatio));
                    break;

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

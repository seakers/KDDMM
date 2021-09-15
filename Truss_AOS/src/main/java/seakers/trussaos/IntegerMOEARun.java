package seakers.trussaos;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.AggregateConstraintComparator;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.real.UM;
import org.moeaframework.problem.AbstractProblem;
import org.moeaframework.util.TypedProperties;
import seakers.aos.aos.AOSMOEA;
import seakers.aos.creditassignment.offspringparent.OffspringParentDomination;
import seakers.aos.operator.AOSVariation;
import seakers.aos.operator.AOSVariationOP;
import com.mathworks.engine.*;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.*;

import org.moeaframework.algorithm.EpsilonMOEA;
import seakers.aos.operatorselectors.OperatorSelector;
import seakers.aos.operatorselectors.ProbabilityMatching;
import seakers.trussaos.initialization.BiasedInitializationIntegerRadii;
import seakers.trussaos.constrainthandling.KnowledgeStochasticRanking;
import seakers.trussaos.operators.integerradii.AddDiagonalMemberIntegerRadii;
import seakers.trussaos.operators.integerradii.AddMemberIntegerRadii;
import seakers.trussaos.operators.integerradii.ImproveOrientationIntegerRadii;
import seakers.trussaos.operators.integerradii.RemoveIntersectionIntegerRadii;
import seakers.trussaos.problems.IntegerTrussProblem;
import seakers.trussaos.initialization.SynchronizedMersenneTwister;

/**
 * Executable class for different eMOEA run experiments for the Truss Optimization problem (Integer encoded radii).
 *
 * @author roshan94
 */

public class IntegerMOEARun {

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

    public static void main(String[] args) throws InterruptedException, EngineException {
        String csvPath = System.getProperty("user.dir");
        double targetStiffnessRatio = 1;

        /**
         * modelChoice = 0 --> Fibre Stiffness Model
         *             = 1 --> Truss Stiffness Model
         *             = 2 --> Beam Model
         */
        int modelChoice = 1;

        // Heuristic Enforcement Methods
        /**
         * partialCollapsibilityConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * nodalPropertiesConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * orientationConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         * intersectionConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
         *
         * heuristicsConstrained = [partialCollapsibilityConstrained, nodalPropertiesConstrained, orientationConstrained, intersectionConstrained]
         */
        boolean[] partialCollapsibilityConstrained = {false, false, false, false, false, false};
        boolean[] nodalPropertiesConstrained = {false, false, false, false, false, false};
        boolean[] orientationConstrained = {false, false, false, false, false, false};
        boolean[] intersectionConstrained = {false, false, false, false, false, false};

        boolean[][] heuristicsConstrained = new boolean[4][6];
        for (int i = 0; i < 6; i++) {
            heuristicsConstrained[0][i] = partialCollapsibilityConstrained[i];
            heuristicsConstrained[1][i] = nodalPropertiesConstrained[i];
            heuristicsConstrained[2][i] = orientationConstrained[i];
            heuristicsConstrained[3][i] = intersectionConstrained[i];
        }

        int numberOfHeuristicConstraints = 0;
        int numberOfHeuristicObjectives = 0;
        for (int i = 0; i < 4; i++) {
            if (heuristicsConstrained[i][5]) {
                numberOfHeuristicConstraints++;
            }
            if (heuristicsConstrained[i][4]) {
                numberOfHeuristicObjectives++;
            }
        }

        int numCPU = 1;
        int numRuns = 1;
        pool = Executors.newFixedThreadPool(numCPU);
        ecs = new ExecutorCompletionService<>(pool);
        engine = MatlabEngine.startMatlab();

        //String userPathOutput = "";
        //userPathOutput = engine.feval("userpath",csvPath);

        // Create the desired algorithm parameters and operators for search
        double[] epsilonDouble = new double[]{0.01, 0.01};
        TypedProperties properties = new TypedProperties();
        // Search paramaters set here
        int popSize = 100;
        int maxEvals = 5000;
        properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);

        Variation crossOver;
        Variation mutation;
        Initialization initialization;

        double crossoverProbability = 1.0;

        // For AOS MOEA Run
        boolean maintainStability = false;
        boolean maintainFeasibility = false;

        String fileSaveNameModel = "_beam";
        if (modelChoice == 0) {
            fileSaveNameModel = "_fibre";
        } else if (modelChoice == 1){
            fileSaveNameModel = "_truss";
        }

        String[] heuristicAbbreviations = {"p","n","o","i"};
        StringBuilder fileSaveNameConstraint = new StringBuilder();
        for (int i = 0; i < heuristicsConstrained[0].length; i++) {
            StringBuilder enforcedHeuristics = new StringBuilder();
            int heuristicCount = 0;
            for (int j = 0; j < heuristicsConstrained.length; j++) {
                if (heuristicsConstrained[j][i]) {
                    enforcedHeuristics.append(heuristicAbbreviations[j]);
                    heuristicCount++;
                }
            }
            if (heuristicCount > 0) {
                fileSaveNameConstraint.append(enforcedHeuristics.toString()).append("con").append(Integer.toString(i)).append("_");
            }
        }

        // New dimensions for printable solutions
        double printableSideLength = 10e-3; // in m
        double sideNodeNumber = 3.0D;
        int nucFactor = 3; // Not used if PBC model is used
        double[] radiusChoices = new double[]{0, 250e-6, 500e-6, 750e-6}; // in m
        double[] YoungsModulusChoices = new double[]{1.8162e6, 3.5e6, 4.5e6, 7.5e6}; // in Pa

        int totalNumberOfMembers;
        if (sideNodeNumber >= 5) {
            int sidenumSquared = (int) (sideNodeNumber*sideNodeNumber);
            totalNumberOfMembers =  sidenumSquared * (sidenumSquared - 1)/2;
        }
        else {
            totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sideNodeNumber*sideNodeNumber))/(CombinatoricsUtils.factorial((int) ((sideNodeNumber*sideNodeNumber) - 2)) * CombinatoricsUtils.factorial(2)));
        }
        int numberOfRepeatableMembers = (int) (2 * (CombinatoricsUtils.factorial((int) sideNodeNumber)/(CombinatoricsUtils.factorial((int) (sideNodeNumber - 2)) * CombinatoricsUtils.factorial(2))));
        int numTrussVariables = totalNumberOfMembers - numberOfRepeatableMembers;

        double mutationProbability = 1. / (numTrussVariables + 1);
        properties.setDouble("mutationProbability", mutationProbability);

        PRNG.setRandom(new SynchronizedMersenneTwister());

        for (int i = 0; i < numRuns; i++) {

            // Problem class for printable solutions
            AbstractProblem trussProblem;
            double[][] globalNodePositions;

            trussProblem = new IntegerTrussProblem(csvPath, modelChoice, numTrussVariables, numberOfHeuristicObjectives, numberOfHeuristicConstraints, targetStiffnessRatio, engine, radiusChoices, printableSideLength, YoungsModulusChoices, sideNodeNumber, nucFactor, heuristicsConstrained);
            globalNodePositions = ((IntegerTrussProblem) trussProblem).getNodalConnectivityArray();

            //String fileSaveNameInit;
            if (partialCollapsibilityConstrained[2] || nodalPropertiesConstrained[2] || orientationConstrained[2] || intersectionConstrained[2]) { // Initialization object
                initialization = new BiasedInitializationIntegerRadii(trussProblem, popSize, radiusChoices, YoungsModulusChoices, partialCollapsibilityConstrained[2], nodalPropertiesConstrained[2], orientationConstrained[2], intersectionConstrained[2], globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
            } else {
                initialization = new RandomInitialization(trussProblem, popSize);
            }

            // Initialize population structure for algorithm
            Population population = new Population();

            EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);

            Algorithm moeaObj;
            ChainedComparator comp = null;
            TournamentSelection selection;
            CompoundVariation var = null;
            AOSVariation aosStrategy = null;

            // Initialize heuristic operators
            AddMemberIntegerRadii addMember = new AddMemberIntegerRadii(maintainFeasibility, radiusChoices, YoungsModulusChoices, engine, globalNodePositions, sideNodeNumber, printableSideLength, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
            RemoveIntersectionIntegerRadii removeIntersection = new RemoveIntersectionIntegerRadii(maintainStability, engine, globalNodePositions, sideNodeNumber, printableSideLength, radiusChoices, YoungsModulusChoices, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
            AddDiagonalMemberIntegerRadii addDiagonalMember = new AddDiagonalMemberIntegerRadii(maintainFeasibility, radiusChoices, YoungsModulusChoices, engine, globalNodePositions, sideNodeNumber, printableSideLength, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
            ImproveOrientationIntegerRadii improveOrientation = new ImproveOrientationIntegerRadii(globalNodePositions, radiusChoices, YoungsModulusChoices, targetStiffnessRatio, (int) sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);

            Variation[] heuristicOperators = {addDiagonalMember, addMember, improveOrientation, removeIntersection};
            String[] heuristicAttributes = {"PartialCollapsibilityViolation","NodalPropertiesViolation","OrientationViolation","IntersectionViolation"};

            if (partialCollapsibilityConstrained[3] || nodalPropertiesConstrained[3] || orientationConstrained[3] || intersectionConstrained[3]) { // Adaptive Constraint Handling objects
                KnowledgeStochasticRanking ksr;

                HashMap<Variation, String> constraintOperatorMap = new HashMap<>();

                for (int j = 0; j < heuristicsConstrained.length; j++) {
                    if (heuristicsConstrained[j][3]) {
                        constraintOperatorMap.put(heuristicOperators[j], heuristicAttributes[j]);
                    }
                }

                ksr = new KnowledgeStochasticRanking(constraintOperatorMap.size(), constraintOperatorMap.values(), archive);
                comp = new ChainedComparator(new AggregateConstraintComparator(), ksr, new ParetoObjectiveComparator());

                selection = new TournamentSelection(2, comp);

                //moeaObj = new EpsilonMOEA(trussProblem, population, archive, selection, var, initialization, comp);
            }
            else { // Epsilon MOEA objects
                comp = new ChainedComparator(new AggregateConstraintComparator(), new ParetoObjectiveComparator());
                selection = new TournamentSelection(2, comp);
            }

            if (partialCollapsibilityConstrained[1] || nodalPropertiesConstrained[1] || orientationConstrained[1] || intersectionConstrained[1]) { // AOS objects
                //comp = new ChainedComparator(new ParetoObjectiveComparator());
                //selection = new TournamentSelection(2, comp);

                // Setup for saving results
                properties.setBoolean("saveQuality", true);
                properties.setBoolean("saveCredits", true);
                properties.setBoolean("saveSelection", true);

                // IMPLEMENTATION WITH ACTUAL REPAIR OPERATORS

                ArrayList<Variation> operators = new ArrayList<>();
                for (int k = 0; k < heuristicsConstrained.length; k++) {
                    if (heuristicsConstrained[k][1]) {
                        operators.add(heuristicOperators[k]);
                    }
                }

                operators.add(new CompoundVariation(new UniformCrossover(crossoverProbability), new UM(mutationProbability)));

                properties.setDouble("pmin", 0.03);

                // Create operator selector
                //OperatorSelector operatorSelector = new AdaptivePursuit(operators, 0.8, 0.8, 0.1);
                OperatorSelector operatorSelector = new ProbabilityMatching(operators, 0.6, 0.03);

                // Create credit assignment
                //SetImprovementDominance creditAssignment = new SetImprovementDominance(archive, 1, 0);
                OffspringParentDomination creditAssignment = new OffspringParentDomination(1.0, 0.5, 0.0);

                // Create AOS
                //aosStrategy = new AOSVariationSI(operatorSelector, creditAssignment, popSize);
                aosStrategy = new AOSVariationOP(operatorSelector, creditAssignment, popSize);

            }
            else { // Epsilon MOEA objects
                properties.setDouble("crossoverProbability", crossoverProbability);
                //crossOver = new OnePointCrossover(crossoverProbability);
                crossOver = new UniformCrossover(crossoverProbability);
                mutation = new UM(mutationProbability);
                var = new CompoundVariation(crossOver, mutation);
            }

            // Creating AOS MOEA object if needed
            if (partialCollapsibilityConstrained[1] || nodalPropertiesConstrained[1] || orientationConstrained[1] || intersectionConstrained[1]) {
                EpsilonMOEA emoea = new EpsilonMOEA(trussProblem, population, archive, selection, aosStrategy, initialization, comp);
                moeaObj = new AOSMOEA(emoea, aosStrategy, true);
            }
            else {
                moeaObj = new EpsilonMOEA(trussProblem, population, archive, selection, var, initialization, comp);
            }

            ecs.submit(new IntegerEvolutionarySearch(moeaObj, properties, csvPath + File.separator + "result", "emoea_" + String.valueOf(i) + fileSaveNameConstraint.toString() + fileSaveNameModel, sideNodeNumber, radiusChoices, YoungsModulusChoices, numberOfHeuristicObjectives, numberOfHeuristicConstraints));
        }

    }
}

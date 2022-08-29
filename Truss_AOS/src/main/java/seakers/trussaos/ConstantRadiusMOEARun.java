package seakers.trussaos;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.AggregateConstraintComparator;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.problem.AbstractProblem;
import org.moeaframework.util.TypedProperties;
import seakers.ahs.AHSMOEA;
import seakers.ahs.AdaptiveHeuristicSelection;
import seakers.aos.aos.AOSMOEA;
import seakers.aos.creditassignment.offspringparent.OffspringParentDomination;
import seakers.aos.creditassignment.setcontribution.SetContributionDominance;
import seakers.aos.creditassignment.setimprovement.SetImprovementDominance;
import seakers.aos.operator.AOSVariation;
import seakers.aos.operator.AOSVariationOP;
import com.mathworks.engine.*;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.moeaframework.algorithm.EpsilonMOEA;
import seakers.aos.operator.AOSVariationSC;
import seakers.aos.operator.AOSVariationSI;
import seakers.aos.operatorselectors.OperatorSelector;
import seakers.aos.operatorselectors.ProbabilityMatching;
import seakers.trussaos.initialization.BiasedInitialization;
import seakers.trussaos.initialization.SynchronizedMersenneTwister;
//import seakers.trussaos.operators.constantradii.RemoveIntersection;
import seakers.trussaos.operators.constantradii.RemoveIntersection2;
import seakers.trussaos.operators.constantradii.AddDiagonalMember;
import seakers.trussaos.operators.constantradii.AddMember;
//import seakers.trussaos.operators.constantradii.ImproveOrientation;
import seakers.trussaos.operators.constantradii.ImproveOrientation2;
import seakers.trussaos.constrainthandling.KnowledgeStochasticRanking;
import seakers.trussaos.problems.ConstantRadiusArteryProblem;
import seakers.trussaos.problems.ConstantRadiusTrussProblem;
import seakers.trussaos.problems.ConstantRadiusTrussProblem2;

/**
 * Executable class for different eMOEA run experiments for the Truss and Artery Optimization problems.
 *
 * @author roshan94
 */

public class ConstantRadiusMOEARun {

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

        /**
         * modelChoice = 0 --> Fibre Stiffness Model
         *             = 1 --> Truss Stiffness Model
         *             = 2 --> Beam Model
         */
        int modelChoice = 1; // Fibre stiffness model cannot be used for the artery problem

        boolean arteryProblem = true; // Solve the artery optimization (otherwise the original truss problem is solved)
        boolean useOptimizationProblem2 = true; // Use ConstantRadiusTrussProblem2 as problem class (instead of ConstantRadiusTrussProblem)

        double targetStiffnessRatio = 1;
        if (arteryProblem) {
            targetStiffnessRatio = 0.421;
        }

        // Heuristic Enforcement Methods
        /**
         * partialCollapsibilityConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint, AHS]
         * nodalPropertiesConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint, AHS]
         * orientationConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint, AHS]
         * intersectionConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint, AHS]
         *
         * heuristicsConstrained = [partialCollapsibilityConstrained, nodalPropertiesConstrained, orientationConstrained, intersectionConstrained]
         */
        boolean[] partialCollapsibilityConstrained = {true, false, false, false, false, false, false};
        boolean[] nodalPropertiesConstrained = {true, false, false, false, false, false, false};
        boolean[] orientationConstrained = {true, false, false, false, false, false, false};
        boolean[] intersectionConstrained = {true, false, false, false, false, false, false};

        // Bias initial population with low number of members
        boolean useLowMemberBiasing = false;

        boolean[][] heuristicsConstrained = new boolean[4][7];
        for (int i = 0; i < 7; i++) {
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
        int numRuns = 30;
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
        int maxEvals = 6000;
        properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);

        Variation crossOver;
        Variation bitMutation;
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

        String fileSaveNameProblem  = "";
        if (arteryProblem) {
            fileSaveNameProblem = "_artery";
        } else {
            if (useOptimizationProblem2) {
                fileSaveNameProblem = "_prob2";
            }
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
        double printableRadius = 250e-6; // in m
        double printableSideLength = 10e-3; // in m
        double printableModulus = 1.8162e6; // in Pa
        double sideNodeNumber = 3.0D;
        int nucFactor = 3; // Not used if PBC model is used

        int totalNumberOfMembers;
        if (sideNodeNumber >= 5) {
            int sidenumSquared = (int) (sideNodeNumber*sideNodeNumber);
            totalNumberOfMembers =  sidenumSquared * (sidenumSquared - 1)/2;
        }
        else {
            totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sideNodeNumber*sideNodeNumber))/(CombinatoricsUtils.factorial((int) ((sideNodeNumber*sideNodeNumber) - 2)) * CombinatoricsUtils.factorial(2)));
        }
        int numberOfRepeatableMembers = (int) (2 * (CombinatoricsUtils.factorial((int) sideNodeNumber)/(CombinatoricsUtils.factorial((int) (sideNodeNumber - 2)) * CombinatoricsUtils.factorial(2))));
        int numVariables = totalNumberOfMembers - numberOfRepeatableMembers;

        double mutationProbability = 1. / numVariables;
        properties.setDouble("mutationProbability", mutationProbability);

        PRNG.setRandom(new SynchronizedMersenneTwister());

        AdaptiveHeuristicSelection ahs = null;

        for (int i = 0; i < numRuns; i++) {
            AtomicInteger nfeval = new AtomicInteger(); // Used in the Adaptive Heuristic Handling method, updated in the Evolutionary Search algorithm
            HashSet<Solution> solutions = new HashSet<>(); // Used in the Adaptive Heuristic Handling method, updated in the Evolutionary Search algorithm

            // Problem class for printable solutions
            AbstractProblem trussProblem;
            double[][] globalNodePositions;
            if (arteryProblem) {
                trussProblem = new ConstantRadiusArteryProblem(csvPath, modelChoice, numVariables, numberOfHeuristicObjectives, numberOfHeuristicConstraints, printableRadius, printableSideLength, printableModulus, sideNodeNumber, nucFactor, targetStiffnessRatio, engine, heuristicsConstrained);
                globalNodePositions = ((ConstantRadiusArteryProblem) trussProblem).getNodalConnectivityArray();
            } else {
                if (useOptimizationProblem2) {
                    trussProblem = new ConstantRadiusTrussProblem2(csvPath, modelChoice, numVariables, numberOfHeuristicObjectives, numberOfHeuristicConstraints, printableRadius, printableSideLength, printableModulus, sideNodeNumber, nucFactor, targetStiffnessRatio, engine, heuristicsConstrained);
                    globalNodePositions = ((ConstantRadiusTrussProblem2) trussProblem).getNodalConnectivityArray();
                } else {
                    trussProblem = new ConstantRadiusTrussProblem(csvPath, modelChoice, numVariables, numberOfHeuristicObjectives, numberOfHeuristicConstraints, printableRadius, printableSideLength, printableModulus, sideNodeNumber, nucFactor, targetStiffnessRatio, engine, partialCollapsibilityConstrained[0], nodalPropertiesConstrained[0], orientationConstrained[0]);
                    globalNodePositions = ((ConstantRadiusTrussProblem) trussProblem).getNodalConnectivityArray();
                }
            }

            //String fileSaveNameInit;
            if (partialCollapsibilityConstrained[2] || nodalPropertiesConstrained[2] || orientationConstrained[2] || intersectionConstrained[2] || useLowMemberBiasing) { // Initialization object
                initialization = new BiasedInitialization(trussProblem, arteryProblem, popSize, useLowMemberBiasing, partialCollapsibilityConstrained[2], nodalPropertiesConstrained[2], orientationConstrained[2], intersectionConstrained[2], globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
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
            Variation addDiagonalMember;
            Variation addMember;
            Variation improveOrientation;
            Variation removeIntersection;

            addMember = new CompoundVariation(new AddMember(maintainFeasibility, arteryProblem, engine, globalNodePositions, sideNodeNumber, printableSideLength, numberOfHeuristicObjectives, numberOfHeuristicConstraints), new BitFlip(mutationProbability));
            //removeIntersection = new CompoundVariation(new RemoveIntersection(maintainStability, arteryProblem, engine, globalNodePositions, sideNodeNumber, printableSideLength, numberOfHeuristicObjectives, numberOfHeuristicConstraints), new BitFlip(mutationProbability));
            removeIntersection = new CompoundVariation(new RemoveIntersection2(arteryProblem, engine, globalNodePositions, sideNodeNumber, printableSideLength, numberOfHeuristicObjectives, numberOfHeuristicConstraints), new BitFlip(mutationProbability));
            addDiagonalMember = new CompoundVariation(new AddDiagonalMember(maintainFeasibility, arteryProblem, engine, globalNodePositions, sideNodeNumber, printableSideLength, numberOfHeuristicObjectives, numberOfHeuristicConstraints), new BitFlip(mutationProbability));
            //improveOrientation = new CompoundVariation(new ImproveOrientation(arteryProblem, globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints), new BitFlip(mutationProbability));
            improveOrientation = new CompoundVariation(new ImproveOrientation2(arteryProblem, globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints), new BitFlip(mutationProbability));

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
            } else if (partialCollapsibilityConstrained[6] || nodalPropertiesConstrained[6] || orientationConstrained[6] || intersectionConstrained[6]) { // Adaptive Heuristic Handling objects
                int numHeurs = 0;

                // Setup for saving results
                properties.setBoolean("saveQuality", true);
                properties.setBoolean("saveCredits", true);

                ArrayList<String> enforcedHeuristicAttributes = new ArrayList<>();
                for (int j = 0; j < heuristicsConstrained.length; j++) {
                    if (heuristicsConstrained[j][6]) {
                        enforcedHeuristicAttributes.add(heuristicAttributes[j]);
                        numHeurs++;
                    }
                }
                int[][] expectedCorrelationParameters; // each row corresponds to a heuristic, for artery problem the columns represent {feas, conn, d_min} and for truss problem the columns represent {feas, conn, stiffrat, d_min}
                if (arteryProblem) { // Corresponding to heuristic and constraint violations
                    expectedCorrelationParameters = new int[][]{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
                } else {
                    expectedCorrelationParameters = new int[][]{{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}};
                }
                HashMap<String, int[]> expectedCorrelationParamsList = new HashMap<>();
                for (int j = 0; j < heuristicsConstrained.length; j++) {
                    expectedCorrelationParamsList.put(heuristicAttributes[j], expectedCorrelationParameters[j]);
                }

                ahs = new AdaptiveHeuristicSelection(numHeurs, nfeval, enforcedHeuristicAttributes, solutions, expectedCorrelationParamsList, 0.6, 0.03);

                comp = new ChainedComparator(new AggregateConstraintComparator(), ahs, new ParetoObjectiveComparator());

                selection = new TournamentSelection(2, comp);
            } else { // Epsilon MOEA objects
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

                ////// NEW
                ArrayList<Variation> operators = new ArrayList<>();
                for (int k = 0; k < heuristicsConstrained.length; k++) {
                    if (heuristicsConstrained[k][1]) {
                        operators.add(heuristicOperators[k]);
                    }
                }

                ////// OLD
                //if (intersectionConstrained[1])
                    //Variation removeIntersection = new RemoveIntersection(false,engine,globalNodePositions,sideNodeNumber,printableSideLength);
                    //operators.add(removeIntersection);
                //}

                //if (partialCollapsibilityConstrained[1] && !nodalPropertiesConstrained[1] && !orientationConstrained[1]) {
                    //Variation addDiagonalMember = new AddDiagonalMember(maintainFeasibility, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    //operators.add(addDiagonalMember);
                //} else if (partialCollapsibilityConstrained[1] && !nodalPropertiesConstrained[1] && orientationConstrained[1]) {
                    //Variation addDiagonalMember = new AddDiagonalMember(maintainFeasibility, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    //Variation improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber);
                    //operators.add(addDiagonalMember);
                    //operators.add(improveOrientation);
                //} else if (nodalPropertiesConstrained[1] && !partialCollapsibilityConstrained[1] && !orientationConstrained[1]) {
                    //Variation addMember = new AddMember(maintainFeasibility, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    //operators.add(addMember);
                //} else if (nodalPropertiesConstrained[1] && !partialCollapsibilityConstrained[1] && orientationConstrained[1]) {
                    //Variation addMember = new AddMember(maintainFeasibility, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    //Variation improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber);
                    //operators.add(addMember);
                    //operators.add(improveOrientation);
                //} else if (partialCollapsibilityConstrained[1] && nodalPropertiesConstrained[1] && !orientationConstrained[1]) {
                    //Variation addDiagonalMember = new AddDiagonalMember(maintainFeasibility, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    //Variation addMember = new AddMember(maintainFeasibility, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    //operators.add(addMember);
                    //operators.add(addDiagonalMember);
                //} else if (partialCollapsibilityConstrained[1] && nodalPropertiesConstrained[1] && orientationConstrained[1]) {
                    //Variation addDiagonalMember = new AddDiagonalMember(maintainFeasibility, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    //Variation addMember = new AddMember(maintainFeasibility, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    //Variation improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber);
                    //operators.add(addMember);
                    //operators.add(addDiagonalMember);
                    //operators.add(improveOrientation);
                //} else if (!partialCollapsibilityConstrained[1] && !nodalPropertiesConstrained[1] && orientationConstrained[1]) {
                    //Variation improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber);
                    //operators.add(improveOrientation);
                //}

                operators.add(new CompoundVariation(new OnePointCrossover(crossoverProbability), new BitFlip(mutationProbability)));

                properties.setDouble("pmin", 0.03);

                // Create operator selector
                //OperatorSelector operatorSelector = new AdaptivePursuit(operators, 0.8, 0.8, 0.1);
                OperatorSelector operatorSelector = new ProbabilityMatching(operators, 0.6, 0.03);

                // Create credit assignment
                //SetImprovementDominance creditAssignment = new SetImprovementDominance(archive, 1, 0);
                SetContributionDominance creditAssignment = new SetContributionDominance(archive, 1, 0);
                //OffspringParentDomination creditAssignment = new OffspringParentDomination(1.0, 0.5, 0.0, comp);

                // Create AOS
                //aosStrategy = new AOSVariationSI(operatorSelector, creditAssignment, popSize);
                aosStrategy = new AOSVariationSC(operatorSelector, creditAssignment, popSize);
                //aosStrategy = new AOSVariationOP(operatorSelector, creditAssignment, popSize);

            } else { // Epsilon MOEA objects
                properties.setDouble("crossoverProbability", crossoverProbability);
                crossOver = new OnePointCrossover(crossoverProbability);
                //crossOver = new TwoPointCrossover(crossoverProbability);
                bitMutation = new BitFlip(mutationProbability);
                var = new CompoundVariation(crossOver, bitMutation);
            }

            // Creating AOS MOEA or AHS MOEA object if needed
            if (partialCollapsibilityConstrained[1] || nodalPropertiesConstrained[1] || orientationConstrained[1] || intersectionConstrained[1]) {
                EpsilonMOEA emoea = new EpsilonMOEA(trussProblem, population, archive, selection, aosStrategy, initialization, comp);
                moeaObj = new AOSMOEA(emoea, aosStrategy, true);
            } else if (partialCollapsibilityConstrained[6] || nodalPropertiesConstrained[6] || orientationConstrained[6] || intersectionConstrained[6]) {
                EpsilonMOEA emoea = new EpsilonMOEA(trussProblem, population, archive, selection, var, initialization, comp);
                moeaObj = new AHSMOEA(emoea, true, ahs);
            } else {
                moeaObj = new EpsilonMOEA(trussProblem, population, archive, selection, var, initialization, comp);
            }

            ecs.submit(new EvolutionarySearch(moeaObj, nfeval, properties, csvPath + File.separator + "result", "emoea_" + String.valueOf(i) + fileSaveNameConstraint.toString() + fileSaveNameProblem + fileSaveNameModel, sideNodeNumber, solutions, numberOfHeuristicObjectives, numberOfHeuristicConstraints));
        }

        for (int i = 0; i < numRuns; i++) {
            try {
                Algorithm alg = ecs.take().get();
            } catch (InterruptedException | ExecutionException ex) {
                Logger.getLogger(ConstantRadiusMOEARun.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        pool.shutdown();
        engine.close();

    }

}



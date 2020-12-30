package seakers.trussaos;

import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.CompoundVariation;
import org.moeaframework.core.operator.OnePointCrossover;
import org.moeaframework.core.operator.binary.BitFlip;
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
import seakers.trussaos.operators.ImproveOrientation;
import seakers.trussaos.operators.RemoveIntersection;
//import seakers.trussaos.constraints.DisjunctiveNormalForm;
import seakers.trussaos.constraints.KnowledgeStochasticRanking;
import seakers.trussaos.problems.ConstantRadiusTrussProblem;
import seakers.trussaos.problems.ConstantRadiusTrussProblem2;

/**
 * Executable class for different eMOEA run experiments in Orthogonal Arrays for the Truss Optimization problem.
 *
 * @author roshan94
 */

public class OrthogonalArrayMOEARun {

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
        boolean useOptimizationProblem2 = false; // Use ConstantRadiusTrussProblem2 as problem class (instead of ConstantRadiusTrussProblem)

        // Constraint parameters
        /**
         * feasibilityConstrained = [interior_penalty, AOS, biased_init, ACH]
         * stabilityConstrained = [interior_penalty, AOS, biased_init, ACH]
         * orientationConstrained = [interior_penalty, AOS, biased_init, ACH]
         */
        boolean[] feasibilityConstrained = {true, false, false, false};
        boolean[] stabilityConstrained = {false, false, false, false};
        boolean[] orientationConstrained = {false, false, true, false};

        int numCPU = 4;
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
        boolean maintainStability = false;
        boolean maintainFeasibility = false;

        String fileSaveNameModel;
        if (useFibreStiffness) {
            fileSaveNameModel = "_fibre";
        } else {
            fileSaveNameModel = "_truss";
        }

        String fileSaveNameProblem  = "";
        if (useOptimizationProblem2) {
            fileSaveNameProblem = "_prob2";
        }

        //String fileSaveNamePenalty;
        //if (feasibilityConstrained[0] || stabilityConstrained[0]) {
            //fileSaveNamePenalty = "penalized";
        //} else {
            //fileSaveNamePenalty = "true";
        //}

        StringBuilder fileSaveNameConstraint = new StringBuilder();
        for (int i = 0; i < 4; i++) {
            if (feasibilityConstrained[i] && !stabilityConstrained[i] && !orientationConstrained[i]) {
                fileSaveNameConstraint.append("fcon" + Integer.toString(i) + "_");
            } else if (feasibilityConstrained[i] && !stabilityConstrained[i] && orientationConstrained[i]) {
                fileSaveNameConstraint.append("focon" + Integer.toString(i) + "_");
            } else if (stabilityConstrained[i] && !feasibilityConstrained[i] && !orientationConstrained[i]) {
                fileSaveNameConstraint.append("scon" + Integer.toString(i) + "_");
            } else if (stabilityConstrained[i] && !feasibilityConstrained[i] && orientationConstrained[i]) {
                fileSaveNameConstraint.append("socon" + Integer.toString(i) + "_");
            } else if (feasibilityConstrained[i] && stabilityConstrained[i] && !orientationConstrained[i]) {
                fileSaveNameConstraint.append("fscon" + Integer.toString(i) + "_");
            } else if (feasibilityConstrained[i] && stabilityConstrained[i] && orientationConstrained[i]) {
                fileSaveNameConstraint.append("fsocon" + Integer.toString(i) + "_");
            } else if (!feasibilityConstrained[i] && !stabilityConstrained[i] && orientationConstrained[i]) {
                fileSaveNameConstraint.append("ocon" + Integer.toString(i) + "_");
            }
        }

        // New dimensions for printable solutions
        double printableRadius = 250e-6; // in m
        double printableSideLength = 10e-3; // in m
        double printableModulus = 1.8162e6; // in Pa
        int sideNodeNumber = 3;
        int nucFactor = 1;

        for (int i = 0; i < numRuns; i++) {
            // Create a new problem class
            //ConstantRadiusTrussProblem trussProblem = new ConstantRadiusTrussProblem(csvPath, useFibreStiffness, targetStiffnessRatio, engine, feasibilityConstrained[0], stabilityConstrained[0]);

            // Problem class for printable solutions
            AbstractProblem trussProblem;
            double[][] globalNodePositions;
            if (useOptimizationProblem2) {
                trussProblem = new ConstantRadiusTrussProblem2(csvPath, useFibreStiffness, 32, printableRadius, printableSideLength, printableModulus, sideNodeNumber, nucFactor, targetStiffnessRatio, engine, feasibilityConstrained[0], stabilityConstrained[0], orientationConstrained[0]);
                globalNodePositions = ((ConstantRadiusTrussProblem2) trussProblem).getNodalConnectivityArray();
            } else {
                trussProblem = new ConstantRadiusTrussProblem(csvPath, useFibreStiffness, 32, printableRadius, printableSideLength, printableModulus, sideNodeNumber, nucFactor, targetStiffnessRatio, engine, feasibilityConstrained[0], stabilityConstrained[0], orientationConstrained[0]);
                globalNodePositions = ((ConstantRadiusTrussProblem) trussProblem).getNodalConnectivityArray();
            }

            //String fileSaveNameInit;
            if (feasibilityConstrained[2] || stabilityConstrained[2]) { // Initialization object
                initialization = new BiasedInitialization(trussProblem, popSize, feasibilityConstrained[2], stabilityConstrained[2], orientationConstrained[2], globalNodePositions, targetStiffnessRatio, sideNodeNumber);
                //fileSaveNameInit = "biased_";
            } else {
                initialization = new RandomInitialization(trussProblem, popSize);
                //fileSaveNameInit = "";
            }

            // Initialize population structure for algorithm
            Population population = new Population();

            EpsilonBoxDominanceArchive archive = new EpsilonBoxDominanceArchive(epsilonDouble);

            Algorithm moeaObj;
            ChainedComparator comp = null;
            TournamentSelection selection;
            CompoundVariation var = null;
            AOSVariation aosStrategy = null;

            //moeaObj = new EpsilonMOEA(trussProblem, population, archive, selection, var, initialization);

            if (feasibilityConstrained[3] || stabilityConstrained[3] || orientationConstrained[3]) { // Adaptive Constraint Handling objects
                RemoveIntersection removeIntersection;
                AddMember addMember;
                ImproveOrientation improveOrientation;

                addMember = new AddMember(maintainFeasibility, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                removeIntersection = new RemoveIntersection(maintainStability, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, sideNodeNumber);

                //HashSet<String> constraints;
                KnowledgeStochasticRanking ksr;

                HashMap<Variation, String> constraintOperatorMap = new HashMap<>();
                if (feasibilityConstrained[3] && !stabilityConstrained[3] && !orientationConstrained[3]) {
                    constraintOperatorMap.put(removeIntersection, "FeasibilityViolation");

                    //constraints = new HashSet<>(constraintOperatorMap.values());

                    ksr = new KnowledgeStochasticRanking(constraintOperatorMap.size(), constraintOperatorMap.values(), archive);
                    comp = new ChainedComparator(ksr, new ParetoObjectiveComparator());
                } else if (feasibilityConstrained[3] && !stabilityConstrained[3] && orientationConstrained[3]) {
                    constraintOperatorMap.put(removeIntersection, "FeasibilityViolation");
                    constraintOperatorMap.put(improveOrientation, "OrientationViolation");

                   //constraints = new HashSet<>(constraintOperatorMap.values());

                    ksr = new KnowledgeStochasticRanking(constraintOperatorMap.size(), constraintOperatorMap.values(), archive);
                    comp = new ChainedComparator(ksr, new ParetoObjectiveComparator());
                } else if (stabilityConstrained[3] && !feasibilityConstrained[3] && !orientationConstrained[3]) {
                    constraintOperatorMap.put(addMember, "StabilityViolation");

                    //constraints = new HashSet<>(constraintOperatorMap.values());

                    ksr = new KnowledgeStochasticRanking(constraintOperatorMap.size(), constraintOperatorMap.values(), archive);
                    comp = new ChainedComparator(ksr, new ParetoObjectiveComparator());
                } else if (stabilityConstrained[3] && !feasibilityConstrained[3] && orientationConstrained[3]) {
                    constraintOperatorMap.put(addMember, "StabilityViolation");
                    constraintOperatorMap.put(improveOrientation, "OrientationViolation");

                    //constraints = new HashSet<>(constraintOperatorMap.values());

                    ksr = new KnowledgeStochasticRanking(constraintOperatorMap.size(), constraintOperatorMap.values(), archive);
                    comp = new ChainedComparator(ksr, new ParetoObjectiveComparator());
                } else if (feasibilityConstrained[3] && stabilityConstrained[3] && !orientationConstrained[3]) {
                    constraintOperatorMap.put(removeIntersection, "FeasibilityViolation");
                    constraintOperatorMap.put(addMember, "StabilityViolation");

                    //constraints = new HashSet<>(constraintOperatorMap.values());

                    ksr = new KnowledgeStochasticRanking(constraintOperatorMap.size(), constraintOperatorMap.values(), archive);
                    comp = new ChainedComparator(ksr, new ParetoObjectiveComparator());
                } else if (feasibilityConstrained[3] && stabilityConstrained[3] && orientationConstrained[3]) {
                    constraintOperatorMap.put(removeIntersection, "FeasibilityViolation");
                    constraintOperatorMap.put(addMember, "StabilityViolation");
                    constraintOperatorMap.put(improveOrientation, "OrientationViolation");

                    //constraints = new HashSet<>(constraintOperatorMap.values());

                    ksr = new KnowledgeStochasticRanking(constraintOperatorMap.size(), constraintOperatorMap.values(), archive);
                    comp = new ChainedComparator(ksr, new ParetoObjectiveComparator());
                } else if (!feasibilityConstrained[3] && !stabilityConstrained[3] && orientationConstrained[3]) {
                    constraintOperatorMap.put(improveOrientation, "OrientationViolation");

                    //constraints = new HashSet<>(constraintOperatorMap.values());

                    ksr = new KnowledgeStochasticRanking(constraintOperatorMap.size(), constraintOperatorMap.values(), archive);
                    comp = new ChainedComparator(ksr, new ParetoObjectiveComparator());
                }

                selection = new TournamentSelection(2, comp);

                //moeaObj = new EpsilonMOEA(trussProblem, population, archive, selection, var, initialization, comp);
            }
            else { // Epsilon MOEA objects
                comp = new ChainedComparator(new ParetoObjectiveComparator());
                selection = new TournamentSelection(2, comp);
            }

            if (feasibilityConstrained[1] || stabilityConstrained[1] || orientationConstrained[1]) { // AOS objects
                //comp = new ChainedComparator(new ParetoObjectiveComparator());
                //selection = new TournamentSelection(2, comp);

                // Setup for saving results
                properties.setBoolean("saveQuality", true);
                properties.setBoolean("saveCredits", true);
                properties.setBoolean("saveSelection", true);

                // IMPLEMENTATION WITH ACTUAL REPAIR OPERATORS
                ArrayList<Variation> operators = new ArrayList<>();

                if (feasibilityConstrained[1] && !stabilityConstrained[1] && !orientationConstrained[1]) {
                    Variation removeIntersection = new RemoveIntersection(maintainStability, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    operators.add(removeIntersection);
                } else if (feasibilityConstrained[1] && !stabilityConstrained[1] && orientationConstrained[1]) {
                    Variation removeIntersection = new RemoveIntersection(maintainStability, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    Variation improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, sideNodeNumber);
                    operators.add(removeIntersection);
                    operators.add(improveOrientation);
                } else if (stabilityConstrained[1] && !feasibilityConstrained[1] && !orientationConstrained[1]) {
                    Variation addMember = new AddMember(maintainFeasibility, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    operators.add(addMember);
                } else if (stabilityConstrained[1] && !feasibilityConstrained[1] && orientationConstrained[1]) {
                    Variation addMember = new AddMember(maintainFeasibility, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    Variation improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, sideNodeNumber);
                    operators.add(addMember);
                    operators.add(improveOrientation);
                } else if (feasibilityConstrained[1] && stabilityConstrained[1] && !orientationConstrained[1]) {
                    Variation removeIntersection = new RemoveIntersection(maintainStability, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    Variation addMember = new AddMember(maintainFeasibility, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    operators.add(addMember);
                    operators.add(removeIntersection);
                } else if (feasibilityConstrained[1] && stabilityConstrained[1] && orientationConstrained[1]) {
                    Variation removeIntersection = new RemoveIntersection(maintainStability, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    Variation addMember = new AddMember(maintainFeasibility, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    Variation improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, sideNodeNumber);
                    operators.add(addMember);
                    operators.add(removeIntersection);
                    operators.add(improveOrientation);
                } else if (!feasibilityConstrained[1] && !stabilityConstrained[1] && orientationConstrained[1]) {
                    Variation improveOrientation = new ImproveOrientation(globalNodePositions, targetStiffnessRatio, sideNodeNumber);
                    operators.add(improveOrientation);
                }

                operators.add(new CompoundVariation(new OnePointCrossover(crossoverProbability), new BitFlip(mutationProbability)));

                properties.setDouble("pmin", 0.1);

                // Create operator selector
                //OperatorSelector operatorSelector = new AdaptivePursuit(operators, 0.8, 0.8, 0.1);
                OperatorSelector operatorSelector = new ProbabilityMatching(operators, 0.6, 0.1);

                // Create credit assignment
                //SetImprovementDominance creditAssignment = new SetImprovementDominance(archive, 1, 0);
                OffspringParentDomination creditAssignment = new OffspringParentDomination(1.0, 0.5, 0.0);

                // Create AOS
                //aosStrategy = new AOSVariationSI(operatorSelector, creditAssignment, popSize);
                aosStrategy = new AOSVariationOP(operatorSelector, creditAssignment, popSize);

            }
            else { // Epsilon MOEA objects
                properties.setDouble("crossoverProbability", crossoverProbability);
                singlecross = new OnePointCrossover(crossoverProbability);
                bitFlip = new BitFlip(mutationProbability);
                var = new CompoundVariation(singlecross, bitFlip);
            }

            // Creating AOS MOEA object if needed
            if (feasibilityConstrained[1] || stabilityConstrained[1] || orientationConstrained[1]) {
                EpsilonMOEA emoea = new EpsilonMOEA(trussProblem, population, archive, selection, aosStrategy, initialization, comp);
                moeaObj = new AOSMOEA(emoea, aosStrategy, true);
            }
            else {
                moeaObj = new EpsilonMOEA(trussProblem, population, archive, selection, var, initialization, comp);
            }

            ecs.submit(new EvolutionarySearch(moeaObj, properties, csvPath + File.separator + "result", "emoea_" + String.valueOf(i) + fileSaveNameConstraint.toString() + fileSaveNameProblem + fileSaveNameModel , engine, useFibreStiffness, targetStiffnessRatio));
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



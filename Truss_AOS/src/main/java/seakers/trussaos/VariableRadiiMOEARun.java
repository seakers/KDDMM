package seakers.trussaos;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.core.*;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.ParetoObjectiveComparator;
import org.moeaframework.core.operator.*;
import org.moeaframework.core.operator.real.UM;
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
//import org.moeaframework.core.operator.real.PM;
//import org.moeaframework.core.operator.real.SBX;
import seakers.aos.operatorselectors.OperatorSelector;
import seakers.aos.operatorselectors.ProbabilityMatching;
import seakers.trussaos.initialization.BiasedInitializationVariableRadii;
//import seakers.trussaos.constrainthandling.DisjunctiveNormalForm;
import seakers.trussaos.constrainthandling.KnowledgeStochasticRanking;
import seakers.trussaos.operators.variableradii.AddDiagonalMemberVariableRadii;
import seakers.trussaos.operators.variableradii.AddMemberVariableRadii;
import seakers.trussaos.operators.variableradii.ImproveOrientationVariableRadii;
import seakers.trussaos.operators.variableradii.RemoveIntersectionVariableRadii;
import seakers.trussaos.problems.VariableRadiusTrussProblem;

/**
 * Executable class for different eMOEA run experiments in Orthogonal Arrays for the Truss Optimization problem.
 *
 * @author roshan94
 */

public class VariableRadiiMOEARun {

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

        // Constraint parameters
        /**
         * partialCollapsibilityConstrained = [interior_penalty, AOS, biased_init, ACH]
         * nodalPropertiesConstrained = [interior_penalty, AOS, biased_init, ACH]
         * orientationConstrained = [interior_penalty, AOS, biased_init, ACH]
         * feasibilityConstrained = [biased_init, AOS]
         */
        boolean[] partialCollapsibilityConstrained = {false, false, false, false};
        boolean[] nodalPropertiesConstrained = {false, false, false, false};
        boolean[] orientationConstrained = {false, true, false, false};
        boolean[] feasibilityConstrained = {true, true};

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
        int maxEvals = 10000;
        properties.setInt("maxEvaluations", maxEvals);
        properties.setInt("populationSize", popSize);

        Variation uniformCross;
        Variation uniformMutation;
        Initialization initialization;

        double crossoverProbability = 1.0;

        // For AOS MOEA Run
        boolean maintainStability = false;
        boolean maintainFeasibility = false;

        String fileSaveNameModel;
        if (useFibreStiffness) {
            fileSaveNameModel = "_fibre_varrad";
        } else {
            fileSaveNameModel = "_truss_varrad";
        }

        //String fileSaveNamePenalty;
        //if (feasibilityConstrained[0] || stabilityConstrained[0]) {
        //fileSaveNamePenalty = "penalized";
        //} else {
        //fileSaveNamePenalty = "true";
        //}

        StringBuilder fileSaveNameConstraint = new StringBuilder();
        for (int i = 0; i < 4; i++) {
            if (partialCollapsibilityConstrained[i] && !nodalPropertiesConstrained[i] && !orientationConstrained[i]) {
                fileSaveNameConstraint.append("pcon" + Integer.toString(i) + "_");
            } else if (partialCollapsibilityConstrained[i] && !nodalPropertiesConstrained[i] && orientationConstrained[i]) {
                fileSaveNameConstraint.append("pocon" + Integer.toString(i) + "_");
            } else if (nodalPropertiesConstrained[i] && !partialCollapsibilityConstrained[i] && !orientationConstrained[i]) {
                fileSaveNameConstraint.append("ncon" + Integer.toString(i) + "_");
            } else if (nodalPropertiesConstrained[i] && !partialCollapsibilityConstrained[i] && orientationConstrained[i]) {
                fileSaveNameConstraint.append("nocon" + Integer.toString(i) + "_");
            } else if (partialCollapsibilityConstrained[i] && nodalPropertiesConstrained[i] && !orientationConstrained[i]) {
                fileSaveNameConstraint.append("pncon" + Integer.toString(i) + "_");
            } else if (partialCollapsibilityConstrained[i] && nodalPropertiesConstrained[i] && orientationConstrained[i]) {
                fileSaveNameConstraint.append("pnocon" + Integer.toString(i) + "_");
            } else if (!partialCollapsibilityConstrained[i] && !nodalPropertiesConstrained[i] && orientationConstrained[i]) {
                fileSaveNameConstraint.append("ocon" + Integer.toString(i) + "_");
            }
        }
        if (feasibilityConstrained[0] && !feasibilityConstrained[1]) {
            fileSaveNameConstraint.append("feasbias_");
        } else if (!feasibilityConstrained[0] && feasibilityConstrained[1]) {
            fileSaveNameConstraint.append("feasaos_");
        } else if (feasibilityConstrained[0] && feasibilityConstrained[1]) {
            fileSaveNameConstraint.append("feasbiasaos_");
        }

        // New dimensions for printable solutions
        double printableRadius = 250e-6; // in m
        double printableSideLength = 10e-3; // in m
        double printableModulus = 1.8162e6; // in Pa
        double smallestAcceptableRadiusFactor = 0.5;
        int sideNodeNumber = 3;
        int nucFactor = 3;

        int totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sideNodeNumber*sideNodeNumber))/(CombinatoricsUtils.factorial((int) ((sideNodeNumber*sideNodeNumber) - 2)) * CombinatoricsUtils.factorial(2)));
        int numberOfRepeatableMembers = (int) (2 * (CombinatoricsUtils.factorial((int) sideNodeNumber)/(CombinatoricsUtils.factorial((int) (sideNodeNumber - 2)) * CombinatoricsUtils.factorial(2))));
        int numVariables = totalNumberOfMembers - numberOfRepeatableMembers;

        double[] printableRadiusLowerBounds = new double[totalNumberOfMembers];
        double[] printableRadiusUpperBounds = new double[totalNumberOfMembers];
        for (int i = 0; i < totalNumberOfMembers; i++) {
            printableRadiusLowerBounds[i] = 0;
            printableRadiusUpperBounds[i] = printableRadius;
        }

        double mutationProbability = 1. / numVariables;
        properties.setDouble("mutationProbability", mutationProbability);

        for (int i = 0; i < numRuns; i++) {
            // Create a new problem class
            //ConstantRadiusTrussProblem trussProblem = new ConstantRadiusTrussProblem(csvPath, useFibreStiffness, targetStiffnessRatio, engine, feasibilityConstrained[0], stabilityConstrained[0]);

            // Problem class for printable solutions
            //AbstractProblem trussProblem;
            double[][] globalNodePositions;
            VariableRadiusTrussProblem trussProblem = new VariableRadiusTrussProblem(printableRadiusLowerBounds,printableRadiusUpperBounds,csvPath,useFibreStiffness,numVariables,printableSideLength,printableModulus,nucFactor,sideNodeNumber,targetStiffnessRatio,smallestAcceptableRadiusFactor,engine,partialCollapsibilityConstrained[0], nodalPropertiesConstrained[0], orientationConstrained[0]);
            globalNodePositions = trussProblem.getNodalConnectivityArrayVariableRadii();

            //String fileSaveNameInit;
            if (partialCollapsibilityConstrained[2] || nodalPropertiesConstrained[2] || orientationConstrained[2] || feasibilityConstrained[0]) { // Initialization object
                initialization = new BiasedInitializationVariableRadii(trussProblem, popSize, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, partialCollapsibilityConstrained[2], nodalPropertiesConstrained[2], orientationConstrained[2], feasibilityConstrained[0], globalNodePositions, targetStiffnessRatio, sideNodeNumber);
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

            if (partialCollapsibilityConstrained[3] || nodalPropertiesConstrained[3] || orientationConstrained[3]) { // Adaptive Constraint Handling objects
                AddDiagonalMemberVariableRadii addDiagonalMember;
                AddMemberVariableRadii addMember;
                ImproveOrientationVariableRadii improveOrientation;

                addMember = new AddMemberVariableRadii(maintainFeasibility, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                //removeIntersection = new RemoveIntersection(maintainStability, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                addDiagonalMember = new AddDiagonalMemberVariableRadii(maintainFeasibility, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                improveOrientation = new ImproveOrientationVariableRadii(globalNodePositions, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, targetStiffnessRatio, sideNodeNumber);

                //HashSet<String> constrainthandling;
                KnowledgeStochasticRanking ksr;

                HashMap<Variation, String> constraintOperatorMap = new HashMap<>();
                if (partialCollapsibilityConstrained[3] && !nodalPropertiesConstrained[3] && !orientationConstrained[3]) {
                    //constraintOperatorMap.put(removeIntersection, "FeasibilityViolation");
                    constraintOperatorMap.put(addDiagonalMember, "PartialCollapsibilityViolation");

                    //constrainthandling = new HashSet<>(constraintOperatorMap.values());
                } else if (partialCollapsibilityConstrained[3] && !nodalPropertiesConstrained[3] && orientationConstrained[3]) {
                    constraintOperatorMap.put(addDiagonalMember, "PartialCollapsibilityViolation");
                    constraintOperatorMap.put(improveOrientation, "OrientationViolation");

                    //constrainthandling = new HashSet<>(constraintOperatorMap.values());
                } else if (nodalPropertiesConstrained[3] && !partialCollapsibilityConstrained[3] && !orientationConstrained[3]) {
                    constraintOperatorMap.put(addMember, "NodalPropertiesViolation");

                    //constrainthandling = new HashSet<>(constraintOperatorMap.values());
                } else if (nodalPropertiesConstrained[3] && !partialCollapsibilityConstrained[3] && orientationConstrained[3]) {
                    constraintOperatorMap.put(addMember, "NodalPropertiesViolation");
                    constraintOperatorMap.put(improveOrientation, "OrientationViolation");

                    //constrainthandling = new HashSet<>(constraintOperatorMap.values());
                } else if (partialCollapsibilityConstrained[3] && nodalPropertiesConstrained[3] && !orientationConstrained[3]) {
                    constraintOperatorMap.put(addDiagonalMember, "PartialCollapsibilityViolation");
                    constraintOperatorMap.put(addMember, "StabilityViolation");

                    //constrainthandling = new HashSet<>(constraintOperatorMap.values());
                } else if (partialCollapsibilityConstrained[3] && nodalPropertiesConstrained[3] && orientationConstrained[3]) {
                    constraintOperatorMap.put(addDiagonalMember, "PartialCollapsibilityViolation");
                    constraintOperatorMap.put(addMember, "StabilityViolation");
                    constraintOperatorMap.put(improveOrientation, "OrientationViolation");

                    //constrainthandling = new HashSet<>(constraintOperatorMap.values());
                } else if (!partialCollapsibilityConstrained[3] && !nodalPropertiesConstrained[3] && orientationConstrained[3]) {
                    constraintOperatorMap.put(improveOrientation, "OrientationViolation");

                    //constrainthandling = new HashSet<>(constraintOperatorMap.values());
                }
                ksr = new KnowledgeStochasticRanking(constraintOperatorMap.size(), constraintOperatorMap.values(), archive);
                comp = new ChainedComparator(ksr, new ParetoObjectiveComparator());

                selection = new TournamentSelection(2, comp);

                //moeaObj = new EpsilonMOEA(trussProblem, population, archive, selection, var, initialization, comp);
            }
            else { // Epsilon MOEA objects
                comp = new ChainedComparator(new ParetoObjectiveComparator());
                selection = new TournamentSelection(2, comp);
            }

            if (partialCollapsibilityConstrained[1] || nodalPropertiesConstrained[1] || orientationConstrained[1]) { // AOS objects
                //comp = new ChainedComparator(new ParetoObjectiveComparator());
                //selection = new TournamentSelection(2, comp);

                // Setup for saving results
                properties.setBoolean("saveQuality", true);
                properties.setBoolean("saveCredits", true);
                properties.setBoolean("saveSelection", true);

                // IMPLEMENTATION WITH ACTUAL REPAIR OPERATORS
                ArrayList<Variation> operators = new ArrayList<>();
                if (feasibilityConstrained[1]) {
                    Variation removeIntersection = new RemoveIntersectionVariableRadii(false, engine, globalNodePositions, sideNodeNumber, printableSideLength, printableRadiusLowerBounds,  printableRadiusUpperBounds, smallestAcceptableRadiusFactor);
                    operators.add(removeIntersection);
                }

                if (partialCollapsibilityConstrained[1] && !nodalPropertiesConstrained[1] && !orientationConstrained[1]) {
                    Variation addDiagonalMember = new AddDiagonalMemberVariableRadii(maintainFeasibility, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    operators.add(addDiagonalMember);
                } else if (partialCollapsibilityConstrained[1] && !nodalPropertiesConstrained[1] && orientationConstrained[1]) {
                    Variation addDiagonalMember = new AddDiagonalMemberVariableRadii(maintainFeasibility, printableRadiusLowerBounds,  printableRadiusUpperBounds, smallestAcceptableRadiusFactor, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    Variation improveOrientation = new ImproveOrientationVariableRadii(globalNodePositions, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, targetStiffnessRatio, sideNodeNumber);
                    operators.add(addDiagonalMember);
                    operators.add(improveOrientation);
                } else if (nodalPropertiesConstrained[1] && !partialCollapsibilityConstrained[1] && !orientationConstrained[1]) {
                    Variation addMember = new AddMemberVariableRadii(maintainFeasibility, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    operators.add(addMember);
                } else if (nodalPropertiesConstrained[1] && !partialCollapsibilityConstrained[1] && orientationConstrained[1]) {
                    Variation addMember = new AddMemberVariableRadii(maintainFeasibility, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    Variation improveOrientation = new ImproveOrientationVariableRadii(globalNodePositions, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, targetStiffnessRatio, sideNodeNumber);
                    operators.add(addMember);
                    operators.add(improveOrientation);
                } else if (partialCollapsibilityConstrained[1] && nodalPropertiesConstrained[1] && !orientationConstrained[1]) {
                    Variation addDiagonalMember = new AddDiagonalMemberVariableRadii(maintainFeasibility, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    Variation addMember = new AddMemberVariableRadii(maintainFeasibility, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor,engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    operators.add(addMember);
                    operators.add(addDiagonalMember);
                } else if (partialCollapsibilityConstrained[1] && nodalPropertiesConstrained[1] && orientationConstrained[1]) {
                    Variation addDiagonalMember = new AddDiagonalMemberVariableRadii(maintainFeasibility, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    Variation addMember = new AddMemberVariableRadii(maintainFeasibility, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, engine, globalNodePositions, sideNodeNumber, printableSideLength);
                    Variation improveOrientation = new ImproveOrientationVariableRadii(globalNodePositions, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, targetStiffnessRatio, sideNodeNumber);
                    operators.add(addMember);
                    operators.add(addDiagonalMember);
                    operators.add(improveOrientation);
                } else if (!partialCollapsibilityConstrained[1] && !nodalPropertiesConstrained[1] && orientationConstrained[1]) {
                    Variation improveOrientation = new ImproveOrientationVariableRadii(globalNodePositions, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor, targetStiffnessRatio, sideNodeNumber);
                    operators.add(improveOrientation);
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
                uniformCross = new UniformCrossover(crossoverProbability);
                uniformMutation = new UM(mutationProbability);
                var = new CompoundVariation(uniformCross, uniformMutation);
            }

            // Creating AOS MOEA object if needed
            if (partialCollapsibilityConstrained[1] || nodalPropertiesConstrained[1] || orientationConstrained[1] || feasibilityConstrained[1]) {
                EpsilonMOEA emoea = new EpsilonMOEA(trussProblem, population, archive, selection, aosStrategy, initialization, comp);
                moeaObj = new AOSMOEA(emoea, aosStrategy, true);
            }
            else {
                moeaObj = new EpsilonMOEA(trussProblem, population, archive, selection, var, initialization, comp);
            }

            ecs.submit(new VariableRadiiEvolutionarySearch(moeaObj, properties, csvPath + File.separator + "result", "emoea_" + String.valueOf(i) + fileSaveNameConstraint.toString() + fileSaveNameModel , engine, useFibreStiffness, targetStiffnessRatio, sideNodeNumber, printableRadiusLowerBounds, printableRadiusUpperBounds, smallestAcceptableRadiusFactor));
        }

        for (int i = 0; i < numRuns; i++) {
            try {
                Algorithm alg = ecs.take().get();
            } catch (InterruptedException | ExecutionException ex) {
                Logger.getLogger(VariableRadiiMOEARun.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        pool.shutdown();
        engine.close();

    }

}



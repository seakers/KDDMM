package seakers.trussaos;

import com.mathworks.engine.EngineException;
import com.mathworks.engine.MatlabEngine;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.core.operator.CompoundVariation;
import org.moeaframework.core.operator.RandomInitialization;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.problem.AbstractProblem;
import org.moeaframework.util.TypedProperties;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;
import seakers.trussaos.initialization.SynchronizedMersenneTwister;
import seakers.trussaos.operators.constantradii.AddDiagonalMember;
import seakers.trussaos.operators.constantradii.AddMember;
import seakers.trussaos.operators.constantradii.ImproveOrientation2;
import seakers.trussaos.operators.constantradii.RemoveIntersection2;
import seakers.trussaos.problems.ConstantRadiusArteryProblem;
import seakers.trussaos.problems.ConstantRadiusTrussProblem;
import seakers.trussaos.problems.ConstantRadiusTrussProblem2;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.ExecutionException;

public class GenerateOperatorIndexDataset {

    /**
     * Matlab Engine for function evaluation
     */
    private static MatlabEngine engine;

    public static void main(String[] args) throws InterruptedException, ExecutionException, EngineException, IOException {

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
        boolean[] partialCollapsibilityConstrained = {false, false, false, false, false, false, false};
        boolean[] nodalPropertiesConstrained = {false, false, false, false, false, false, false};
        boolean[] orientationConstrained = {false, true, false, false, false, false, false};
        boolean[] intersectionConstrained = {false, true, false, false, false, false, false};

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

        engine = MatlabEngine.startMatlab();

        // Create the desired algorithm parameters and operators for search
        TypedProperties properties = new TypedProperties();

        Variation bitMutation;
        Initialization initialization;

        // For AOS MOEA Run
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

        int numRuns = 10;
        int numDesigns = 1000;

        for (int i = 0; i < numRuns; i++) {

            System.out.println("Generating data for run " + i);

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

            initialization = new RandomInitialization(trussProblem, numDesigns);
            Solution[] population = initialization.initialize();

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

            String csvFilename = csvPath + File.separator + "result" + File.separator + "random_operator_index_data" + fileSaveNameProblem + fileSaveNameModel + i + ".csv";
            File csvFile = new File(csvFilename);
            FileWriter csvWrite = new FileWriter(csvFile);

            // CSV file headers
            csvWrite.append("Full Design - Initial");
            csvWrite.append(",");
            csvWrite.append("Pen. Obj. 1");
            csvWrite.append(",");
            csvWrite.append("Pen. Obj. 2");
            csvWrite.append(",");
            csvWrite.append("Feasibility");
            csvWrite.append(",");
            csvWrite.append("Connectivity");
            csvWrite.append(",");
            if (!arteryProblem) {
                csvWrite.append("Stiffness Ratio");
                csvWrite.append(",");
            }
            csvWrite.append("Full Design - Partial Collapsibility");
            csvWrite.append(",");
            csvWrite.append("Pen. Obj. 1");
            csvWrite.append(",");
            csvWrite.append("Pen. Obj. 2");
            csvWrite.append(",");
            csvWrite.append("Feasibility");
            csvWrite.append(",");
            csvWrite.append("Connectivity");
            csvWrite.append(",");
            if (!arteryProblem) {
                csvWrite.append("Stiffness Ratio");
                csvWrite.append(",");
            }
            csvWrite.append("Full Design - Nodal Properties");
            csvWrite.append(",");
            csvWrite.append("Pen. Obj. 1");
            csvWrite.append(",");
            csvWrite.append("Pen. Obj. 2");
            csvWrite.append(",");
            csvWrite.append("Feasibility");
            csvWrite.append(",");
            csvWrite.append("Connectivity");
            csvWrite.append(",");
            if (!arteryProblem) {
                csvWrite.append("Stiffness Ratio");
                csvWrite.append(",");
            }
            csvWrite.append("Full Design - Orientation");
            csvWrite.append(",");
            csvWrite.append("Pen. Obj. 1");
            csvWrite.append(",");
            csvWrite.append("Pen. Obj. 2");
            csvWrite.append(",");
            csvWrite.append("Feasibility");
            csvWrite.append(",");
            csvWrite.append("Connectivity");
            csvWrite.append(",");
            if (!arteryProblem) {
                csvWrite.append("Stiffness Ratio");
                csvWrite.append(",");
            }
            csvWrite.append("Full Design - Intersection");
            csvWrite.append(",");
            csvWrite.append("Pen. Obj. 1");
            csvWrite.append(",");
            csvWrite.append("Pen. Obj. 2");
            csvWrite.append(",");
            csvWrite.append("Feasibility");
            csvWrite.append(",");
            csvWrite.append("Connectivity");
            csvWrite.append(",");
            if (!arteryProblem) {
                csvWrite.append("Stiffness Ratio");
                csvWrite.append(",");
            }
            csvWrite.append("\n");

            for (int j = 0; j < population.length; j++) {
                // Evaluate current solution
                trussProblem.evaluate(population[j]);

                // Operate using Partial Collapsibility Repair Operator and evaluate
                Solution partialCollapsibilitySolution = addDiagonalMember.evolve(new Solution[]{population[j]})[0];
                trussProblem.evaluate(partialCollapsibilitySolution);

                // Operate using Nodal Properties Repair Operator and evaluate
                Solution nodalPropertiesSolution = addMember.evolve(new Solution[]{population[j]})[0];
                trussProblem.evaluate(nodalPropertiesSolution);

                // Operate using Orientation Repair Operator and evaluate
                Solution orientationSolution = improveOrientation.evolve(new Solution[]{population[j]})[0];
                trussProblem.evaluate(orientationSolution);

                // Operate using Intersection Repair Operator and evaluate
                Solution intersectionSolution = removeIntersection.evolve(new Solution[]{population[j]})[0];
                trussProblem.evaluate(intersectionSolution);

                // Populate row of CSV file
                TrussRepeatableArchitecture currentArch = new TrussRepeatableArchitecture(population[j], sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
                boolean[] currentBooleanDesign = currentArch.getCompleteBooleanArrayFromSolutionNxN(population[j]);
                csvWrite.append(convertBooleanArrayToBitstring(currentBooleanDesign));
                csvWrite.append(",");

                csvWrite.append(Double.toString(population[j].getObjective(0)));
                csvWrite.append(",");

                csvWrite.append(Double.toString(population[j].getObjective(1)));
                csvWrite.append(",");

                csvWrite.append(Double.toString(1.0 - (Double) population[j].getAttribute("FeasibilityViolation")));
                csvWrite.append(",");

                csvWrite.append(Double.toString(1.0 - (Double) population[j].getAttribute("ConnectivityViolation")));
                csvWrite.append(",");

                if (!arteryProblem) {
                    csvWrite.append(Double.toString((Double) population[j].getAttribute("StiffnessRatioViolation")));
                    csvWrite.append(",");
                }

                TrussRepeatableArchitecture partCollArch = new TrussRepeatableArchitecture(partialCollapsibilitySolution, sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
                boolean[] partialCollapsibilityBooleanDesign = partCollArch.getCompleteBooleanArrayFromSolutionNxN(partialCollapsibilitySolution);
                csvWrite.append(convertBooleanArrayToBitstring(partialCollapsibilityBooleanDesign));
                csvWrite.append(",");

                csvWrite.append(Double.toString(partialCollapsibilitySolution.getObjective(0)));
                csvWrite.append(",");

                csvWrite.append(Double.toString(partialCollapsibilitySolution.getObjective(1)));
                csvWrite.append(",");

                csvWrite.append(Double.toString(1.0 - (Double) partialCollapsibilitySolution.getAttribute("FeasibilityViolation")));
                csvWrite.append(",");

                csvWrite.append(Double.toString(1.0 - (Double) partialCollapsibilitySolution.getAttribute("ConnectivityViolation")));
                csvWrite.append(",");

                if (!arteryProblem) {
                    csvWrite.append(Double.toString((Double) partialCollapsibilitySolution.getAttribute("StiffnessRatioViolation")));
                    csvWrite.append(",");
                }

                TrussRepeatableArchitecture nodalPropArch = new TrussRepeatableArchitecture(nodalPropertiesSolution, sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
                boolean[] nodalPropertiesBooleanDesign = nodalPropArch.getCompleteBooleanArrayFromSolutionNxN(nodalPropertiesSolution);
                csvWrite.append(convertBooleanArrayToBitstring(nodalPropertiesBooleanDesign));
                csvWrite.append(",");

                csvWrite.append(Double.toString(nodalPropertiesSolution.getObjective(0)));
                csvWrite.append(",");

                csvWrite.append(Double.toString(nodalPropertiesSolution.getObjective(1)));
                csvWrite.append(",");

                csvWrite.append(Double.toString(1.0 - (Double) nodalPropertiesSolution.getAttribute("FeasibilityViolation")));
                csvWrite.append(",");

                csvWrite.append(Double.toString(1.0 - (Double) nodalPropertiesSolution.getAttribute("ConnectivityViolation")));
                csvWrite.append(",");

                if (!arteryProblem) {
                    csvWrite.append(Double.toString((Double) nodalPropertiesSolution.getAttribute("StiffnessRatioViolation")));
                    csvWrite.append(",");
                }

                TrussRepeatableArchitecture orientArch = new TrussRepeatableArchitecture(orientationSolution, sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
                boolean[] orientationBooleanDesign = orientArch.getCompleteBooleanArrayFromSolutionNxN(orientationSolution);
                csvWrite.append(convertBooleanArrayToBitstring(orientationBooleanDesign));
                csvWrite.append(",");

                csvWrite.append(Double.toString(orientationSolution.getObjective(0)));
                csvWrite.append(",");

                csvWrite.append(Double.toString(orientationSolution.getObjective(1)));
                csvWrite.append(",");

                csvWrite.append(Double.toString(1.0 - (Double) orientationSolution.getAttribute("FeasibilityViolation")));
                csvWrite.append(",");

                csvWrite.append(Double.toString(1.0 - (Double) orientationSolution.getAttribute("ConnectivityViolation")));
                csvWrite.append(",");

                if (!arteryProblem) {
                    csvWrite.append(Double.toString((Double) orientationSolution.getAttribute("StiffnessRatioViolation")));
                    csvWrite.append(",");
                }

                TrussRepeatableArchitecture intersArch = new TrussRepeatableArchitecture(intersectionSolution, sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
                boolean[] intersectionBooleanDesign = intersArch.getCompleteBooleanArrayFromSolutionNxN(intersectionSolution);
                csvWrite.append(convertBooleanArrayToBitstring(intersectionBooleanDesign));
                csvWrite.append(",");

                csvWrite.append(Double.toString(intersectionSolution.getObjective(0)));
                csvWrite.append(",");

                csvWrite.append(Double.toString(intersectionSolution.getObjective(1)));
                csvWrite.append(",");

                csvWrite.append(Double.toString(1.0 - (Double) intersectionSolution.getAttribute("FeasibilityViolation")));
                csvWrite.append(",");

                csvWrite.append(Double.toString(1.0 - (Double) intersectionSolution.getAttribute("ConnectivityViolation")));

                if (!arteryProblem) {
                    csvWrite.append(",");
                    csvWrite.append(Double.toString((Double) intersectionSolution.getAttribute("StiffnessRatioViolation")));
                }
                csvWrite.append("\n");

                System.out.println("Evaluated design no. " + j);
            }
            csvWrite.flush();
            csvWrite.close();
        }
        engine.close();
    }

    private static String convertBooleanArrayToBitstring (boolean[] design) {
        StringBuilder designString = new StringBuilder();
        for (boolean b : design) {
            if (b) {
                designString.append(Integer.toString(1));
            } else {
                designString.append(Integer.toString(0));
            }
        }
        return designString.toString();
    }
}

package seakers.trussaos;

import com.mathworks.engine.EngineException;
import com.mathworks.engine.MatlabEngine;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.core.*;
import org.moeaframework.core.operator.CompoundVariation;
import org.moeaframework.core.operator.RandomInitialization;
import org.moeaframework.core.operator.binary.BitFlip;
import org.moeaframework.problem.AbstractProblem;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;
import seakers.trussaos.initialization.SynchronizedMersenneTwister;
import seakers.trussaos.operators.constantradii.AddDiagonalMember;
import seakers.trussaos.operators.constantradii.AddMember;
import seakers.trussaos.operators.constantradii.ImproveOrientation2;
import seakers.trussaos.operators.constantradii.RemoveIntersection2;
import seakers.trussaos.problems.ConstantRadiusArteryProblem;
import seakers.trussaos.problems.ConstantRadiusTrussProblem2;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Checks to make sure each heuristic operator improves the heuristic satisfaction of a design. This is done by generating random designs for each run and
 * computing heuristic function value before and after manipulation by the heuristic operator.
 */

public class CheckOperatorHeuristicImprovement {

    /**
     * Matlab Engine for function evaluation
     */
    private static MatlabEngine engine;

    public static void main(String[] args) throws InterruptedException, EngineException, IOException {

        // Define problem parameters
        //String csvPath = "C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\src\\main\\java\\seakers\\trussaos";
        String csvPath = System.getProperty("user.dir");

        boolean arteryProblem = true; // True -> artery problem, False -> truss problem (Keep as Truss Problem)

        /**
         * modelChoice = 0 --> Fibre Stiffness Model
         *             = 1 --> Truss Stiffness Model
         *             = 2 --> Beam Model
         */
        int modelChoice = 1; // Fibre stiffness model cannot be used for the artery problem

        double targetStiffnessRatio = 1;
        if (arteryProblem) {
            targetStiffnessRatio = 0.421;
        }

        engine = MatlabEngine.startMatlab();

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

        // Bias initial population with low number of members
        boolean useLowMemberBiasing = false;

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

        String fileSaveNameProblem  = "";
        if (arteryProblem) {
            fileSaveNameProblem = "_artery";
        } else {
            fileSaveNameProblem = "_truss";
        }

        int numDesigns = 100;

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

        // For AOS MOEA Run
        boolean maintainFeasibility = false;

        PRNG.setRandom(new SynchronizedMersenneTwister());

        Initialization initialization;

        // Problem class
        AbstractProblem trussProblem;
        double[][] globalNodePositions;
        if (arteryProblem) {
            trussProblem = new ConstantRadiusArteryProblem(csvPath, modelChoice, numVariables, numberOfHeuristicObjectives, numberOfHeuristicConstraints, printableRadius, printableSideLength, printableModulus, sideNodeNumber, nucFactor, targetStiffnessRatio, engine, heuristicsConstrained);
            globalNodePositions = ((ConstantRadiusArteryProblem) trussProblem).getNodalConnectivityArray();
        } else {
            trussProblem = new ConstantRadiusTrussProblem2(csvPath, modelChoice, numVariables, numberOfHeuristicObjectives, numberOfHeuristicConstraints, printableRadius, printableSideLength, printableModulus, sideNodeNumber, nucFactor, targetStiffnessRatio, engine, heuristicsConstrained);
            globalNodePositions = ((ConstantRadiusTrussProblem2) trussProblem).getNodalConnectivityArray();
        }

        // Initialize heuristic operators
        //Variation addMember = new CompoundVariation(new AddMember(maintainFeasibility, arteryProblem, engine, globalNodePositions, sideNodeNumber, printableSideLength, numberOfHeuristicObjectives, numberOfHeuristicConstraints), new BitFlip(mutationProbability));
        //Variation removeIntersection = new CompoundVariation(new RemoveIntersection2(arteryProblem, engine, globalNodePositions, sideNodeNumber, printableSideLength, numberOfHeuristicObjectives, numberOfHeuristicConstraints), new BitFlip(mutationProbability));
        //Variation addDiagonalMember = new CompoundVariation(new AddDiagonalMember(maintainFeasibility, arteryProblem, engine, globalNodePositions, sideNodeNumber, printableSideLength, numberOfHeuristicObjectives, numberOfHeuristicConstraints), new BitFlip(mutationProbability));
        //Variation improveOrientation = new CompoundVariation(new ImproveOrientation2(arteryProblem, globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints), new BitFlip(mutationProbability));

        Variation addMember = new AddMember(maintainFeasibility, arteryProblem, engine, globalNodePositions, sideNodeNumber, printableSideLength, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
        Variation removeIntersection = new RemoveIntersection2(arteryProblem, engine, globalNodePositions, sideNodeNumber, printableSideLength, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
        Variation addDiagonalMember = new AddDiagonalMember(maintainFeasibility, arteryProblem, engine, globalNodePositions, sideNodeNumber, printableSideLength, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
        Variation improveOrientation = new ImproveOrientation2(arteryProblem, globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);

        String csvFilename = csvPath + File.separator + "result" + File.separator + "operator_heuristic_satisfaction" + fileSaveNameProblem + ".csv";
        File csvFile = new File(csvFilename);
        FileWriter csvWrite = new FileWriter(csvFile);

        // CSV file headers
        csvWrite.append("Full Design - Initial");
        csvWrite.append(",");
        csvWrite.append("Full Design - Partial Collapsibility");
        csvWrite.append(",");
        csvWrite.append("Partial Collapsibility - Before");
        csvWrite.append(",");
        csvWrite.append("Partial Collapsibility - After");
        csvWrite.append(",");
        csvWrite.append("Full Design - Nodal Properties");
        csvWrite.append(",");
        csvWrite.append("Nodal Properties - Before");
        csvWrite.append(",");
        csvWrite.append("Nodal Properties - After");
        csvWrite.append(",");
        csvWrite.append("Full Design - Orientation");
        csvWrite.append(",");
        csvWrite.append("Orientation - Before");
        csvWrite.append(",");
        csvWrite.append("Orientation - After");
        csvWrite.append(",");
        csvWrite.append("Full Design - Intersection");
        csvWrite.append(",");
        csvWrite.append("Intersection - Before");
        csvWrite.append(",");
        csvWrite.append("Intersection - After");
        csvWrite.append("\n");

        // Randomly generate test designs
        initialization = new RandomInitialization(trussProblem, 1);
        //Solution[] designSolutions = initialization.initialize();

        int designCounter = 0;

        while (designCounter < numDesigns) {
            Solution currentSolution = initialization.initialize()[0];

            // Evaluate current solution to obtain prior heuristic values
            trussProblem.evaluate(currentSolution);
            double currentPartialCollapsibilityScore = 1.0d - (double)currentSolution.getAttribute("PartialCollapsibilityViolation");
            double currentNodalPropertiesScore = 1.0d - (double)currentSolution.getAttribute("NodalPropertiesViolation");
            double currentOrientationScore = 1.0d - (double)currentSolution.getAttribute("OrientationViolation");
            double currentIntersectionScore = 1.0d - (double)currentSolution.getAttribute("IntersectionViolation");

            // Only accept designs that can be improved by the heuristic operator
            if ((currentPartialCollapsibilityScore == 1) || (currentNodalPropertiesScore == 1) || (currentOrientationScore == 1) || (currentIntersectionScore == 1)) {
                continue;
            } else {
                // Pass through partial collapsibility operator and evaluate new design to obtain post heuristic values
                Solution partialCollapsibilitySolution = addDiagonalMember.evolve(new Solution[]{currentSolution})[0];
                trussProblem.evaluate(partialCollapsibilitySolution);

                // Pass through nodal properties operator and evaluate new design to obtain post heuristic values
                Solution nodalPropertiesSolution = addMember.evolve(new Solution[]{currentSolution})[0];
                trussProblem.evaluate(nodalPropertiesSolution);

                // Pass through orientation operator and evaluate new design to obtain post heuristic values
                Solution orientationSolution = improveOrientation.evolve(new Solution[]{currentSolution})[0];
                trussProblem.evaluate(orientationSolution);

                // Pass through intersection operator and evaluate new design to obtain post heuristic values
                Solution intersectionSolution = removeIntersection.evolve(new Solution[]{currentSolution})[0];
                trussProblem.evaluate(intersectionSolution);

                // Populate row of CSV file
                TrussRepeatableArchitecture currentArch = new TrussRepeatableArchitecture(currentSolution, sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
                boolean[] currentBooleanDesign = currentArch.getCompleteBooleanArrayFromSolutionNxN(currentSolution);
                csvWrite.append(convertBooleanArrayToBitstring(currentBooleanDesign));
                csvWrite.append(",");

                TrussRepeatableArchitecture partCollArch = new TrussRepeatableArchitecture(partialCollapsibilitySolution, sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
                boolean[] partialCollapsibilityBooleanDesign = partCollArch.getCompleteBooleanArrayFromSolutionNxN(partialCollapsibilitySolution);
                csvWrite.append(convertBooleanArrayToBitstring(partialCollapsibilityBooleanDesign));
                csvWrite.append(",");
                double newPartialCollapsibilityScore = 1.0d - (double)partialCollapsibilitySolution.getAttribute("PartialCollapsibilityViolation");
                csvWrite.append(Double.toString(currentPartialCollapsibilityScore));
                csvWrite.append(",");
                csvWrite.append(Double.toString(newPartialCollapsibilityScore));
                csvWrite.append(",");

                TrussRepeatableArchitecture nodalPropArch = new TrussRepeatableArchitecture(nodalPropertiesSolution, sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
                boolean[] nodalPropertiesBooleanDesign = nodalPropArch.getCompleteBooleanArrayFromSolutionNxN(nodalPropertiesSolution);
                csvWrite.append(convertBooleanArrayToBitstring(nodalPropertiesBooleanDesign));
                csvWrite.append(",");
                double newNodalPropertiesScore = 1.0d - (double)nodalPropertiesSolution.getAttribute("NodalPropertiesViolation");
                csvWrite.append(Double.toString(currentNodalPropertiesScore));
                csvWrite.append(",");
                csvWrite.append(Double.toString(newNodalPropertiesScore));
                csvWrite.append(",");

                TrussRepeatableArchitecture orientArch = new TrussRepeatableArchitecture(orientationSolution, sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
                boolean[] orientationBooleanDesign = orientArch.getCompleteBooleanArrayFromSolutionNxN(orientationSolution);
                csvWrite.append(convertBooleanArrayToBitstring(orientationBooleanDesign));
                csvWrite.append(",");
                double newOrientationScore = 1.0d - (double)orientationSolution.getAttribute("OrientationViolation");
                csvWrite.append(Double.toString(currentOrientationScore));
                csvWrite.append(",");
                csvWrite.append(Double.toString(newOrientationScore));
                csvWrite.append(",");

                TrussRepeatableArchitecture intersArch = new TrussRepeatableArchitecture(intersectionSolution, sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
                boolean[] intersectionBooleanDesign = intersArch.getCompleteBooleanArrayFromSolutionNxN(intersectionSolution);
                csvWrite.append(convertBooleanArrayToBitstring(intersectionBooleanDesign));
                csvWrite.append(",");
                double newIntersectionScore = 1.0d - (double)intersectionSolution.getAttribute("IntersectionViolation");
                csvWrite.append(Double.toString(currentIntersectionScore));
                csvWrite.append(",");
                csvWrite.append(Double.toString(newIntersectionScore));
                csvWrite.append("\n");

                designCounter++;
                System.out.println("Evaluated design no. " + designCounter);
            }
        }
        csvWrite.flush();
        csvWrite.close();
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

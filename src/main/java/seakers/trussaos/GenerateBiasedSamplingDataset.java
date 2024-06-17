package seakers.trussaos;

import com.mathworks.engine.EngineException;
import com.mathworks.engine.MatlabEngine;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.core.operator.RandomInitialization;
import org.moeaframework.problem.AbstractProblem;
import org.moeaframework.util.TypedProperties;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;
import seakers.trussaos.initialization.BiasedInitialization;
import seakers.trussaos.problems.ConstantRadiusArteryProblem;
import seakers.trussaos.problems.ConstantRadiusTrussProblem2;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class GenerateBiasedSamplingDataset {

    /**
     * Matlab Engine for function evaluation
     */
    private static MatlabEngine engine;

    public static void main(String[] args) throws InterruptedException, EngineException, IOException {

        String csvPath = System.getProperty("user.dir");

        /**
         * modelChoice = 0 --> Fibre Stiffness Model
         *             = 1 --> Truss Stiffness Model
         *             = 2 --> Beam Model
         */
        int modelChoice = 1; // Fibre stiffness model cannot be used for the artery problem

        boolean arteryProblem = true; // Solve the artery optimization (otherwise the original truss problem is solved)

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
        boolean[] orientationConstrained = {false, false, true, false, false, false, false};
        boolean[] intersectionConstrained = {false, false, false, false, false, false, false};

        // Bias initial population with low number of members
        boolean useLowMemberBiasing = false;

        boolean[][] heuristicsConstrained = new boolean[4][7];
        for (int i = 0; i < 7; i++) {
            heuristicsConstrained[0][i] = partialCollapsibilityConstrained[i];
            heuristicsConstrained[1][i] = nodalPropertiesConstrained[i];
            heuristicsConstrained[2][i] = orientationConstrained[i];
            heuristicsConstrained[3][i] = intersectionConstrained[i];
        }

        int numberOfBiasedSampleHeuristics = 0;
        int numberOfHeuristicConstraints = 0;
        int numberOfHeuristicObjectives = 0;
        boolean[] biasedSamplingEnforcement = new boolean[4];
        for (int i = 0; i < 4; i++) {
            biasedSamplingEnforcement[i] = heuristicsConstrained[i][2];
            if (heuristicsConstrained[i][2]) {
                numberOfBiasedSampleHeuristics++;
            }
            if (heuristicsConstrained[i][5]) {
                numberOfHeuristicConstraints++;
            }
            if (heuristicsConstrained[i][4]) {
                numberOfHeuristicObjectives++;
            }
        }

        engine = MatlabEngine.startMatlab();

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
            fileSaveNameProblem = "_prob2";
        }

        String[] heuristicNames = new String[]{"pc", "np", "or", "is"};
        StringBuilder fileSaveNameHeuristic = new StringBuilder();
        for (int i = 0; i < biasedSamplingEnforcement.length; i++) {
            if (biasedSamplingEnforcement[i]) {
                fileSaveNameHeuristic.append("_").append(heuristicNames[i]);
            }
        }

        int numDesigns = 300;
        int numRuns = 10;

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

        Initialization initialization;
        if (numberOfBiasedSampleHeuristics == 0) {
            initialization = new RandomInitialization(trussProblem, numDesigns);
        } else {
            initialization = new BiasedInitialization(trussProblem, arteryProblem, numDesigns, useLowMemberBiasing, partialCollapsibilityConstrained[2], nodalPropertiesConstrained[2], orientationConstrained[2], intersectionConstrained[2], globalNodePositions, targetStiffnessRatio, (int) sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
        }

        String csvFilename = "";

        for (int i = 0; i < numRuns; i++) {
            System.out.println("Generating random data for run " + i);

            Solution[] population = initialization.initialize();
            csvFilename = csvPath + File.separator + "result" + File.separator + "random_biased_sampling_index_data" + fileSaveNameProblem + fileSaveNameModel + fileSaveNameHeuristic.toString() + i + ".csv";
            computeAndStoreIntoCsv(population, csvFilename, trussProblem, arteryProblem, sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
        }
        engine.close();

    }

    private static void computeAndStoreIntoCsv(Solution[] designs, String fileSaveName, AbstractProblem problem, boolean arteryProb, double sideNum, int numHeurObjs, int numHeurConstr) throws IOException {

        File csvFile = new File(fileSaveName);
        FileWriter csvWrite = new FileWriter(csvFile);

        // CSV file headers
        csvWrite.append("Full Design");
        csvWrite.append(",");
        csvWrite.append("Norm. Obj. 1");
        csvWrite.append(",");
        csvWrite.append("Norm. Obj. 2");
        csvWrite.append(",");
        csvWrite.append("Feasibility");
        csvWrite.append(",");
        csvWrite.append("Connectivity");
        csvWrite.append(",");
        if (!arteryProb) {
            csvWrite.append("Stiffness Ratio");
            csvWrite.append(",");
        }
        csvWrite.append("\n");

        for (int j = 0; j < designs.length; j++) {
            // Evaluate current solution
            problem.evaluate(designs[j]);

            TrussRepeatableArchitecture currentArch = new TrussRepeatableArchitecture(designs[j], sideNum, numHeurObjs, numHeurConstr);
            boolean[] currentBooleanDesign = currentArch.getCompleteBooleanArrayFromSolutionNxN(designs[j]);
            csvWrite.append(convertBooleanArrayToBitstring(currentBooleanDesign));
            csvWrite.append(",");

            csvWrite.append(Double.toString(designs[j].getObjective(0)));
            csvWrite.append(",");

            csvWrite.append(Double.toString(designs[j].getObjective(1)));
            csvWrite.append(",");

            csvWrite.append(Double.toString(1.0 - (Double) designs[j].getAttribute("FeasibilityViolation")));
            csvWrite.append(",");

            csvWrite.append(Double.toString(1.0 - (Double) designs[j].getAttribute("ConnectivityViolation")));
            csvWrite.append(",");

            if (!arteryProb) {
                csvWrite.append(Double.toString((Double) designs[j].getAttribute("StiffnessRatioViolation")));
                csvWrite.append(",");
            }

            csvWrite.append("\n");

            System.out.println("Evaluated design no. " + j);
        }
        csvWrite.flush();
        csvWrite.close();

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

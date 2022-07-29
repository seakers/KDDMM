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
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.RealVariable;
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
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
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

        RunMode runMode = RunMode.EpsilonMOEA;

        boolean finalPopulationOnly = false; // used only for epsilon MOEA

        /**
         * modelChoice = 0 --> Fibre Stiffness Model
         *             = 1 --> Truss Stiffness Model
         *             = 2 --> Beam Model
         */
        int modelChoice = 1; // Fibre stiffness model cannot be used for the artery problem

        boolean arteryProblem = false; // Solve the artery optimization (otherwise the original truss problem is solved)
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
        int numDesigns = 100; // used only for random data generation

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

        String csvFilename = "";

        switch (runMode) {
            case RandomPopulation:
                for (int i = 0; i < numRuns; i++) {
                    System.out.println("Generating random data for run " + i);

                    initialization = new RandomInitialization(trussProblem, numDesigns);
                    Solution[] population = initialization.initialize();

                    csvFilename = csvPath + File.separator + "result" + File.separator + "random_operator_index_data" + fileSaveNameProblem + fileSaveNameModel + i + ".csv";

                    int[] nfes = new int[population.length];

                    computeAndStoreIntoCsv(population, nfes, csvFilename, heuristicOperators, trussProblem, arteryProblem, true, sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
                }
                engine.close();
                break;

            case EpsilonMOEA:
                for (int i = 0; i < numRuns; i++) {
                    System.out.println("Generating Epsilon MOEA data for run " + i);

                    String fileProbPath = "Truss Problem";
                    if (arteryProblem) {
                        fileProbPath = "Artery Problem";
                    }

                    String fileModelPath = "Fibre Model";
                    if (modelChoice == 1) {
                        fileModelPath = "Truss Model";
                    } else if (modelChoice == 2) {
                        fileModelPath = "Beam Model";
                    }

                    String fileNamePop = "_fullpop";
                    if (finalPopulationOnly) {
                        fileNamePop = "";
                    }

                    String dataFilePath = csvPath + File.separator + "result" + File.separator + fileProbPath + File.separator + "Constant Radii" + File.separator + fileModelPath + File.separator + "Epsilon MOEA - Metrics\\";
                    String dataFileName = "EpsilonMOEA_emoea_" + i + fileSaveNameProblem + fileSaveNameModel + fileNamePop + ".csv";

                    ArrayList<Solution> populationList = new ArrayList<>();
                    ArrayList<Integer> nfesList = new ArrayList<>();

                    // Read appropriate csv file
                    List<List<String>> rows = new ArrayList<>();
                    try (Scanner scanner = new Scanner(new File(dataFilePath + dataFileName))) {
                        while (scanner.hasNextLine()) {
                            rows.add(getRecordFromLine(scanner.nextLine()));
                        }
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }

                    boolean header = true;
                    for (List<String> row : rows) {
                        if (header) {
                            header = false;
                            continue;
                        }
                        Solution newSolution = trussProblem.newSolution();
                        TrussRepeatableArchitecture newDesign = new TrussRepeatableArchitecture(newSolution, sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
                        populateVariablesFromString(row.get(0), newDesign);

                        populationList.add(newDesign);
                        if (!finalPopulationOnly) {
                            nfesList.add(Integer.parseInt(row.get(1)));
                        }
                    }

                    Solution[] solutions = new Solution[populationList.size()];
                    int[] nfeVals = new int[nfesList.size()];
                    for (int k = 0; k < populationList.size(); k++) {
                        solutions[k] = populationList.get(k);
                        if (!finalPopulationOnly) {
                            nfeVals[k] = nfesList.get(k);
                        }
                    }

                    csvFilename = csvPath + File.separator + "result" + File.separator + "EpsilonMOEA_emoea_" + fileSaveNameProblem + fileSaveNameModel + i + fileNamePop + ".csv";

                    computeAndStoreIntoCsv(solutions, nfeVals, csvFilename, heuristicOperators, trussProblem, arteryProblem, finalPopulationOnly, sideNodeNumber, numberOfHeuristicObjectives, numberOfHeuristicConstraints);
                }
                engine.close();
                break;

            default :
                throw new IllegalStateException("Unrecognized run mode");
        }
    }

    private static void computeAndStoreIntoCsv(Solution[] designs, int[] nfes, String fileSaveName, Variation[] heuristicOps, AbstractProblem problem, boolean arteryProb, boolean finalPop, double sideNum, int numHeurObjs, int numHeurConstr) throws IOException {

        Variation repairPartialCollapsibility = heuristicOps[0];
        Variation repairNodalProperties = heuristicOps[1];
        Variation repairOrientation = heuristicOps[2];
        Variation repairIntersection = heuristicOps[3];

        File csvFile = new File(fileSaveName);
        FileWriter csvWrite = new FileWriter(csvFile);

        // CSV file headers
        if (!finalPop) {
            csvWrite.append("NFE");
            csvWrite.append(",");
        }
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
        if (!arteryProb) {
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
        if (!arteryProb) {
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
        if (!arteryProb) {
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
        if (!arteryProb) {
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
        if (!arteryProb) {
            csvWrite.append("Stiffness Ratio");
            csvWrite.append(",");
        }
        csvWrite.append("\n");

        for (int j = 0; j < designs.length; j++) {
            // Evaluate current solution
            problem.evaluate(designs[j]);

            // Operate using Partial Collapsibility Repair Operator and evaluate
            Solution partialCollapsibilitySolution = repairPartialCollapsibility.evolve(new Solution[]{designs[j]})[0];
            problem.evaluate(partialCollapsibilitySolution);

            // Operate using Nodal Properties Repair Operator and evaluate
            Solution nodalPropertiesSolution = repairNodalProperties.evolve(new Solution[]{designs[j]})[0];
            problem.evaluate(nodalPropertiesSolution);

            // Operate using Orientation Repair Operator and evaluate
            Solution orientationSolution = repairOrientation.evolve(new Solution[]{designs[j]})[0];
            problem.evaluate(orientationSolution);

            // Operate using Intersection Repair Operator and evaluate
            Solution intersectionSolution = repairIntersection.evolve(new Solution[]{designs[j]})[0];
            problem.evaluate(intersectionSolution);

            // Populate row of CSV file
            if (!finalPop) {
                csvWrite.append(Integer.toString(nfes[j]));
                csvWrite.append(",");
            }

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

            TrussRepeatableArchitecture partCollArch = new TrussRepeatableArchitecture(partialCollapsibilitySolution, sideNum, numHeurObjs, numHeurConstr);
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

            if (!arteryProb) {
                csvWrite.append(Double.toString((Double) partialCollapsibilitySolution.getAttribute("StiffnessRatioViolation")));
                csvWrite.append(",");
            }

            TrussRepeatableArchitecture nodalPropArch = new TrussRepeatableArchitecture(nodalPropertiesSolution, sideNum, numHeurObjs, numHeurConstr);
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

            if (!arteryProb) {
                csvWrite.append(Double.toString((Double) nodalPropertiesSolution.getAttribute("StiffnessRatioViolation")));
                csvWrite.append(",");
            }

            TrussRepeatableArchitecture orientArch = new TrussRepeatableArchitecture(orientationSolution, sideNum, numHeurObjs, numHeurConstr);
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

            if (!arteryProb) {
                csvWrite.append(Double.toString((Double) orientationSolution.getAttribute("StiffnessRatioViolation")));
                csvWrite.append(",");
            }

            TrussRepeatableArchitecture intersArch = new TrussRepeatableArchitecture(intersectionSolution, sideNum, numHeurObjs, numHeurConstr);
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

            if (!arteryProb) {
                csvWrite.append(",");
                csvWrite.append(Double.toString((Double) intersectionSolution.getAttribute("StiffnessRatioViolation")));
            }
            csvWrite.append("\n");

            System.out.println("Evaluated design no. " + j);
        }
        csvWrite.flush();
        csvWrite.close();

    }

    private static List<String> getRecordFromLine(String line) {
        List<String> architectures = new ArrayList<>();
        try (Scanner rowScanner = new Scanner(line)) {
            rowScanner.useDelimiter(",");
            while (rowScanner.hasNext()) {
                architectures.add(rowScanner.next());
            }
        }
        return architectures;
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

    private static void populateVariablesFromString (String desString, TrussRepeatableArchitecture solution) {
        boolean[] completeBoolean = new boolean[desString.length()];
        for (int i = 0; i < desString.length(); i++) {
            String bit = desString.substring(i,i+1);
            completeBoolean[i] = bit.equalsIgnoreCase("1");
        }
        boolean[] repeatableBoolean = solution.getRepeatableBooleanArrayFromCompleteArray(completeBoolean);

        // Populate solution with bits from desString
        for (int i = 0; i < repeatableBoolean.length; i++) {
            BinaryVariable var = new BinaryVariable(1);
            if (repeatableBoolean[i]) {
                var.set(0, true);
            } else {
                var.set(0, false);
            }
            solution.setVariable(i, var);
        }
    }

    public enum RunMode{
        RandomPopulation,
        EpsilonMOEA,
    }
}

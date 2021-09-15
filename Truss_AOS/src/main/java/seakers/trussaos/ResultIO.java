package seakers.trussaos;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.moeaframework.analysis.collector.Accumulator;
import org.moeaframework.analysis.collector.InstrumentedAlgorithm;
import org.moeaframework.core.Population;
import org.moeaframework.core.PopulationIO;
import org.moeaframework.core.Solution;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;

public class ResultIO implements Serializable {
    private static final long serialVersionUID = -2048768998854760056L;
    private static String problemClassName;
    private final double sidenum;
    private final int numHeurObjectives;
    private final int numHeurConstraints;
    //private static AbstractProblem problem;

    public ResultIO(String problemClass, double sidenum, int numHeuristicObjectives, int numHeuristicConstraints) {
        problemClassName = problemClass;
        this.sidenum = sidenum;
        this.numHeurConstraints = numHeuristicConstraints;
        this.numHeurObjectives = numHeuristicObjectives;
    }

    public static void saveSearchMetrics(InstrumentedAlgorithm instAlgorithm, String filename) {
        Accumulator accum = instAlgorithm.getAccumulator();
        File results = new File(filename + ".res");
        System.out.println("Saving metrics");
        try {
            FileWriter writer = new FileWriter(results);

            try {
                Set<String> keys = accum.keySet();
                Iterator keyIter = keys.iterator();

                while(true) {
                    if (!keyIter.hasNext()) {
                        writer.flush();
                        break;
                    }

                    String key = (String)keyIter.next();
                    int dataSize = accum.size(key);
                    writer.append(key).append(",");

                    for(int i = 0; i < dataSize; ++i) {
                        writer.append(accum.get(key, i).toString());
                        if (i + 1 < dataSize) {
                            writer.append(",");
                        }
                    }

                    writer.append("\n");
                }
            } catch (Throwable var11) {
                try {
                    writer.close();
                } catch (Throwable var10) {
                    var11.addSuppressed(var10);
                }

                throw var11;
            }

            writer.close();
        } catch (IOException var12) {
            Logger.getLogger(ResultIO.class.getName()).log(Level.SEVERE, (String)null, var12);
        }

    }

    public static void saveObjectives(Population pop, String filename) {
        System.out.println("Saving objectives");

        try {
            PopulationIO.writeObjectives(new File(filename + ".obj"), pop);
        } catch (IOException var3) {
            Logger.getLogger(ResultIO.class.getName()).log(Level.SEVERE, (String)null, var3);
        }

    }

    public static Population readObjectives(File file) throws IOException {
        return PopulationIO.readObjectives(file);
    }

    public static void savePopulation(Population pop, String filename) {
        System.out.println("Saving population");

        try {
            PopulationIO.write(new File(filename + ".pop"), pop);
        } catch (IOException var3) {
            Logger.getLogger(ResultIO.class.getName()).log(Level.SEVERE, (String)null, var3);
        }

    }

    public static Population loadPopulation(String filename) throws IOException {
        return PopulationIO.read(new File(filename));
    }

    public void saveFinalResultToCsv(Population pop, String filename) throws IOException, ExecutionException, InterruptedException {
        System.out.println("Saving final population to csv file");
        File csvFile = new File(filename);
        //ConstantRadiusTrussProblem aosProblem = new ConstantRadiusTrussProblem(filename,useFibreStiffness,targetStiffnessRatio,engine);
        FileWriter csvWrite = new FileWriter(csvFile);
        csvWrite.append("Full Design");
        csvWrite.append(",");
        if (problemClassName.equals("ConstantRadiusTrussProblem2")){
            csvWrite.append("Penalized Stiffness Objective");
        } else if (problemClassName.equals("ConstantRadiusTrussProblem")){
            csvWrite.append("Penalized Target Ratio Objective");
        } else if (problemClassName.equals("ConstantRadiusArteryProblem")) {
            csvWrite.append("Penalized Stiffness Objective");
        }
        csvWrite.append(",");
        if (problemClassName.equals("ConstantRadiusTrussProblem2")){
            csvWrite.append("Penalized Volume Fraction Objective");
        } else if (problemClassName.equals("ConstantRadiusTrussProblem")){
            csvWrite.append("Penalized Stiffness Objective");
        } else if (problemClassName.equals("ConstantRadiusArteryProblem")) {
            csvWrite.append("Penalized Deviation Objective");
        }
        csvWrite.append(",");
        if (problemClassName.equals("ConstantRadiusTrussProblem2")){
            csvWrite.append("True Stiffness Objective");
        } else if (problemClassName.equals("ConstantRadiusTrussProblem")){
            csvWrite.append("True Target Ratio Objective");
        } else if (problemClassName.equals("ConstantRadiusArteryProblem")) {
            csvWrite.append("True Stiffness Objective");
        }
        csvWrite.append(",");
        if (problemClassName.equals("ConstantRadiusTrussProblem2")){
            csvWrite.append("True Volume Fraction Objective");
        } else if (problemClassName.equals("ConstantRadiusTrussProblem")){
            csvWrite.append("True Stiffness Objective");
        } else if (problemClassName.equals("ConstantRadiusArteryProblem")) {
            csvWrite.append("True Deviation Objective");
        }
        csvWrite.append(",");
        csvWrite.append("Feasibility Score");
        csvWrite.append(",");
        csvWrite.append("Connectivity Score");
        csvWrite.append(",");
        csvWrite.append("Absolute Stiffness Ratio Difference");
        csvWrite.append(",");
        csvWrite.append("Partial Collapsibility Score");
        csvWrite.append(",");
        csvWrite.append("Nodal Properties Score");
        csvWrite.append(",");
        csvWrite.append("Orientation Score");
        csvWrite.append(",");
        csvWrite.append("Intersection Score");
        csvWrite.append("\n");
        //String[] designPopulation = new String[pop.size()];
        //double[][] objectives = new double[1][2];
        //double[] feasibilityScores = new double[pop.size()];
        //double[] stabilityScores = new double[pop.size()];

        for (int i = 0; i < pop.size(); i++) {
            Solution currentSltn = pop.get(i);

            TrussRepeatableArchitecture currentArch = new TrussRepeatableArchitecture(currentSltn,sidenum,numHeurObjectives,numHeurConstraints);
            boolean[] currentBooleanDesign = currentArch.getCompleteBooleanArrayFromSolutionNxN(currentSltn);
            //designPopulation[i] = Arrays.toString(currentBooleanDesign);
            double[] currentObjectives = currentSltn.getObjectives();

            //System.arraycopy(currentObjectives, 0, objectives[i], 0, 2);
            //int[][] currentConnectivityArray = currentArch.getConnectivityArrayFromSolution(currentSltn);
            //double currentFeasibilityScore = problem.getFeasibilityScore(currentConnectivityArray);
            //double currentStabilityScore = problem.getStabilityScore(currentConnectivityArray);

            double currentFeasibilityScore = 1.0d - (double)currentSltn.getAttribute("FeasibilityViolation");
            double currentConnectivityScore = 1.0d - (double)currentSltn.getAttribute("ConnectivityViolation");
            double currentStiffnessRatioDifference = (double)currentSltn.getAttribute("StiffnessRatioViolation");
            double currentPartialCollapsibilityScore = 1.0d - (double)currentSltn.getAttribute("PartialCollapsibilityViolation");
            double currentNodalPropertiesScore = 1.0d - (double)currentSltn.getAttribute("NodalPropertiesViolation");
            double currentOrientationScore = 1.0d - (double)currentSltn.getAttribute("OrientationViolation");
            double currentIntersectionScore = 1.0d - (double)currentSltn.getAttribute("IntersectionViolation");
            double currentTrueStiffnessObjective = (double)currentSltn.getAttribute("TrueObjective1");
            double currentTrueVolumeFractionObjective = (double)currentSltn.getAttribute("TrueObjective2");
            csvWrite.append(convertBooleanArrayToBitstring(currentBooleanDesign));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentObjectives[0]));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentObjectives[1]));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentTrueStiffnessObjective));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentTrueVolumeFractionObjective));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentFeasibilityScore));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentConnectivityScore));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentStiffnessRatioDifference));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentPartialCollapsibilityScore));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentNodalPropertiesScore));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentOrientationScore));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentIntersectionScore));
            csvWrite.append("\n");
        }
        csvWrite.flush();
        csvWrite.close();
    }

    public void savePopulationHistoryToCsv (HashSet<Solution> solutionSet, String filename) throws IOException, ExecutionException, InterruptedException {
        System.out.println("Saving population history to csv file");
        File csvFile = new File(filename);
        FileWriter csvWrite = new FileWriter(csvFile);
        csvWrite.append("Full Design");
        csvWrite.append(",");
        csvWrite.append("NFE");
        csvWrite.append(",");
        if (problemClassName.equals("ConstantRadiusTrussProblem2")){
            csvWrite.append("Penalized Stiffness Objective");
        } else if (problemClassName.equals("ConstantRadiusTrussProblem")){
            csvWrite.append("Penalized Target Ratio Objective");
        } else if (problemClassName.equals("ConstantRadiusArteryProblem")) {
            csvWrite.append("Penalized Stiffness Objective");
        }
        csvWrite.append(",");
        if (problemClassName.equals("ConstantRadiusTrussProblem2")){
            csvWrite.append("Penalized Volume Fraction Objective");
        } else if (problemClassName.equals("ConstantRadiusTrussProblem")){
            csvWrite.append("Penalized Stiffness Objective");
        } else if (problemClassName.equals("ConstantRadiusArteryProblem")) {
            csvWrite.append("Penalized Deviation Objective");
        }
        csvWrite.append(",");
        if (problemClassName.equals("ConstantRadiusTrussProblem2")){
            csvWrite.append("True Stiffness Objective");
        } else if (problemClassName.equals("ConstantRadiusTrussProblem")){
            csvWrite.append("True Target Ratio Objective");
        } else if (problemClassName.equals("ConstantRadiusArteryProblem")) {
            csvWrite.append("True Stiffness Objective");
        }
        csvWrite.append(",");
        if (problemClassName.equals("ConstantRadiusTrussProblem2")){
            csvWrite.append("True Volume Fraction Objective");
        } else if (problemClassName.equals("ConstantRadiusTrussProblem")){
            csvWrite.append("True Stiffness Objective");
        } else if (problemClassName.equals("ConstantRadiusArteryProblem")) {
            csvWrite.append("True Deviation Objective");
        }
        csvWrite.append(",");
        csvWrite.append("Feasibility Score");
        csvWrite.append(",");
        csvWrite.append("Connectivity Score");
        csvWrite.append(",");
        csvWrite.append("Absolute Stiffness Ratio Difference");
        csvWrite.append(",");
        csvWrite.append("Partial Collapsibility Score");
        csvWrite.append(",");
        csvWrite.append("Nodal Properties Score");
        csvWrite.append(",");
        csvWrite.append("Orientation Score");
        csvWrite.append(",");
        csvWrite.append("Intersection Score");
        csvWrite.append("\n");

        Iterator iter = solutionSet.iterator();

        while (iter.hasNext()){
            Solution currentSolution = (Solution) iter.next();
            TrussRepeatableArchitecture currentArch = new TrussRepeatableArchitecture(currentSolution,sidenum,numHeurObjectives,numHeurConstraints);
            boolean[] currentBooleanDesign = currentArch.getCompleteBooleanArrayFromSolutionNxN(currentSolution);
            csvWrite.append(convertBooleanArrayToBitstring(currentBooleanDesign));
            csvWrite.append(",");
            int numFunctionEvals = (int) currentSolution.getAttribute("NFE");
            csvWrite.append(Integer.toString(numFunctionEvals));
            csvWrite.append(",");
            double[] currentObjectives = currentSolution.getObjectives();

            //System.arraycopy(currentObjectives, 0, objectives[i], 0, 2);
            //int[][] currentConnectivityArray = currentArch.getConnectivityArrayFromSolution(currentSolution);
            //double currentFeasibilityScore = problem.getFeasibilityScore(currentConnectivityArray);
            //double currentStabilityScore = problem.getStabilityScore(currentConnectivityArray);

            double currentFeasibilityScore = 1.0d - (double)currentSolution.getAttribute("FeasibilityViolation");
            double currentConnectivityScore = 1.0d - (double)currentSolution.getAttribute("ConnectivityViolation");
            double currentStiffnessRatioDifference = (double)currentSolution.getAttribute("StiffnessRatioViolation");
            double currentPartialCollapsibilityScore = 1.0d - (double)currentSolution.getAttribute("PartialCollapsibilityViolation");
            double currentNodalPropertiesScore = 1.0d - (double)currentSolution.getAttribute("NodalPropertiesViolation");
            double currentOrientationScore = 1.0d - (double)currentSolution.getAttribute("OrientationViolation");
            double currentIntersectionScore = 1.0d - (double)currentSolution.getAttribute("IntersectionViolation");
            double currentTrueStiffnessObjective = (double)currentSolution.getAttribute("TrueObjective1");
            double currentTrueVolumeFractionObjective = (double)currentSolution.getAttribute("TrueObjective2");

            csvWrite.append(Double.toString(currentObjectives[0]));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentObjectives[1]));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentTrueStiffnessObjective));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentTrueVolumeFractionObjective));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentFeasibilityScore));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentConnectivityScore));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentStiffnessRatioDifference));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentPartialCollapsibilityScore));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentNodalPropertiesScore));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentOrientationScore));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentIntersectionScore));
            csvWrite.append("\n");
        }
        csvWrite.flush();
        csvWrite.close();
    }

    private String convertBooleanArrayToBitstring (boolean[] design) {
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

    //public void saveCreditHistoryCsv(CreditHistory creditHistory, String filename) throws IOException {
        //System.out.println("Saving credit history to csv file");
        //File csvFile = new File(filename);
        //FileWriter csvWrite = new FileWriter(csvFile);

    //}
}

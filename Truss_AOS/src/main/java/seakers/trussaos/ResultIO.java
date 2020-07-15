package seakers.trussaos;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;
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
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import com.mathworks.engine.*;
import seakers.aos.history.CreditHistory;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;

public class ResultIO implements Serializable {
    private static final long serialVersionUID = -2048768998854760056L;

    private static MatlabEngine engine;

    private static boolean useFibreStiffness;

    private static double targetStiffnessRatio;

    private static TrussAOSProblem problem;

    public ResultIO(TrussAOSProblem runProblem, MatlabEngine eng, boolean fibreStiffness, double targetRatio) {
        engine = eng;
        useFibreStiffness = fibreStiffness;
        targetStiffnessRatio = targetRatio;
        problem = runProblem;
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
        System.out.println("Saving results to csv file");
        File csvFile = new File(filename);
        //TrussAOSProblem aosProblem = new TrussAOSProblem(filename,useFibreStiffness,targetStiffnessRatio,engine);
        FileWriter csvWrite = new FileWriter(csvFile);
        csvWrite.append("Full Design");
        csvWrite.append(",");
        csvWrite.append("Penalized Target Ratio Objective");
        csvWrite.append(",");
        csvWrite.append("Penalized Stiffness Objective");
        csvWrite.append(",");
        csvWrite.append("Feasibility Score");
        csvWrite.append(",");
        csvWrite.append("Stability Score");
        csvWrite.append("\n");
        //String[] designPopulation = new String[pop.size()];
        double[][] objectives = new double[1][2];
        //double[] feasibilityScores = new double[pop.size()];
        //double[] stabilityScores = new double[pop.size()];

        for (int i = 0; i < pop.size(); i++) {
            Solution currentSltn = pop.get(i);
            TrussRepeatableArchitecture currentArch = new TrussRepeatableArchitecture(currentSltn);
            boolean[] currentBooleanDesign = currentArch.getBooleanDesignFromSolution(currentSltn);
            //designPopulation[i] = Arrays.toString(currentBooleanDesign);
            double[] currentObjectives = currentSltn.getObjectives();
            //System.arraycopy(currentObjectives, 0, objectives[i], 0, 2);
            int[][] currentConnectivityArray = currentArch.getConnectivityArrayFromSolution(currentSltn);
            double currentFeasibilityScore = problem.getFeasibilityScore(currentConnectivityArray,engine);
            double currentStabilityScore = problem.getStabilityScore(currentConnectivityArray,engine);
            csvWrite.append(convertBooleanArrayToBitstring(currentBooleanDesign));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentObjectives[0]));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentObjectives[1]));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentFeasibilityScore));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentStabilityScore));
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
        csvWrite.append("Penalized Target Ratio Objective");
        csvWrite.append(",");
        csvWrite.append("Penalized Stiffness Objective");
        csvWrite.append(",");
        csvWrite.append("Feasibility Score");
        csvWrite.append(",");
        csvWrite.append("Stability Score");
        csvWrite.append("\n");

        Iterator iter = solutionSet.iterator();

        while (iter.hasNext()){
            Solution currentSolution = (Solution) iter.next();
            TrussRepeatableArchitecture currentArch = new TrussRepeatableArchitecture(currentSolution);
            boolean[] currentBooleanDesign = currentArch.getBooleanDesignFromSolution(currentSolution);
            csvWrite.append(convertBooleanArrayToBitstring(currentBooleanDesign));
            csvWrite.append(",");
            int numFunctionEvals = (int) currentSolution.getAttribute("NFE");
            csvWrite.append(Integer.toString(numFunctionEvals));
            csvWrite.append(",");
            double[] currentObjectives = currentSolution.getObjectives();
            //System.arraycopy(currentObjectives, 0, objectives[i], 0, 2);
            int[][] currentConnectivityArray = currentArch.getConnectivityArrayFromSolution(currentSolution);
            double currentFeasibilityScore = problem.getFeasibilityScore(currentConnectivityArray,engine);
            double currentStabilityScore = problem.getStabilityScore(currentConnectivityArray,engine);
            csvWrite.append(Double.toString(currentObjectives[0]));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentObjectives[1]));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentFeasibilityScore));
            csvWrite.append(",");
            csvWrite.append(Double.toString(currentStabilityScore));
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

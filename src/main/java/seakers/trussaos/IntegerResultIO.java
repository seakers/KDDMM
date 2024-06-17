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

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.moeaframework.analysis.collector.Accumulator;
import org.moeaframework.analysis.collector.InstrumentedAlgorithm;
import org.moeaframework.core.Population;
import org.moeaframework.core.PopulationIO;
import org.moeaframework.core.Solution;
import com.mathworks.engine.*;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.problem.AbstractProblem;
import seakers.trussaos.architecture.IntegerRepeatableArchitecture;

public class IntegerResultIO implements Serializable{
    private static final long serialVersionUID = -2048768998946870386L;
    private static String problemClassName;
    private final double sidenum;
    private final double[] radii;
    private final double[] youngsModulii;
    private final int numHeurObjectives;
    private final int numHeurConstraints;

    public IntegerResultIO(String problemClass, double sidenum, double[] radii, double[] youngsModulii, int numHeurObjectives, int numHeurConstraints) {
        problemClassName = problemClass;
        this.sidenum = sidenum;
        this.radii = radii;
        this.youngsModulii = youngsModulii;
        this.numHeurObjectives = numHeurObjectives;
        this.numHeurConstraints = numHeurConstraints;
    }

    public void saveFinalResultToCsvInteger(Population pop, String filename) throws IOException, ExecutionException, InterruptedException {
        System.out.println("Saving final population to csv file");
        File csvFile = new File(filename);
        FileWriter csvWrite = new FileWriter(csvFile);
        int totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sidenum*sidenum))/(CombinatoricsUtils.factorial((int) ((sidenum*sidenum) - 2)) * CombinatoricsUtils.factorial(2)));;
        for (int i = 0; i < totalNumberOfMembers+1; i++) {
            String header = "Variable " + String.valueOf(i+1);
            csvWrite.append(header);
            csvWrite.append(",");
        }
        csvWrite.append("Penalized Stiffness Objective");
        csvWrite.append(",");
        csvWrite.append("Penalized Volume Fraction Objective");
        csvWrite.append(",");
        csvWrite.append("True Stiffness Objective");
        csvWrite.append(",");
        csvWrite.append("True Volume Fraction Objective");
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
        csvWrite.append("\n");
        for (int j = 0; j < pop.size(); j++) {
            Solution currentSltn = pop.get(j);
            IntegerRepeatableArchitecture architecture = new IntegerRepeatableArchitecture(currentSltn, sidenum, numHeurObjectives, numHeurConstraints, radii, youngsModulii);
            double[] completeRadiusArray = architecture.getCompleteIntegerRadiusArrayFromSolutionNxN(currentSltn);
            double[] currentObjectives = currentSltn.getObjectives();
            double currentFeasibilityScore = 1.0d - (double)currentSltn.getAttribute("FeasibilityViolation");
            double currentConnectivityScore = 1.0d - (double)currentSltn.getAttribute("ConnectivityViolation");
            double currentStiffnessRatioDifference = (double)currentSltn.getAttribute("StiffnessRatioViolation");
            double currentPartialCollapsibilityScore = 1.0d - (double)currentSltn.getAttribute("PartialCollapsibilityViolation");
            double currentNodalPropertiesScore = 1.0d - (double)currentSltn.getAttribute("NodalPropertiesViolation");
            double currentOrientationScore = 1.0d - (double)currentSltn.getAttribute("OrientationViolation");
            double currentTrueStiffnessObjective = (double)currentSltn.getAttribute("TrueObjective1");
            double currentTrueVolumeFractionObjective = (double)currentSltn.getAttribute("TrueObjective2");
            for (int k = 0; k < totalNumberOfMembers; k++) {
                csvWrite.append(Double.toString(completeRadiusArray[k]));
                csvWrite.append(",");
            }
            csvWrite.append(Double.toString(youngsModulii[EncodingUtils.getInt(currentSltn.getVariable(totalNumberOfMembers+1))]));
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
            csvWrite.append("\n");
        }
        csvWrite.flush();
        csvWrite.close();
    }

    public void savePopulationHistoryToCsvInteger (HashSet<Solution> solutionSet, String filename) throws IOException, ExecutionException, InterruptedException {
        System.out.println("Saving population history to csv file");
        int totalNumberOfMembers = (int) (CombinatoricsUtils.factorial((int) (sidenum*sidenum))/(CombinatoricsUtils.factorial((int) ((sidenum*sidenum) - 2)) * CombinatoricsUtils.factorial(2)));;
        File csvFile = new File(filename);
        FileWriter csvWrite = new FileWriter(csvFile);
        for (int i = 0; i < totalNumberOfMembers+1; i++) {
            String header = "Variable " + String.valueOf(i+1);
            csvWrite.append(header);
            csvWrite.append(",");
        }
        csvWrite.append("NFE");
        csvWrite.append(",");
        csvWrite.append("Penalized Stiffness Objective");
        csvWrite.append(",");
        csvWrite.append("Penalized Volume Fraction Objective");
        csvWrite.append(",");
        csvWrite.append("True Stiffness Objective");
        csvWrite.append(",");
        csvWrite.append("True Volume Fraction Objective");
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
        csvWrite.append("\n");

        Iterator iter = solutionSet.iterator();

        while (iter.hasNext()) {
            Solution currentSltn = (Solution) iter.next();
            IntegerRepeatableArchitecture architecture = new IntegerRepeatableArchitecture(currentSltn, sidenum, numHeurObjectives, numHeurConstraints, radii, youngsModulii);
            double[] completeRadiusArray = architecture.getCompleteIntegerRadiusArrayFromSolutionNxN(currentSltn);
            double[] currentObjectives = currentSltn.getObjectives();
            double currentFeasibilityScore = 1.0d - (double)currentSltn.getAttribute("FeasibilityViolation");
            double currentConnectivityScore = 1.0d - (double)currentSltn.getAttribute("ConnectivityViolation");
            double currentStiffnessRatioDifference = (double)currentSltn.getAttribute("StiffnessRatioViolation");
            double currentPartialCollapsibilityScore = 1.0d - (double)currentSltn.getAttribute("PartialCollapsibilityViolation");
            double currentNodalPropertiesScore = 1.0d - (double)currentSltn.getAttribute("NodalPropertiesViolation");
            double currentOrientationScore = 1.0d - (double)currentSltn.getAttribute("OrientationViolation");
            double currentTrueStiffnessObjective = (double)currentSltn.getAttribute("TrueObjective1");
            double currentTrueVolumeFractionObjective = (double)currentSltn.getAttribute("TrueObjective2");
            for (int k = 0; k < totalNumberOfMembers; k++) {
                csvWrite.append(Double.toString(completeRadiusArray[k]));
                csvWrite.append(",");
            }
            csvWrite.append(Double.toString(youngsModulii[EncodingUtils.getInt(currentSltn.getVariable(totalNumberOfMembers+1))]));
            csvWrite.append(",");
            int numFunctionEvals = (int) currentSltn.getAttribute("NFE");
            csvWrite.append(Integer.toString(numFunctionEvals));
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
            csvWrite.append("\n");
        }
        csvWrite.flush();
        csvWrite.close();
    }

}

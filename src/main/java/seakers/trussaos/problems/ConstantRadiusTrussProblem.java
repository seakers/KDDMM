package seakers.trussaos.problems;

import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.problem.AbstractProblem;
import com.mathworks.engine.*;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.lang.Math.*;

/**
 * The problem class for the Truss Optimization problem for the 2D 3x3 nodal grid case with constant and equal radius for all members.
 * The objectives are to minimize [|c22/c11 - c_target|, -c22/v_f] with feasibility constraint
 *
 * @author roshan94
 */

public class ConstantRadiusTrussProblem extends AbstractProblem {

    private final int modelSelection;

    private final String csvSavePath;

    private final MatlabEngine engine;

    private final double sel;

    private final double sidenum;

    private final double nucFac;

    private final double radius;

    private final double YoungsModulus;

    private final double targetStiffnessRatio;

    private final boolean constrainFeasibility;

    private final boolean constrainStability;

    private final boolean constrainOrientation;

    private final int numHeurObjectives;

    private final int numHeurConstraints;

    private final double[][] NodalPositionArray = new double[9][2];

    public ConstantRadiusTrussProblem(String savePath, int modelSelection, double targetCRatio, MatlabEngine eng, boolean constrainFeasibility, boolean constrainStability, boolean constrainOrientation, int numHeurObjectives, int numHeurConstraints) {

        super(32,2);

        this.modelSelection = modelSelection;
        this.constrainFeasibility = constrainFeasibility;
        this.constrainStability = constrainStability;
        this.constrainOrientation = constrainOrientation;
        this.numHeurObjectives = numHeurObjectives;
        this.numHeurConstraints = numHeurConstraints;

        this.radius = 0.00005;
        this.YoungsModulus = 10000.0;
        this.sel = 0.05;
        this.sidenum = 3.0;
        this.nucFac = 1.0;
        this.targetStiffnessRatio = targetCRatio;
        this.csvSavePath = savePath;
        this.engine = eng;

        double[][] RelativeNodalPositions = {{0,0},{0,0.5},{0,1},{0.5,0},{0.5,0.5},{0.5,1},{1,0},{1,0.5},{1,1}};
        for (int i = 0; i < RelativeNodalPositions.length; i++){
            for (int j = 0; j < RelativeNodalPositions[0].length; j++){
            this.NodalPositionArray[i][j] = RelativeNodalPositions[i][j] * sel;
            }
        }
    }

    public ConstantRadiusTrussProblem(String savePath, int modelSelection, int numVariables, int numHeurObjectives, int numHeurConstraints, double rad, double sideLength, double E, double sideNodeNum, double nucFac, double targetCRatio, MatlabEngine eng, boolean constrainFeasibility, boolean constrainStability, boolean constrainOrientation) {

        super(numVariables,2);

        this.modelSelection = modelSelection;
        this.constrainFeasibility = constrainFeasibility;
        this.constrainStability = constrainStability;
        this.constrainOrientation = constrainOrientation;
        this.numHeurObjectives = numHeurObjectives;
        this.numHeurConstraints = numHeurConstraints;

        this.radius = rad;
        this.sel = sideLength;
        this.YoungsModulus = E;
        this.sidenum = sideNodeNum;
        this.nucFac = nucFac;
        this.targetStiffnessRatio = targetCRatio;
        this.csvSavePath = savePath;
        this.engine = eng;

        double[][] RelativeNodalPositions = {{0,0},{0,0.5},{0,1},{0.5,0},{0.5,0.5},{0.5,1},{1,0},{1,0.5},{1,1}};
        for (int i = 0; i < RelativeNodalPositions.length; i++){
            for (int j = 0; j < RelativeNodalPositions[0].length; j++){
                this.NodalPositionArray[i][j] = RelativeNodalPositions[i][j] * sel;
            }
        }
    }

    @Override
    public void evaluate(Solution sltn) {
        TrussRepeatableArchitecture trussArch = new TrussRepeatableArchitecture(sltn, sidenum, numHeurObjectives, numHeurConstraints);
        double[][] designConnArray = trussArch.getConnectivityArrayFromSolution(sltn);
        //int designTrussCount = trussArch.getTrussCountFromSolution(sltn);
        //System.out.println(Arrays.deepToString(designConnArray));
        //System.out.println(designTrussCount);

        //MatlabEngine eng = null;
        double C11;
        double C22;
        // int C11Integer;
        // int C22Integer;
        double penaltyFactor = 1;
        boolean useVariableRadiiModels = true;

        //try {
            //engine.startMatlab();
        //} catch (EngineException | InterruptedException e) {
            //e.printStackTrace();
        //}
        //MatlabEngine eng = MatlabEngine.startMatlab();

        double volFrac = 0.0;
        if (modelSelection == 0) {
            Object[] outputs = null;
            if (useVariableRadiiModels) {
                double[] radiusArray = new double[designConnArray.length];
                for (int i = 0; i < designConnArray.length; i++) {
                    radiusArray[i] = radius;
                }
                try {
                    //assert eng != null;
                    outputs = engine.feval(3,"fiberStiffnessModel_rVar_V3_mod",sel,radiusArray,YoungsModulus,designConnArray,nucFac,sidenum);
                } catch (InterruptedException | ExecutionException | NullPointerException e) {
                    e.printStackTrace();
                }
                C11 = (double)outputs[0];
                C22 = (double)outputs[1];
                volFrac = (double)outputs[2];
            } else {
                try {
                    //assert eng != null;
                    outputs = engine.feval(2,"fiberStiffnessModel",sel,radius,YoungsModulus,designConnArray,sidenum,nucFac);
                } catch (InterruptedException | ExecutionException | NullPointerException e) {
                    e.printStackTrace();
                }
                C11 = (double)outputs[0];
                C22 = (double)outputs[1];
                try {
                    volFrac = getVolumeFraction(designConnArray);
                } catch (ExecutionException | InterruptedException | NullPointerException e) {
                    e.printStackTrace();
                }
            }

            // C11Integer = (int)outputs[0];
            // C22Integer = (int)outputs[1];
            // C11 = (double)C11Integer;
            // C22 = (double)C22Integer;

            penaltyFactor = 1.5;
        }
        else if (modelSelection == 1) {
            double[][] stiffnessMatrix;
            Object[] outputs = null;
            if (useVariableRadiiModels) {
                double[] radiusArray = new double[designConnArray.length];
                for (int i = 0; i < designConnArray.length; i++) {
                    radiusArray[i] = radius;
                }
                try {
                    outputs = engine.feval(2,"trussMetaCalc_NxN_rVar_AVar_mod",nucFac,sidenum,sel,radiusArray,YoungsModulus,designConnArray);
                } catch (InterruptedException | ExecutionException | NullPointerException e) {
                    e.printStackTrace();
                }
                stiffnessMatrix = (double[][])outputs[0];
                C11 = stiffnessMatrix[0][0];
                C22 = stiffnessMatrix[1][1];
                if (Double.isNaN(C11) | Double.isNaN(C22)) {
                    C11 = 1e-6;
                    C22 = 1e-3;
                }
                volFrac = (double)outputs[1];
            } else {
                double area = Math.PI * radius * radius;
                stiffnessMatrix = new double[3][3];
                try {
                    outputs = engine.feval(3,"generateC",sel,radius,NodalPositionArray,designConnArray,area,YoungsModulus,stiffnessMatrix);
                } catch (InterruptedException | ExecutionException | NullPointerException e) {
                    e.printStackTrace();
                }
                stiffnessMatrix = (double[][])outputs[0];
                C11 = stiffnessMatrix[0][0];
                C22 = stiffnessMatrix[1][1];
                if (Double.isNaN(C11) | Double.isNaN(C22)) {
                    C11 = 1e-6;
                    C22 = 1e-3;
                }
                try {
                    volFrac = getVolumeFraction(designConnArray);
                } catch (ExecutionException | InterruptedException | NullPointerException e) {
                    e.printStackTrace();
                }
            }

        } else { // ANSYS APDL Beam Model
            double[][] stiffnessMatrix;
            Object[] outputs = null;
            try {
                outputs = engine.feval("APDL_2D_NxN_V1",sel,sidenum,radius,YoungsModulus,designConnArray);
            } catch (InterruptedException | ExecutionException | NullPointerException e) {
                e.printStackTrace();
            }
            stiffnessMatrix = (double[][])outputs[0];
            C11 = stiffnessMatrix[0][0];
            C22 = stiffnessMatrix[1][1];
            if (Double.isNaN(C11) | Double.isNaN(C22) | (C22 >= YoungsModulus) | (C11 >= YoungsModulus)) {
                C11 = 1e-6;
                C22 = 1e-3;
            }
            try {
                volFrac = getVolumeFraction(designConnArray);
            } catch (ExecutionException | InterruptedException | NullPointerException e) {
                e.printStackTrace();
            }
        }

        double designFeasibilityScore = 0.0;
        try {
            designFeasibilityScore = getFeasibilityScore(designConnArray);
        } catch (ExecutionException | InterruptedException | NullPointerException e) {
            e.printStackTrace();
        }
        double designStabilityScore = 0.0;
        try {
            designStabilityScore = getStabilityScore(designConnArray);
        } catch (ExecutionException | InterruptedException | NullPointerException e) {
            e.printStackTrace();
        }

        double designOrientationScore = 0.0;
        try {
            designOrientationScore = getOrientationScore(designConnArray);
        } catch (ExecutionException | InterruptedException e) {
            e.printStackTrace();
        }

        double penaltyFeasibility = 0.0;
        if (constrainFeasibility) {
            penaltyFeasibility = Math.log10(Math.abs(designFeasibilityScore))/20;
        }
        double penaltyStability = 0.0;
        if (constrainStability) {
            penaltyStability = Math.log10(Math.abs(designStabilityScore))/20;
        }

        double penaltyOrientation = 0.0;
        if (constrainOrientation) {
            penaltyOrientation = Math.log10(Math.abs(designOrientationScore))/20;
        }

        double[] objectives = new double[2];
        double penalty = 0;
        if (constrainFeasibility && constrainStability && !constrainOrientation) {
            penalty = (penaltyFeasibility + penaltyStability)/2;
        } else if (constrainFeasibility && constrainOrientation && !constrainStability) {
            penalty = (penaltyFeasibility + penaltyOrientation)/2;
        } else if (constrainOrientation && constrainStability && !constrainFeasibility) {
            penalty = (penaltyStability + penaltyOrientation)/2;
        } else if (constrainFeasibility && !constrainStability && !constrainOrientation) {
            penalty = penaltyFeasibility;
        } else if (constrainStability && !constrainFeasibility && !constrainOrientation) {
            penalty = penaltyStability;
        } else if (constrainOrientation && !constrainFeasibility && !constrainStability) {
            penalty = penaltyOrientation;
        } else if (constrainFeasibility && constrainStability && constrainOrientation) {
            penalty = (penaltyStability + penaltyOrientation + penaltyFeasibility)/3;
        }

        double[] trueObjectives = new double[2];
        trueObjectives[0] = Math.abs((C11/C22) - targetStiffnessRatio);
        trueObjectives[1] = C22/volFrac;

        objectives[0] = trueObjectives[0]/25 - penaltyFactor*penalty;
        objectives[1] = -trueObjectives[1]/(200*YoungsModulus) - penaltyFactor*penalty;

        sltn.setObjectives(objectives);

        //HashMap<String, Object> constrainthandling = new HashMap<String, Object>();
        //constrainthandling.put("FeasibilityViolation", 1 - designFeasibilityScore);
        //constrainthandling.put("StabilityViolation", 1 - designStabilityScore);
        //sltn.addAttributes(constrainthandling);

        sltn.setAttribute("FeasibilityViolation", 1.0 - designFeasibilityScore);
        sltn.setAttribute("StabilityViolation", 1.0 - designStabilityScore);
        sltn.setAttribute("OrientationViolation", 1.0 - designOrientationScore);
        sltn.setAttribute("TrueObjective1", trueObjectives[0]);
        sltn.setAttribute("TrueObjective2", trueObjectives[1]);
    }

    public double getFeasibilityScore (double[][] designConnectivityArray) throws ExecutionException, InterruptedException, NullPointerException {
        return (double)engine.feval("feasibility_checker_nonbinary_V2",NodalPositionArray,designConnectivityArray);
    }

    public double getStabilityScore (double[][] designConnectivityArray) throws ExecutionException, InterruptedException, NullPointerException {
        // Object stabilityOutput;
        // stabilityOutput = eng.feval("stabilityTester_2D_V6_1",sidenum,designConnectivityArray,NodalPositionArray,sel);
        // return (double)stabilityOutput;
        return (double)engine.feval("stabilityTester_2D_V7",sidenum,designConnectivityArray,NodalPositionArray,sel);
        // return (double)eng.feval("stabilityTester_2D_updated",sidenum,designConnectivityArray,NodalPositionArray);
    }

    public double getOrientationScore (double[][] designConnectivityArray) throws ExecutionException, InterruptedException {
        //Object orientationOutput;
        //orientationOutput = engine.feval("orientationHeuristic_V2",NodalPositionArray,designConnectivityArray,targetStiffnessRatio);
        return (double) engine.feval("orientationHeuristic",NodalPositionArray,designConnectivityArray,targetStiffnessRatio);
    }

    private double getVolumeFraction (double[][] designConnectivityArray) throws ExecutionException, InterruptedException, NullPointerException {
        // This function does not consider edge-correction. Edge-corrected volume fraction obtained from models.
        return (double)engine.feval("calcVF",NodalPositionArray,designConnectivityArray,radius,sel);
    }

    public double[][] getNodalConnectivityArray () {
        return NodalPositionArray;
    }

    @Override
    public Solution newSolution() {
        synchronized (PRNG.getRandom()) {
            Solution newSol = new Solution(this.numberOfVariables, 2);
            for (int i = 0; i < this.numberOfVariables; i++) {
                BinaryVariable newVar = new BinaryVariable(1);
                newVar.set(0, PRNG.nextBoolean());
                newSol.setVariable(i, newVar);
            }
            return newSol;
        }
    }
}

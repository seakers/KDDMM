package seakers.trussaos.problems;

import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.problem.AbstractProblem;
import org.apache.commons.math3.util.CombinatoricsUtils;
import com.mathworks.engine.*;
import seakers.trussaos.architecture.IntegerRepeatableArchitecture;
import seakers.architecture.util.IntegerVariable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.ExecutionException;

import static java.lang.Double.NaN;
import static java.lang.Double.min;

/**
 * The alternate problem class for the variable radius truss optimization problem (where each design variable
 * is an integer indicating choice of radius value for the corresponding truss member) using the 2D NxN lattice unit cell.
 * The objectives are to minimize [-c22, v_f] with feasbility and |c22/c11 - c_target|, feasibility and connectivity
 * constraint handling
 * NOTE: the variable radii truss model assumes a 3x3 node grid (sidenum = 3) and can take different values for cell repetitions
 *
 * @author roshan94
 */

public class IntegerTrussProblem extends AbstractProblem {

    private final int modelSelection;
    private final String csvSavePath;
    private final MatlabEngine engine;
    private final double sel;
    private final double sidenum;
    private final double nucFac;
    private final double[] radii;
    private final double[] YoungsModulii;
    private final double targetStiffnessRatio;
    private final boolean[][] heuristicsConstrainedBooleans;
    private final int numHeurObjectives;
    private final int numHeurConstraints;
    private double[][] NodalPositionArray;

    public IntegerTrussProblem(String savePath, int modelSelection, int numTrussVariables, int numHeurObjectives, int numHeurConstraints, double targetCRatio, MatlabEngine eng, double[] radii, double sideLength, double[] YoungsModulii, double sideNodeNum, double nucFac, boolean[][] constrainHeuristics) {
        // numVariables is decided based on the side node number and is computed in the main class
        super(numTrussVariables+1,2+numHeurObjectives,3+numHeurConstraints);
        this.radii = radii;
        this.modelSelection = modelSelection;
        this.numHeurObjectives = numHeurObjectives;
        this.numHeurConstraints = numHeurConstraints;
        this.heuristicsConstrainedBooleans = constrainHeuristics;

        this.sel = sideLength;
        this.YoungsModulii = YoungsModulii;
        this.sidenum = sideNodeNum;
        this.nucFac = nucFac;
        this.targetStiffnessRatio = targetCRatio;
        this.csvSavePath = savePath;
        this.engine = eng;

        this.NodalPositionArray = getNodalConnectivityArrayFromSidenum((int) sidenum, sel);
    }

    @Override
    public void evaluate(Solution sltn) {
        IntegerRepeatableArchitecture arch = new IntegerRepeatableArchitecture(sltn, sidenum, numHeurObjectives, numHeurConstraints, radii, YoungsModulii);
        double YoungsModulusChoice = YoungsModulii[arch.getYoungsModulusChoice()];
        double[][] designConnArray = arch.getFullConnectivityArrayIntegerRadii();
        double[] designFullRadiusArray = arch.getFullRadiusArray();
        double[] radiusArray = getNonZeroRadiusArray(designFullRadiusArray);

        double C11;
        double C22;
        double penaltyFactor = 1;
        double volFrac = 0.0;
        double heuristicBiasFactor = 1.0; // For certain heuristics (unless specified within function)

        boolean stiffnessIsBad = false;
        boolean noStiffnessIndicated = false;

        if (modelSelection == 0) { // Fibre Stiffness Model

            Object[] outputs = null;
            try {
                //assert eng != null;
                outputs = engine.feval(3,"fiberStiffnessModel_rVar_V3",sel,radiusArray,YoungsModulusChoice,designConnArray,nucFac,sidenum);
            } catch (InterruptedException | ExecutionException | NullPointerException e) {
                e.printStackTrace();
            }
            // C11Integer = (int)outputs[0];
            // C22Integer = (int)outputs[1];
            // C11 = (double)C11Integer;
            // C22 = (double)C22Integer;
            C11 = (double)outputs[0];
            C22 = (double)outputs[1];
            volFrac = (double)outputs[2];

            penaltyFactor = 1.5;
        } else if (modelSelection == 1) { // Truss Stiffness Model
            double[][] stiffnessMatrix;
            Object[] outputs = null;
            // For the model, sidenum = (2*nucFac) + 1
            try {
                //outputs = engine.feval(2,"trussMetaCalc_NxN_rVar_AVar",nucFac,sel,radiusArray,YoungsModulus,designConnArray);
                outputs = engine.feval(2,"trussMetaCalc_NxN_1UC_rVar_AVar",sidenum,sel,radiusArray,YoungsModulusChoice,designConnArray);
            } catch (InterruptedException | ExecutionException | NullPointerException e) {
                e.printStackTrace();
            }
            stiffnessMatrix = (double[][])outputs[0];
            C11 = stiffnessMatrix[0][0];
            C22 = stiffnessMatrix[1][1];
            if (Double.isNaN(C11) | Double.isNaN(C22) | (C22 >= YoungsModulusChoice) | (C11 >= YoungsModulusChoice)) {
                stiffnessIsBad = true;
                C11 = 1e-6;
                C22 = 1e-3;
            }
            if ((C22 < 1) | (C11 < 1)) {
                noStiffnessIndicated = true;
            }
            try {
                volFrac = getVolumeFractionNxN(designConnArray,radiusArray);
            } catch (ExecutionException | InterruptedException | NullPointerException e) {
                e.printStackTrace();
            }
        } else { // ANSYS APDL Beam Model - Not created for the integer radii case therefore don't use!!!
            double[][] stiffnessMatrix;
            Object[] outputs = null;
            try {
                outputs = engine.feval("Beam_2D_NxN_PBC",sel,sidenum,radiusArray[0],YoungsModulusChoice,designConnArray);
            } catch (InterruptedException | ExecutionException | NullPointerException e) {
                e.printStackTrace();
            }
            stiffnessMatrix = (double[][])outputs;
            C11 = stiffnessMatrix[0][0];
            C22 = stiffnessMatrix[1][1];
            if (Double.isNaN(C11) | Double.isNaN(C22) | (C22 >= YoungsModulusChoice) | (C11 >= YoungsModulusChoice)) {
                stiffnessIsBad = true;
                C11 = 1e-6;
                C22 = 1e-3;
            }
            try {
                volFrac = getVolumeFractionNxN(designConnArray,radiusArray);
            } catch (ExecutionException | InterruptedException | NullPointerException e) {
                e.printStackTrace();
            }
        }
        double designFeasibilityScore = 0.0;
        try {
            designFeasibilityScore = getFeasibilityScoreVariableRadii(designConnArray);
        } catch (ExecutionException | InterruptedException | NullPointerException e) {
            e.printStackTrace();
        }

        double designConnectivityScore = 0.0;
        try {
            designConnectivityScore = getConnectivityScoreVariableRadii(designConnArray,heuristicBiasFactor);
        } catch (ExecutionException | InterruptedException | NullPointerException e) {
            e.printStackTrace();
        }

        double designPartialCollapsibilityScore = 0.0;
        try {
            designPartialCollapsibilityScore = getPartialCollapsibilityScoreVariableRadii(designConnArray);
        } catch (ExecutionException | InterruptedException e) {
            e.printStackTrace();
        }

        double designNodalPropertiesScore = 0.0;
        try {
            designNodalPropertiesScore = getNodalPropertiesScoreVariableRadii(designConnArray,heuristicBiasFactor);
        } catch (ExecutionException | InterruptedException e) {
            e.printStackTrace();
        }

        double designOrientationScore = 0.0;
        try {
            designOrientationScore = getOrientationScoreVariableRadii(designConnArray);
        } catch (ExecutionException | InterruptedException e) {
            e.printStackTrace();
        }

        double designIntersectionScore = 0.0;
        try {
            designIntersectionScore = getIntersectionScore(designConnArray);
        } catch (ExecutionException | InterruptedException e) {
            e.printStackTrace();
        }

        double[] heuristicObjectives = {-designPartialCollapsibilityScore, -designNodalPropertiesScore, -designOrientationScore, -designIntersectionScore};

        double penaltyFeasibility = Math.log10(Math.abs(designFeasibilityScore))/16;

        double penaltyConnectivity = Math.log10(Math.abs(designConnectivityScore))/16;

        double penaltyPartialCollapsibility = 0.0;
        if (heuristicsConstrainedBooleans[0][0] || heuristicsConstrainedBooleans[0][4] || heuristicsConstrainedBooleans[0][5]) {
            penaltyPartialCollapsibility = Math.log10(Math.abs(designPartialCollapsibilityScore))/16;
        }

        double penaltyNodalProperties = 0.0;
        if (heuristicsConstrainedBooleans[1][0] || heuristicsConstrainedBooleans[1][4] || heuristicsConstrainedBooleans[1][5]) {
            penaltyNodalProperties = Math.log10(Math.abs(designNodalPropertiesScore))/16;
        }

        double penaltyOrientation = 0.0;
        if (heuristicsConstrainedBooleans[2][0] || heuristicsConstrainedBooleans[2][4] || heuristicsConstrainedBooleans[2][5]) {
            penaltyOrientation = Math.log10(Math.abs(designOrientationScore))/16;
        }

        double penaltyIntersection = 0.0;
        if (heuristicsConstrainedBooleans[3][0] || heuristicsConstrainedBooleans[3][4] || heuristicsConstrainedBooleans[3][5]) {
            penaltyIntersection = Math.log10(Math.abs(designIntersectionScore))/16;
        }

        double[] heuristicPenalties = {penaltyPartialCollapsibility, penaltyNodalProperties, penaltyOrientation, penaltyIntersection};
        double[] heuristicConstraints = {1 - designPartialCollapsibilityScore, 1 - designNodalPropertiesScore, 1 - designOrientationScore, 1 - designIntersectionScore};

        double[] objectives = new double[2+numHeurObjectives];
        double[] constraints = new double[3+numHeurConstraints];

        //double constraintWeight = 10;
        //double stiffnessRatioConstraintWeight = 10;
        double heuristicWeight = 1;

        double absoluteStiffnessRatioDifference = Math.abs((C22/C11) - targetStiffnessRatio);
        //double penalty = -constraintWeight*penaltyConnectivity;
        double penalty = 0;
        double heuristicPenalty = 0;
        double numHeuristicsInteriorPenalty = 0;

        for (int i = 0; i < heuristicsConstrainedBooleans.length; i++) {
            if (heuristicsConstrainedBooleans[i][0]) {
                heuristicPenalty -= heuristicPenalties[i];
                numHeuristicsInteriorPenalty++;
            }
        }
        if (numHeuristicsInteriorPenalty != 0) {
            heuristicPenalty /= numHeuristicsInteriorPenalty;
        }
        penalty += heuristicWeight*heuristicPenalty;

        double[] trueObjectives = new double[2+numHeurObjectives];
        trueObjectives[0] = C22;
        trueObjectives[1] = volFrac;

        objectives[0] = -trueObjectives[0]/YoungsModulusChoice + penaltyFactor*penalty;
        objectives[1] = trueObjectives[1] + penaltyFactor*penalty;

        int heurObjectiveCount = 0;
        if (numHeurObjectives > 0) {
            for (int i = 0; i < 4; i++) {
                if (heuristicsConstrainedBooleans[i][4]) {
                    trueObjectives[2+heurObjectiveCount] = heuristicObjectives[i];
                    heurObjectiveCount++;
                }
            }
            for (int j = 0; j < numHeurObjectives; j++) {
                objectives[2+j] = trueObjectives[2+j] + penaltyFactor*penalty;
            }
        }

        constraints[0] = 1 - designFeasibilityScore;
        constraints[1] = 1 - designConnectivityScore;
        constraints[2] = absoluteStiffnessRatioDifference;
        int heurConstraintCount = 0;
        if (numHeurConstraints > 0) {
            for (int i = 0; i < 4; i++) {
                if (heuristicsConstrainedBooleans[i][5]) {
                    constraints[3+heurConstraintCount] = heuristicConstraints[i];
                    heurConstraintCount++;
                }
            }
        }

        sltn.setObjectives(objectives);
        sltn.setConstraints(constraints);

        sltn.setAttribute("FeasibilityViolation", 1.0 - designFeasibilityScore);
        sltn.setAttribute("ConnectivityViolation", 1.0 - designConnectivityScore);
        sltn.setAttribute("StiffnessRatioViolation", absoluteStiffnessRatioDifference);
        sltn.setAttribute("PartialCollapsibilityViolation", 1.0 - designPartialCollapsibilityScore);
        sltn.setAttribute("NodalPropertiesViolation", 1.0 - designNodalPropertiesScore);
        sltn.setAttribute("OrientationViolation", 1.0 - designOrientationScore);
        sltn.setAttribute("IntersectionViolation", 1.0 - designIntersectionScore);
        if (stiffnessIsBad | noStiffnessIndicated) {
            sltn.setAttribute("TrueObjective1", NaN);
        } else {
            sltn.setAttribute("TrueObjective1", trueObjectives[0]);

        }
        sltn.setAttribute("TrueObjective2", trueObjectives[1]);
    }

    public double getFeasibilityScoreVariableRadii(double[][] designConnectivityArray) throws ExecutionException, InterruptedException, NullPointerException {
        return (double)engine.feval("feasibility_checker_nonbinary_V5",NodalPositionArray,designConnectivityArray,sel,sidenum);
    }

    public double getConnectivityScoreVariableRadii(double[][] designConnectivityArray, double biasFactor) throws ExecutionException, InterruptedException, NullPointerException {
        return (double)engine.feval("connectivityConstraint_PBC_2D",sidenum,NodalPositionArray,designConnectivityArray,sel,biasFactor);
    }

    private double getPartialCollapsibilityScoreVariableRadii(double[][] designConnectivityArray) throws ExecutionException, InterruptedException, NullPointerException {
        double collapsibilityBiasFactor= 0.5;
        return (double)engine.feval("partCollapseHeuristic_2D",sidenum,designConnectivityArray,NodalPositionArray,sel,collapsibilityBiasFactor);
    }

    private double getNodalPropertiesScoreVariableRadii(double[][] designConnectivityArray, double biasFactor) throws ExecutionException, InterruptedException, NullPointerException {
        return (double)engine.feval("connectivityHeuristic_2D",sidenum,NodalPositionArray,designConnectivityArray,sel, biasFactor);
    }

    private double getOrientationScoreVariableRadii(double[][] designConnectivityArray) throws ExecutionException, InterruptedException {
        //Object orientationOutput;
        //orientationOutput = engine.feval("orientationHeuristic_V2",NodalPositionArray,designConnectivityArray,targetStiffnessRatio);
        return (double) engine.feval("orientationHeuristicNorm",NodalPositionArray,designConnectivityArray,sel,targetStiffnessRatio);
    }

    private double getIntersectionScore(double[][] designConnectivityArray) throws ExecutionException, InterruptedException {
        return (double) engine.feval("intersectHeuristic",NodalPositionArray,designConnectivityArray);
    }

    private double getVolumeFractionNxN (double[][] designConnectivityArray, double[] radiusArray) throws ExecutionException, InterruptedException {
        double volFrac = engine.feval("calcVF_NxN_feasOnly",designConnectivityArray,radiusArray,sel,sidenum);
        return min(volFrac,1.0D); // Capping the volume fraction value to 1
    }

    public double[][] getNodalConnectivityArray () {
        return NodalPositionArray;
    }

    private double[][] getNodalConnectivityArrayFromSidenum (int sideNodeNumber, double sideElementLength) {
        double[][] NodalPositionArray = new double[sideNodeNumber*sideNodeNumber][2];

        for (int i = 0; i < NodalPositionArray.length; i++){
            NodalPositionArray[i][0] = (Math.floor(i/sideNodeNumber))/(sideNodeNumber-1) * sideElementLength;
        }

        for (int j = 0; j < NodalPositionArray.length; j++){
            switch (j%sideNodeNumber) {
                case 0:
                    NodalPositionArray[j][1] = 0;
                    break;
                default:
                    double remainder = j%sideNodeNumber;
                    NodalPositionArray[j][1] = (remainder/(sideNodeNumber-1)) * sideElementLength;
                    break;
            }
        }

        return NodalPositionArray;
    }

    private double[] getNonZeroRadiusArray(double[] fullRadiusArray) {
        ArrayList<Double> nonZeroRadiusArrayList = new ArrayList<Double>();
        for (int i = 0; i < fullRadiusArray.length; i++) {
            if (fullRadiusArray[i] != 0D) {
                nonZeroRadiusArrayList.add(fullRadiusArray[i]);
            }
        }
        return nonZeroRadiusArrayList.stream().mapToDouble(Double::doubleValue).toArray();
    }

    @Override
    public Solution newSolution() {
        Solution newSol = new Solution(this.numberOfVariables, 2+numHeurObjectives,3+numHeurConstraints);
        //Random rnd = new Random();
        for (int i = 0; i <= this.numberOfVariables-2; i++) {
            RealVariable newVar = new RealVariable(0,radii.length);
            newVar.setValue(PRNG.nextInt(radii.length));
            //IntegerVariable newVar = new IntegerVariable(PRNG.nextInt(radii.length),0,radii.length);
            newSol.setVariable(i, newVar);
        }
        RealVariable newEVar = new RealVariable(0,YoungsModulii.length);
        newEVar.setValue(PRNG.nextInt(YoungsModulii.length));
        //IntegerVariable newEVar = new IntegerVariable(PRNG.nextInt(YoungsModulii.length),0,YoungsModulii.length);
        newSol.setVariable(this.numberOfVariables-1, newEVar);
        return newSol;
    }
}

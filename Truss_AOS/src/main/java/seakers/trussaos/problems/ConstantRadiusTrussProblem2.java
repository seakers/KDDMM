package seakers.trussaos.problems;

import org.apache.commons.math3.util.CombinatoricsUtils;
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

import static java.lang.Double.NaN;
import static java.lang.Double.min;

/**
 * The alternate problem class for the constant radius truss optimization problem using the 2D NxN lattice unit cell.
 * The objectives are to minimize [-c22, v_f] with feasbility, |c22/c11 - c_target| and connectivity
 * constraint handling
 * NOTE: the variable radii truss model assumes a 3x3 node grid (sidenum = 3) and can take different values for cell repetitions
 *
 * @author roshansuresh
 */

public class ConstantRadiusTrussProblem2 extends AbstractProblem {

    private final int modelSelection;
    private final String csvSavePath;
    private final MatlabEngine engine;
    private final double sel;
    private final double sidenum;
    private final double nucFac;
    private final double radius;
    private final double YoungsModulus;
    private final double targetStiffnessRatio;
    private final boolean[][] heuristicsConstrainedBooleans;
    private final int numHeurObjectives;
    private final int numHeurConstraints;
    private double[][] NodalPositionArray;

    /**
     * partialCollapsibilityConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
     * nodalPropertiesConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
     * orientationConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
     * intersectionConstrained = [interior_penalty, AOS, biased_init, ACH, objective, constraint]
     *
     * constrainHeuristics = [partialCollapsibilityConstrained, nodalPropertiesConstrained, orientationConstrained, intersectionConstrained]
     */

    public ConstantRadiusTrussProblem2(String savePath, int modelSelection, int numVariables, int numHeurObjectives, int numHeurConstraints, double targetCRatio, MatlabEngine eng, boolean[][] constrainHeuristics) {
        // For default NxN node grid
        super(numVariables,2+numHeurObjectives,3+numHeurConstraints);

        this.modelSelection = modelSelection;
        this.numHeurObjectives = numHeurObjectives;
        this.numHeurConstraints = numHeurConstraints;
        this.heuristicsConstrainedBooleans = constrainHeuristics;

        this.radius = 0.00005;
        this.YoungsModulus = 10000.0;
        this.sel = 0.05;
        this.sidenum = 3.0;
        this.nucFac = 1.0;
        this.targetStiffnessRatio = targetCRatio;
        this.csvSavePath = savePath;
        this.engine = eng;

        //double[][] RelativeNodalPositions = {{0,0},{0,0.5},{0,1},{0.5,0},{0.5,0.5},{0.5,1},{1,0},{1,0.5},{1,1}};
        //for (int i = 0; i < RelativeNodalPositions.length; i++){
            //for (int j = 0; j < RelativeNodalPositions[0].length; j++){
                //this.NodalPositionArray[i][j] = RelativeNodalPositions[i][j] * sel;
            //}
        //}
        this.NodalPositionArray = getNodalConnectivityArrayFromSidenum((int) sidenum, sel);
    }

    public ConstantRadiusTrussProblem2(String savePath, int modelSelection, int numVariables, int numHeurObjectives, int numHeurConstraints, double rad, double sideLength, double E, double sideNodeNum, double nucFac, double targetCRatio, MatlabEngine eng, boolean[][] constrainHeuristics) {

        super(numVariables,2+numHeurObjectives,3+numHeurConstraints);

        this.modelSelection = modelSelection;
        this.numHeurObjectives = numHeurObjectives;
        this.numHeurConstraints = numHeurConstraints;
        this.heuristicsConstrainedBooleans = constrainHeuristics;

        this.radius = rad;
        this.sel = sideLength;
        this.YoungsModulus = E;
        this.sidenum = sideNodeNum;
        this.nucFac = nucFac;
        this.targetStiffnessRatio = targetCRatio;
        this.csvSavePath = savePath;
        this.engine = eng;

        //double[][] RelativeNodalPositions = {{0,0},{0,0.5},{0,1},{0.5,0},{0.5,0.5},{0.5,1},{1,0},{1,0.5},{1,1}};
        //for (int i = 0; i < RelativeNodalPositions.length; i++){
            //for (int j = 0; j < RelativeNodalPositions[0].length; j++){
                //this.NodalPositionArray[i][j] = RelativeNodalPositions[i][j] * sel;
            //}
        //}
        this.NodalPositionArray = getNodalConnectivityArrayFromSidenum((int) sidenum, sel);
    }

    @Override
    public void evaluate(Solution sltn) {
        TrussRepeatableArchitecture trussArch = new TrussRepeatableArchitecture(sltn, sidenum, numHeurObjectives, numHeurConstraints);
        double[][] designConnArray = trussArch.getConnectivityArrayFromSolution(sltn);

        double C11;
        double C22;
        double penaltyFactor = 1;
        double volFrac = 0.0;
        double heuristicBiasFactor = 1.0; // For certain heuristics (unless specified within function)

        boolean useVariableRadiiModels = true;
        boolean stiffnessIsBad = false;
        boolean noStiffnessIndicated = false;

        double[] radiusArray = new double[designConnArray.length];
        for (int i = 0; i < designConnArray.length; i++) {
            radiusArray[i] = radius;
        }

        if (modelSelection == 0) { // Fibre Stiffness Model

            Object[] outputs = null;
            if (useVariableRadiiModels) {
                try {
                    //assert eng != null;
                    outputs = engine.feval(3,"fiberStiffnessModel_rVar_V3",sel,radiusArray,YoungsModulus,designConnArray,nucFac,sidenum);
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
            penaltyFactor = 1.5;
        }
        else if (modelSelection == 1) { // Truss Stiffness Model
            double[][] stiffnessMatrix;
            Object[] outputs = null;
            if (useVariableRadiiModels) {

                // For the model, sidenum = (2*nucFac) + 1
                try {
                    //outputs = engine.feval(2,"trussMetaCalc_NxN_rVar_AVar",nucFac,sel,radiusArray,YoungsModulus,designConnArray);
                    outputs = engine.feval(2,"trussMetaCalc_NxN_1UC_rVar_AVar",sidenum,sel,radiusArray,YoungsModulus,designConnArray);
                } catch (InterruptedException | ExecutionException | NullPointerException e) {
                    e.printStackTrace();
                }
                stiffnessMatrix = (double[][])outputs[0];
                C11 = stiffnessMatrix[0][0];
                C22 = stiffnessMatrix[1][1];
                if (Double.isNaN(C11) | Double.isNaN(C22) | (C22 >= YoungsModulus) | (C11 >= YoungsModulus)) {
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
                if (Double.isNaN(C11) | Double.isNaN(C22) | (C22 >= YoungsModulus) | (C11 >= YoungsModulus)) {
                    stiffnessIsBad = true;
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
                outputs = engine.feval("Beam_2D_NxN_PBC",sel,sidenum,radius,YoungsModulus,designConnArray);
            } catch (InterruptedException | ExecutionException | NullPointerException e) {
                e.printStackTrace();
            }
            stiffnessMatrix = (double[][])outputs;
            C11 = stiffnessMatrix[0][0];
            C22 = stiffnessMatrix[1][1];
            if (Double.isNaN(C11) | Double.isNaN(C22) | (C22 >= YoungsModulus) | (C11 >= YoungsModulus)) {
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

        //double penaltyFeasibility = Math.log10(Math.abs(designFeasibilityScore))/16;

        //double penaltyConnectivity = Math.log10(Math.abs(designConnectivityScore))/16;

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

        //if (heuristicsConstrainedBooleans[0][0] && heuristicsConstrainedBooleans[1][0] && !heuristicsConstrainedBooleans[2][0]) {
            //heuristicPenalty = (penaltyPartialCollapsibility + penaltyNodalProperties)/2;
        //} else if (heuristicsConstrainedBooleans[0][0] && heuristicsConstrainedBooleans[2][0] && !heuristicsConstrainedBooleans[1][0]) {
            //heuristicPenalty = (penaltyPartialCollapsibility + penaltyOrientation)/2;
        //} else if (heuristicsConstrainedBooleans[2][0] && heuristicsConstrainedBooleans[1][0] && !heuristicsConstrainedBooleans[0][0]) {
            //heuristicPenalty = (penaltyNodalProperties + penaltyOrientation)/2;
        //} else if (heuristicsConstrainedBooleans[0][0] && !heuristicsConstrainedBooleans[1][0] && !heuristicsConstrainedBooleans[2][0]) {
            //heuristicPenalty = penaltyPartialCollapsibility;
        //} else if (heuristicsConstrainedBooleans[1][0] && !heuristicsConstrainedBooleans[0][0] && !heuristicsConstrainedBooleans[2][0]) {
            //heuristicPenalty = penaltyNodalProperties;
        //} else if (heuristicsConstrainedBooleans[2][0] && !heuristicsConstrainedBooleans[0][0] && !heuristicsConstrainedBooleans[1][0]) {
            //heuristicPenalty = penaltyOrientation;
        //} else if (heuristicsConstrainedBooleans[0][0] && heuristicsConstrainedBooleans[1][0] && heuristicsConstrainedBooleans[2][0]) {
            //heuristicPenalty = (penaltyNodalProperties + penaltyOrientation + penaltyPartialCollapsibility)/3;
        //}

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

        objectives[0] = -trueObjectives[0]/YoungsModulus + penaltyFactor*penalty;
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
        //orientationOutput = engine.feval("orientationHeuristicNorm",NodalPositionArray,designConnectivityArray,sel,targetStiffnessRatio);
        return (double) engine.feval("orientationHeuristic_V2",NodalPositionArray,designConnectivityArray,targetStiffnessRatio);
    }

    private double getIntersectionScore(double[][] designConnectivityArray) throws ExecutionException, InterruptedException {
        return (double) engine.feval("intersectHeuristic",NodalPositionArray,designConnectivityArray);
    }

    private double getVolumeFraction (double[][] designConnectivityArray) throws ExecutionException, InterruptedException, NullPointerException {
        // This function does not consider edge-correction. Edge-corrected volume fraction obtained from models.
        return (double)engine.feval("calcVF",NodalPositionArray,designConnectivityArray,radius,sel);
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

    @Override
    public Solution newSolution() {
        Solution newSol = new Solution(this.numberOfVariables, 2+numHeurObjectives, 3+numHeurConstraints);
        for (int i = 0; i < this.numberOfVariables; i++) {
            BinaryVariable newVar = new BinaryVariable(1);
            newVar.set(0, PRNG.nextBoolean());
            newSol.setVariable(i, newVar);
        }
        return newSol;
    }
}

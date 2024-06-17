package seakers.trussaos.problems;

import com.mathworks.engine.MatlabEngine;
import com.mathworks.engine.MatlabExecutionException;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.problem.AbstractProblem;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;

import java.io.File;
import java.util.concurrent.ExecutionException;

import static java.lang.Double.NaN;
import static java.lang.Double.min;

/**
 * The problem class for the rabbit carotid artery problem.
 * The objectives are {C11/density, (abs(C22/C11 - 0.421) + abs(C12/C11 - 0.0745) + abs(C21/C11 - 0.0745) + abs(C61) + abs(C16) + abs(C62) + abs(C26) + abs(C66/C11 - 5.038))/8}
 * the constraints are {feasibility, connectivity}
 *
 * @author roshansuresh
 */

public class ConstantRadiusArteryProblem extends AbstractProblem {

    private final int modelSelection;
    private final String csvSavePath;
    private final MatlabEngine engine;
    private final double sel;
    private final double sidenum;
    private final double nucFac;
    private final double targetStiffnessRatio;
    private final double radius;
    private final double YoungsModulus;
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

    public ConstantRadiusArteryProblem (String savePath, int modelSelection, int numVariables, int numHeurObjectives, int numHeurConstraints, double targetStiffnessRatio, MatlabEngine eng, boolean[][] constrainHeuristics) {
        // For default NxN node grid
        super(numVariables,2+numHeurObjectives,2+numHeurConstraints);

        this.modelSelection = modelSelection;
        this.numHeurObjectives = numHeurObjectives;
        this.numHeurConstraints = numHeurConstraints;
        this.heuristicsConstrainedBooleans = constrainHeuristics;

        this.radius = 0.00005;
        this.YoungsModulus = 10000.0;
        this.sel = 0.05;
        this.sidenum = 3.0;
        this.nucFac = 1.0;
        this.targetStiffnessRatio = targetStiffnessRatio;
        this.csvSavePath = savePath;
        this.engine = eng;

        this.NodalPositionArray = getNodalConnectivityArrayFromSidenum((int) sidenum, sel);
    }

    public ConstantRadiusArteryProblem (String savePath, int modelSelection, int numVariables, int numHeurObjectives, int numHeurConstraints, double rad, double sideLength, double E, double sideNodeNum, double nucFac, double targetStiffnessRatio, MatlabEngine eng, boolean[][] constrainHeuristics) {

        super(numVariables,2+numHeurObjectives,2+numHeurConstraints);

        this.modelSelection = modelSelection;
        this.numHeurObjectives = numHeurObjectives;
        this.numHeurConstraints = numHeurConstraints;
        this.heuristicsConstrainedBooleans = constrainHeuristics;

        this.radius = rad;
        this.sel = sideLength;
        this.YoungsModulus = E;
        this.sidenum = sideNodeNum;
        this.nucFac = nucFac;
        this.targetStiffnessRatio = targetStiffnessRatio;
        this.csvSavePath = savePath;
        this.engine = eng;

        this.NodalPositionArray = getNodalConnectivityArrayFromSidenum((int) sidenum, sel);
    }

    @Override
    public void evaluate(Solution sltn) {
        TrussRepeatableArchitecture trussArch = new TrussRepeatableArchitecture(sltn, sidenum, numHeurObjectives, numHeurConstraints);
        double[][] designConnArray = trussArch.getConnectivityArrayFromSolution(sltn);
        double[] radiusArray = new double[designConnArray.length];
        for (int i = 0; i < designConnArray.length; i++) {
            radiusArray[i] = radius;
        }

        double C11 = 0;
        double C22 = 0;
        double C12 = 0;
        double C21 = 0;
        double C16 = 0;
        double C26 = 0;
        double C61 = 0;
        double C62 = 0;
        double C66 = 0;
        double penaltyFactor = 1;
        double volFrac = 0.0;
        double heuristicBiasFactor = 1.0; // For certain heuristics (unless specified within function)

        boolean useVariableRadiiModels = true;
        boolean stiffnessIsBad = false;
        boolean noStiffnessIndicated = false;

        // Current version of Fibre Stiffness model does not compute shear elements of stiffness matrix

        if (modelSelection == 1) { // Truss Stiffness Model
            double[][] stiffnessMatrix;
            Object[] outputs = null;
            if (useVariableRadiiModels) {
                try {
                    //outputs = engine.feval(2,"trussMetaCalc_NxN_rVar_AVar",nucFac,sel,radiusArray,YoungsModulus,designConnArray);
                    outputs = engine.feval(2,"trussMetaCalc_NxN_1UC_rVar_AVar",sidenum,sel,radiusArray,YoungsModulus,designConnArray);
                } catch (InterruptedException | ExecutionException | NullPointerException e) {
                    e.printStackTrace();
                }
                stiffnessMatrix = (double[][])outputs[0];
                C11 = stiffnessMatrix[0][0];
                C22 = stiffnessMatrix[1][1];
                C12 = stiffnessMatrix[0][1];
                C21 = stiffnessMatrix[1][0];
                C16 = stiffnessMatrix[0][2];
                C26 = stiffnessMatrix[1][2];
                C61 = stiffnessMatrix[2][0];
                C62 = stiffnessMatrix[2][1];
                C66 = stiffnessMatrix[2][2];
                if (Double.isNaN(C11) | Double.isNaN(C22) | (C22 >= YoungsModulus) | (C11 >= YoungsModulus)) {
                    stiffnessIsBad = true;
                    C11 = 1e-6;
                    C22 = 1e-3;
                }
                if (Double.isNaN(C11) | Double.isNaN(C12) | (C11 >= YoungsModulus) | (C12 >= YoungsModulus)) {
                    stiffnessIsBad = true;
                    C11 = 1e-6;
                    C12 = 1e-3;
                }
                if (Double.isNaN(C11) | Double.isNaN(C21) | (C11 >= YoungsModulus) | (C21 >= YoungsModulus)) {
                    stiffnessIsBad = true;
                    C11 = 1e-6;
                    C21 = 1e-3;
                }
                if (Double.isNaN(C11) | Double.isNaN(C66) | (C11 >= YoungsModulus) | (C66 >= YoungsModulus)) {
                    stiffnessIsBad = true;
                    C11 = 1e-6;
                    C66 = 1e-3;
                }
                if (Double.isNaN(C61) | Double.isNaN(C62) | Double.isNaN(C16) | Double.isNaN(C26)) {
                    stiffnessIsBad = true;
                    C61 = 1e3;
                    C62 = 1e3;
                    C16 = 1e3;
                    C26 = 1e3;
                }
                if ((C61 >= YoungsModulus) | (C62 >= YoungsModulus) | (C16 >= YoungsModulus) | (C26 >= YoungsModulus)) {
                    stiffnessIsBad = true;
                    C61 = 1e3;
                    C62 = 1e3;
                    C16 = 1e3;
                    C26 = 1e3;
                }
                if ((C22 < 1) | (C11 < 1) | (C12 < 1) | (C21 < 1)) {
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
            C12 = stiffnessMatrix[0][1];
            C21 = stiffnessMatrix[1][0];
            C16 = stiffnessMatrix[0][2];
            C26 = stiffnessMatrix[1][2];
            C61 = stiffnessMatrix[2][0];
            C62 = stiffnessMatrix[2][1];
            C66 = stiffnessMatrix[2][2];
            if (Double.isNaN(C11) | Double.isNaN(C22) | (C22 >= YoungsModulus) | (C11 >= YoungsModulus)) {
                stiffnessIsBad = true;
                C11 = 1e-6;
                C22 = 1e-3;
            }
            if (Double.isNaN(C11) | Double.isNaN(C12) | (C11 >= YoungsModulus) | (C12 >= YoungsModulus)) {
                stiffnessIsBad = true;
                C11 = 1e-6;
                C12 = 1e-3;
            }
            if (Double.isNaN(C11) | Double.isNaN(C21) | (C11 >= YoungsModulus) | (C21 >= YoungsModulus)) {
                stiffnessIsBad = true;
                C11 = 1e-6;
                C21 = 1e-3;
            }
            if (Double.isNaN(C11) | Double.isNaN(C66) | (C11 >= YoungsModulus) | (C66 >= YoungsModulus)) {
                stiffnessIsBad = true;
                C11 = 1e-6;
                C66 = 1e-3;
            }
            if (Double.isNaN(C61) | Double.isNaN(C62) | Double.isNaN(C16) | Double.isNaN(C26)) {
                stiffnessIsBad = true;
                C61 = 1e3;
                C62 = 1e3;
                C16 = 1e3;
                C26 = 1e3;
            }
            if ((C61 >= YoungsModulus) | (C62 >= YoungsModulus) | (C16 >= YoungsModulus) | (C26 >= YoungsModulus)) {
                stiffnessIsBad = true;
                C61 = 1e3;
                C62 = 1e3;
                C16 = 1e3;
                C26 = 1e3;
            }
            if ((C22 < 1) | (C11 < 1) | (C12 < 1) | (C21 < 1)) {
                noStiffnessIndicated = true;
            }
            try {
                volFrac = getVolumeFractionNxN(designConnArray, radiusArray);
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
        double[] constraints = new double[2+numHeurConstraints];

        //double constraintWeight = 10;
        //double stiffnessRatioConstraintWeight = 10;
        double heuristicWeight = 1;

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
        trueObjectives[0] = C11/volFrac;
        trueObjectives[1] = (Math.abs(C22/C11 - targetStiffnessRatio) + Math.abs(C12/C11 - 0.0745) + Math.abs(C21/C11 - 0.0745) + Math.abs(C61/YoungsModulus) + Math.abs(C16/YoungsModulus) + Math.abs(C62/YoungsModulus) + Math.abs(C26/YoungsModulus) + Math.abs(C66/C11 - 5.038))/8;

        objectives[0] = -(trueObjectives[0] - 2e5)/1e6 + penaltyFactor*penalty;
        objectives[1] = ((Math.abs(C22/C11 - targetStiffnessRatio))/6 + Math.abs(C12/C11 - 0.0745) + Math.abs(C21/C11 - 0.0745) + Math.abs(C61/1.5e5) + Math.abs(C16/9e4) + Math.abs(C62/1.5e5) + Math.abs(C26/9.5e4) + (Math.abs(C66/C11 - 5.038) - 4.5)/0.5)/8 + penaltyFactor*penalty;

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

        int heurConstraintCount = 0;
        if (numHeurConstraints > 0) {
            for (int i = 0; i < 4; i++) {
                if (heuristicsConstrainedBooleans[i][5]) {
                    constraints[2+heurConstraintCount] = heuristicConstraints[i];
                    heurConstraintCount++;
                }
            }
        }

        sltn.setObjectives(objectives);
        try {
            sltn.setConstraints(constraints);
        } catch (Exception e) {
            e.printStackTrace();
        }


        sltn.setAttribute("FeasibilityViolation", 1.0 - designFeasibilityScore);
        sltn.setAttribute("ConnectivityViolation", 1.0 - designConnectivityScore);
        sltn.setAttribute("PartialCollapsibilityViolation", 1.0 - designPartialCollapsibilityScore);
        sltn.setAttribute("NodalPropertiesViolation", 1.0 - designNodalPropertiesScore);
        sltn.setAttribute("OrientationViolation", 1.0 - designOrientationScore);
        sltn.setAttribute("IntersectionViolation", 1.0 - designIntersectionScore);
        if (stiffnessIsBad | noStiffnessIndicated) {
            sltn.setAttribute("TrueObjective1", NaN);
        } else {
            sltn.setAttribute("TrueObjective1", -(trueObjectives[0] - 2e5)/1e6);
        }
        sltn.setAttribute("TrueObjective2", trueObjectives[1]);
        sltn.setAttribute("alreadyEvaluated", true);
    }

    public double getFeasibilityScoreVariableRadii(double[][] designConnectivityArray) throws ExecutionException, InterruptedException, NullPointerException {
        double feasibilityScore = 0.0;
        try {
            feasibilityScore = (double)engine.feval("feasibility_checker_nonbinary_V5",NodalPositionArray,designConnectivityArray,sel,sidenum);
        } catch (ExecutionException | InterruptedException | NullPointerException e) {
            e.printStackTrace();
        }
        return feasibilityScore;
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
        Solution newSol = new Solution(this.numberOfVariables, 2+numHeurObjectives, 2+numHeurConstraints);
        for (int i = 0; i < this.numberOfVariables; i++) {
            BinaryVariable newVar = new BinaryVariable(1);
            newVar.set(0, PRNG.nextBoolean());
            newSol.setVariable(i, newVar);
        }
        return newSol;
    }
}

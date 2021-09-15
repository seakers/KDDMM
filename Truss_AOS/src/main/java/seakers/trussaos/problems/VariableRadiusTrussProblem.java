package seakers.trussaos.problems;

import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.problem.AbstractProblem;
import org.apache.commons.math3.util.CombinatoricsUtils;
import com.mathworks.engine.*;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;
import seakers.trussaos.architecture.VariableRadiiRepeatableArchitecture;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.lang.Math.*;

import static java.lang.Double.NaN;

/**
 * The alternate problem class for the variable radius truss optimization problem using the 2D NxN lattice unit cell.
 * The objectives are to minimize [-c22, v_f] with feasbility and |c22/c11 - c_target|, feasibility and connectivity
 * constraint handling
 * NOTE: the variable radii truss model assumes a 3x3 node grid (sidenum = 3) and can take different values for cell repetitions
 *
 * @author roshan94
 */

public class VariableRadiusTrussProblem extends AbstractProblem {

    private final boolean UseFibreStiffnessModel;
    private final String csvSavePath;
    private final MatlabEngine engine;
    private final double[] radiusLowerBounds;
    private final double[] radiusUpperBounds;
    private static double lowestRadiusFactor;
    private final double sel;
    private final double nucFac;
    private double sidenum;
    private final double YoungsModulus;
    private final double targetStiffnessRatio;
    private final boolean constrainPartialCollapsibility;
    private final boolean constrainNodalProperties;
    private final boolean constrainOrientation;
    private final double[][] NodalPositionArray;

    public VariableRadiusTrussProblem(int numVariables, double[] radiusLowerBounds, double[] radiusUpperBounds, double sidenum, String savePath, boolean FibreStiffness, double targetCRatio, MatlabEngine eng, boolean constrainPartialCollapsibility, boolean constrainNodalProperties, boolean constrainOrientation){
        // numVariables is decided based on the side node number and is computed in the main class
        super(numVariables,2);
        this.radiusLowerBounds = radiusLowerBounds;
        this.radiusUpperBounds = radiusUpperBounds;
        this.csvSavePath = savePath;
        this.UseFibreStiffnessModel = FibreStiffness;
        this.engine = eng;
        this.constrainPartialCollapsibility = constrainPartialCollapsibility;
        this.constrainNodalProperties = constrainNodalProperties;
        this.constrainOrientation = constrainOrientation;
        this.targetStiffnessRatio = targetCRatio;

        this.YoungsModulus = 1e5;
        this.nucFac = 1;
        this.sel = 0.05;
        lowestRadiusFactor = 5e-2;

        this.sidenum = sidenum;
        //int numberOfVariables = (int) (CombinatoricsUtils.factorial(sidenum*sidenum)/(CombinatoricsUtils.factorial((sidenum*sidenum) - 2) * CombinatoricsUtils.factorial(2)) - (sidenum - 1)*(sidenum  - 1));
        this.NodalPositionArray = getNodalConnectivityArray();
    }

    public VariableRadiusTrussProblem(double[] radiusLowerBounds, double[] radiusUpperBounds, String savePath, boolean FibreStiffness, int numVariables, double sideLength, double E, double nucFac, double sidenum, double targetCRatio, double smallestRadiusFactor, MatlabEngine eng, boolean constrainPartialCollapsibility, boolean constrainNodalProperties, boolean constrainOrientation) {
        // numVariables is decided based on the side node number and is computed in the main class
        super(numVariables,2);
        this.csvSavePath = savePath;
        this.radiusLowerBounds = radiusLowerBounds;
        this.radiusUpperBounds = radiusUpperBounds;
        this.UseFibreStiffnessModel = FibreStiffness;
        this.engine = eng;
        this.constrainPartialCollapsibility = constrainPartialCollapsibility;
        this.constrainNodalProperties = constrainNodalProperties;
        this.constrainOrientation = constrainOrientation;
        this.targetStiffnessRatio = targetCRatio;
        this.YoungsModulus = E;
        this.nucFac = nucFac;
        this.sel = sideLength;
        lowestRadiusFactor = smallestRadiusFactor;

        this.sidenum = sidenum;
        this.NodalPositionArray = getNodalConnectivityArray();
    }

    private double[][] getNodalConnectivityArray() {
        double[][] NodalPositionArray = new double[(int) (sidenum*sidenum)][2];

        for (int i = 0; i < NodalPositionArray.length; i++){
            NodalPositionArray[i][0] = (Math.floor(i/sidenum))/(sidenum-1) * sel;
        }

        for (int j = 0; j < NodalPositionArray.length; j++){
            switch ((int) (j%sidenum)) {
                case 0:
                    NodalPositionArray[j][1] = 0;
                    break;
                default:
                    double remainder = j%sidenum;
                    NodalPositionArray[j][1] = (remainder/(sidenum-1)) * sel;
                    break;
            }
        }

        return NodalPositionArray;
    }

    @Override
    public void evaluate(Solution sltn) {
        VariableRadiiRepeatableArchitecture architecture = new VariableRadiiRepeatableArchitecture(sltn, sidenum, radiusLowerBounds, radiusUpperBounds, lowestRadiusFactor);
        int[][] designConnArray = architecture.getFullCAFromSolutionVariableRadii();
        double[] designFullRadiusArray = architecture.getFullRadiusArray();
        double[] radiusArray = getNonZeroRadiusArray(designFullRadiusArray);

        double C11;
        double C22;
        double penaltyFactor = 1;
        double volFrac = 0.0;
        double heuristicBiasFactor = 1.0; // For certain heuristics (unless specified within function)

        boolean stiffnessIsBad = false;
        boolean noStiffnessIndicated = false;

        if (UseFibreStiffnessModel) {

            Object[] outputs = null;
            try {
                //assert eng != null;
                outputs = engine.feval(3, "fiberStiffnessModel_rVar_V3_mod", sel, radiusArray, YoungsModulus, designConnArray, nucFac, sidenum);
            } catch (InterruptedException | ExecutionException | NullPointerException e) {
                e.printStackTrace();
            }
            // C11Integer = (int)outputs[0];
            // C22Integer = (int)outputs[1];
            // C11 = (double)C11Integer;
            // C22 = (double)C22Integer;
            C11 = (double) outputs[0];
            C22 = (double) outputs[1];
            volFrac = (double) outputs[2];

            penaltyFactor = 1.5;
        } else {
            double[][] stiffnessMatrix;
            Object[] outputs = null;
            if (radiusArray.length == 0) {
                C11 = NaN;
                C22 = NaN;
            } else {
                try {
                    outputs = engine.feval(2,"trussMetaCalc_NxN_rVar_AVar",nucFac,sel,radiusArray,YoungsModulus,designConnArray);
                } catch (InterruptedException | ExecutionException | NullPointerException e) {
                    e.printStackTrace();
                }
                stiffnessMatrix = (double[][])outputs[0];
                C11 = stiffnessMatrix[0][0];
                C22 = stiffnessMatrix[1][1];
                volFrac = (double)outputs[1];
            }
            if (Double.isNaN(C11) | Double.isNaN(C22) | (C22 >= YoungsModulus) | (C11 >= YoungsModulus)) {
                stiffnessIsBad = true;
                C11 = 1e-6;
                C22 = 1e-3;
            }
            if ((C22 < 1) | (C11 < 1)) {
                noStiffnessIndicated = true;
            }
        }

        double designFeasibilityScore = 1e-16;
        double designConnectivityScore = 1e-16;
        double designPartialCollapsibilityScore = 1e-16;
        double designNodalPropertiesScore = 1e-16;
        double designOrientationScore = 1e-16;
        if (radiusArray.length != 0 ) {
            try {
                designFeasibilityScore = getFeasibilityScoreVariableRadii(designConnArray);
            } catch (ExecutionException | InterruptedException | NullPointerException e) {
                e.printStackTrace();
            }

            try {
                designConnectivityScore = getConnectivityScoreVariableRadii(designConnArray,heuristicBiasFactor);
            } catch (ExecutionException | InterruptedException | NullPointerException e) {
                e.printStackTrace();
            }

            try {
                designPartialCollapsibilityScore = getPartialCollapsibilityScoreVariableRadii(designConnArray);
            } catch (ExecutionException | InterruptedException e) {
                e.printStackTrace();
            }

            try {
                designNodalPropertiesScore = getNodalPropertiesScoreVariableRadii(designConnArray,heuristicBiasFactor);
            } catch (ExecutionException | InterruptedException e) {
                e.printStackTrace();
            }

            try {
                designOrientationScore = getOrientationScoreVariableRadii(designConnArray);
            } catch (ExecutionException | InterruptedException e) {
                e.printStackTrace();
            }

        }

        double penaltyFeasibility = Math.log10(Math.abs(designFeasibilityScore))/16;

        double penaltyConnectivity = Math.log10(Math.abs(designConnectivityScore))/16;

        double penaltyPartialCollapsibility = 0.0;
        if (constrainPartialCollapsibility) {
            penaltyPartialCollapsibility = Math.log10(Math.abs(designNodalPropertiesScore))/16;
        }

        double penaltyNodalProperties = 0.0;
        if (constrainNodalProperties) {
            penaltyNodalProperties = Math.log10(Math.abs(designNodalPropertiesScore))/16;
        }

        double penaltyOrientation = 0.0;
        if (constrainOrientation) {
            penaltyOrientation = Math.log10(Math.abs(designOrientationScore))/16;
        }

        double[] objectives = new double[2];

        double constraintWeight = 10;
        double stiffnessRatioConstraintWeight = 10;
        double heuristicWeight = 1;

        double normalizedAbsoluteStiffnessRatioDifference = Math.abs((C22/C11) - targetStiffnessRatio)/10;
        double penalty = stiffnessRatioConstraintWeight*normalizedAbsoluteStiffnessRatioDifference - constraintWeight*(penaltyFeasibility+penaltyConnectivity)/2;
        double heuristicPenalty = 0;
        if (constrainPartialCollapsibility && constrainNodalProperties && !constrainOrientation) {
            heuristicPenalty = (penaltyPartialCollapsibility + penaltyNodalProperties)/2;
        } else if (constrainPartialCollapsibility && constrainOrientation && !constrainNodalProperties) {
            heuristicPenalty = (penaltyPartialCollapsibility + penaltyOrientation)/2;
        } else if (constrainOrientation && constrainNodalProperties && !constrainPartialCollapsibility) {
            heuristicPenalty = (penaltyNodalProperties + penaltyOrientation)/2;
        } else if (constrainPartialCollapsibility && !constrainNodalProperties && !constrainOrientation) {
            heuristicPenalty = penaltyPartialCollapsibility;
        } else if (constrainNodalProperties && !constrainPartialCollapsibility && !constrainOrientation) {
            heuristicPenalty = penaltyNodalProperties;
        } else if (constrainOrientation && !constrainPartialCollapsibility && !constrainNodalProperties) {
            heuristicPenalty = penaltyOrientation;
        } else if (constrainPartialCollapsibility && constrainNodalProperties && constrainOrientation) {
            heuristicPenalty = (penaltyNodalProperties + penaltyOrientation + penaltyPartialCollapsibility)/3;
        }
        penalty += heuristicWeight*heuristicPenalty;

        double[] trueObjectives = new double[2];
        trueObjectives[0] = C22;
        trueObjectives[1] = volFrac;

        objectives[0] = -trueObjectives[0]/YoungsModulus + penaltyFactor*penalty;
        objectives[1] = trueObjectives[1]/0.96 + penaltyFactor*penalty;

        sltn.setObjectives(objectives);

        sltn.setAttribute("FeasibilityViolation", 1.0 - designFeasibilityScore);
        sltn.setAttribute("ConnectivityViolation", 1.0 - designConnectivityScore);
        sltn.setAttribute("StiffnessRatioViolation", normalizedAbsoluteStiffnessRatioDifference*10);
        sltn.setAttribute("PartialCollapsibilityViolation", 1.0 - designPartialCollapsibilityScore);
        sltn.setAttribute("NodalPropertiesViolation", 1.0 - designNodalPropertiesScore);
        sltn.setAttribute("OrientationViolation", 1.0 - designOrientationScore);
        if (stiffnessIsBad | noStiffnessIndicated) {
            sltn.setAttribute("TrueObjective1", NaN);
        } else {
            sltn.setAttribute("TrueObjective1", trueObjectives[0]);
        }
        sltn.setAttribute("TrueObjective2", trueObjectives[1]);
    }

    public double getFeasibilityScoreVariableRadii (int[][] designConnectivityArray) throws ExecutionException, InterruptedException, NullPointerException {
        return (double)engine.feval("feasibility_checker_nonbinary_V2",NodalPositionArray,designConnectivityArray);
    }

    public double getConnectivityScoreVariableRadii (int[][] designConnectivityArray, double biasFactor) throws ExecutionException, InterruptedException, NullPointerException {
        return (double)engine.feval("connectivityConstraint_NPBC_2D_V2",NodalPositionArray,designConnectivityArray,biasFactor);
    }

    public double getPartialCollapsibilityScoreVariableRadii (int[][] designConnectivityArray) throws ExecutionException, InterruptedException, NullPointerException {
        double collapsibilityBiasFactor= 0.5;
        return (double)engine.feval("partCollapseHeuristic_2D",sidenum,designConnectivityArray,NodalPositionArray,sel,collapsibilityBiasFactor);
    }

    public double getNodalPropertiesScoreVariableRadii (int[][] designConnectivityArray, double biasFactor) throws ExecutionException, InterruptedException, NullPointerException {
        return (double)engine.feval("connectivityHeuristic_2D",sidenum,NodalPositionArray,designConnectivityArray,sel,biasFactor);
    }

    public double getOrientationScoreVariableRadii (int[][] designConnectivityArray) throws ExecutionException, InterruptedException {
        //Object orientationOutput;
        //orientationOutput = engine.feval("orientationHeuristic_V2",NodalPositionArray,designConnectivityArray,targetStiffnessRatio);
        return (double) engine.feval("orientationHeuristicNorm",NodalPositionArray,designConnectivityArray,sel,targetStiffnessRatio);
    }

    public double[][] getNodalConnectivityArrayVariableRadii () {
        return NodalPositionArray;
    }

    private double[] getNonZeroRadiusArray(double[] fullRadiusArray) {
        ArrayList<Double> nonZeroRadiusArrayList = new ArrayList<Double>();
        for (int i = 0; i < fullRadiusArray.length; i++) {
            if (fullRadiusArray[i] > lowestRadiusFactor*radiusUpperBounds[i]) {
                nonZeroRadiusArrayList.add(fullRadiusArray[i]);
            }
        }
        return nonZeroRadiusArrayList.stream().mapToDouble(Double::doubleValue).toArray();
    }

    @Override
    public Solution newSolution() {
        synchronized (PRNG.getRandom()) {
            Solution newSol = new Solution(this.numberOfVariables, 2);
            //Random rnd = new Random();
            for (int i = 0; i < this.numberOfVariables; i++) {
                RealVariable newVar = new RealVariable(radiusLowerBounds[i],radiusUpperBounds[i]);
                newVar.setValue(PRNG.nextDouble(radiusLowerBounds[i],radiusUpperBounds[i]));
                newSol.setVariable(i, newVar);
            }
            return newSol;
        }

    }
}

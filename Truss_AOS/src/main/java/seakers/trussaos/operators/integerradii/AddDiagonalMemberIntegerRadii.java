package seakers.trussaos.operators.integerradii;

import com.mathworks.engine.MatlabEngine;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import seakers.trussaos.architecture.IntegerRepeatableArchitecture;

import java.util.*;
import java.util.concurrent.ExecutionException;

public class AddDiagonalMemberIntegerRadii implements Variation {

    private final boolean keepFeasible;
    private static MatlabEngine engine;
    private final double[][] nodalConnectivityArray;
    private final double sidenum;
    private final double sel;
    private final double[] radii;
    private final double[] youngsModulii;
    private final int numHeurObjectives;
    private final int numHeurConstraints;

    public AddDiagonalMemberIntegerRadii(boolean keepFeasible, double[] radii, double[] youngsModulii, MatlabEngine eng, double[][] nodalConnArray, double sidenum, double sel, int numHeurObjectives, int numHeurConstraints) {
        this.keepFeasible = keepFeasible;
        engine = eng;
        this.nodalConnectivityArray = nodalConnArray;
        this.sidenum = sidenum;
        this.sel = sel;
        this.radii = radii;
        this.youngsModulii = youngsModulii;
        this.numHeurObjectives = numHeurObjectives;
        this.numHeurConstraints = numHeurConstraints;
    }

    @Override
    public int getArity() {
        return 1;
    }

    @Override
    public Solution[] evolve(Solution[] solutions) {
        IntegerRepeatableArchitecture architecture = new IntegerRepeatableArchitecture(solutions[0], (int) sidenum, numHeurObjectives, numHeurConstraints, radii, youngsModulii);
        double[][] connectivityArray = architecture.getFullCAFromSolutionIntegerRadii();
        double[][] completeConnectivityArray = architecture.getCompleteConnectivityArrayFromSidenum();
        double[] radiusArray = architecture.getFullRadiusArray();
        int youngsModulusChoice = architecture.getYoungsModulusChoice(); // Use the same material as the parent architecture

        IntegerRepeatableArchitecture newArchitecture = null;
        AddMemberIntegerRadii addMemberObject = new AddMemberIntegerRadii(keepFeasible,radii,youngsModulii,engine,nodalConnectivityArray,sidenum,sel,numHeurObjectives,numHeurConstraints);
        double[][] updatedConnectivityArray;
        double[] updatedRadiusArray;

        if (keepFeasible) {
            boolean feasibilitySameOrBetter = false;
            int numberOfAttempts = 0;
            //ArrayList<ArrayList<Integer>> oldDesignTrussIntersections = new ArrayList<ArrayList<Integer>>();
            //ArrayList<ArrayList<Integer>> newDesignTrussIntersections = new ArrayList<ArrayList<Integer>>();
            int oldDesignTrussIntersections = 0;
            int newDesignTrussIntersections = 0;
            while (!feasibilitySameOrBetter && numberOfAttempts < 5) {
                try {
                    oldDesignTrussIntersections = addMemberObject.findNumberOfIntersectingTrusses(connectivityArray);
                } catch (ExecutionException | InterruptedException e) {
                    e.printStackTrace();
                }
                // Algorithm to find diagonal members not present and add them at random
                double[][] diagonalMembersAbsent = findAbsentDiagonalMembers(connectivityArray);
                if (diagonalMembersAbsent.length == 0) {
                    updatedConnectivityArray = connectivityArray.clone();
                    updatedRadiusArray = radiusArray.clone();
                }
                else if (diagonalMembersAbsent.length == 1) {
                    updatedConnectivityArray = addMemberObject.addMemberToConnectivityArray(connectivityArray,diagonalMembersAbsent[0]);
                    updatedRadiusArray = addRandomRadiusMemberToArray(completeConnectivityArray,diagonalMembersAbsent[0],radiusArray);
                }
                else {
                    int diagonalMemberChoice = PRNG.nextInt(diagonalMembersAbsent.length);
                    updatedConnectivityArray = addMemberObject.addMemberToConnectivityArray(connectivityArray,diagonalMembersAbsent[diagonalMemberChoice]);
                    updatedRadiusArray = addRandomRadiusMemberToArray(completeConnectivityArray,diagonalMembersAbsent[diagonalMemberChoice],radiusArray);
                }
                try {
                    newDesignTrussIntersections = addMemberObject.findNumberOfIntersectingTrusses(updatedConnectivityArray);
                } catch (ExecutionException | InterruptedException e) {
                    e.printStackTrace();
                }
                newArchitecture = architecture.getIntegerArchitectureFromRadiusArray(updatedRadiusArray, youngsModulusChoice);
                if (oldDesignTrussIntersections >= newDesignTrussIntersections) {
                    feasibilitySameOrBetter = true;
                }
                else {
                    numberOfAttempts += 1;
                }
            }
        }
        else {
            double[][] diagonalMembersAbsent = findAbsentDiagonalMembers(connectivityArray);
            if (diagonalMembersAbsent.length == 0) {
                //updatedConnectivityArray = connectivityArray.clone();
                updatedRadiusArray = radiusArray.clone();
            }
            else if (diagonalMembersAbsent.length == 1) {
                //updatedConnectivityArray = addMemberObject.addMemberToConnectivityArray(connectivityArray,diagonalMembersAbsent[0]);
                updatedRadiusArray = addRandomRadiusMemberToArray(completeConnectivityArray,diagonalMembersAbsent[0],radiusArray);
            }
            else {
                int diagonalMemberChoice = PRNG.nextInt(diagonalMembersAbsent.length);
                //updatedConnectivityArray = addMemberObject.addMemberToConnectivityArray(connectivityArray,diagonalMembersAbsent[diagonalMemberChoice]);
                updatedRadiusArray = addRandomRadiusMemberToArray(completeConnectivityArray,diagonalMembersAbsent[diagonalMemberChoice],radiusArray);
            }
            newArchitecture = architecture.getIntegerArchitectureFromRadiusArray(updatedRadiusArray,youngsModulusChoice);
        }
        return new Solution[]{newArchitecture};
    }

    private double[] getCornerNodes () {
        // Array order: {bottom left, top left, bottom right, top right}
        return new double[]{1, sidenum, (sidenum * sidenum - (sidenum - 1)), (sidenum * sidenum)};
    }

    private double[] getEdgeCentreNodes () {
        // Array order: {left, down, top, right} (assuming only odd values for sidenum)
        return new double[]{(sidenum - ((sidenum - 1)/2)), ((((sidenum - 1)/2) * sidenum) + 1), (sidenum * (sidenum + 1)/2), ((sidenum * sidenum) - ((sidenum - 1)/2))};

    }

    private double[][] getAllDiagonalMembers () {
        ArrayList<ArrayList<Double>> diagonalMembers = new ArrayList<ArrayList<Double>>();

        for (int i = 0; i < ((sidenum * sidenum) - 1); i++) {
            double nodeNumber = i+1;
            double closestTopNode = findClosestTopNode(nodeNumber);
            for (int j = (int) (nodeNumber+1); j < ((sidenum * sidenum)+1); j++) {
                if ((j <= closestTopNode) || ((j - nodeNumber)%sidenum == 0)) {
                    continue;
                } else {
                    ArrayList<Double> currentDiagonalMember = new ArrayList<Double>();
                    currentDiagonalMember.add(nodeNumber);
                    currentDiagonalMember.add((double) j);
                    diagonalMembers.add(currentDiagonalMember);
                }
            }
        }
        return convertArrayListTo2DArray(diagonalMembers);
    }

    private double findClosestTopNode (double node) {
        boolean topReached = false;
        double currentNode = node;
        while (!topReached) {
            if (currentNode%sidenum == 0) {
                topReached = true;
            } else {
                currentNode += 1;
            }
        }
        return currentNode;
    }

    private double[][] convertArrayListTo2DArray (ArrayList<ArrayList<Double>> arrayList2D) {
        double[][] array2D = new double[arrayList2D.size()][2];
        for (int i = 0; i < arrayList2D.size(); i++) {
            ArrayList<Double> row = arrayList2D.get(i);
            double[] rowArray = row.stream().mapToDouble(j->j).toArray();
            array2D[i] = rowArray;
        }
        return array2D;
    }

    private double[][] findAbsentDiagonalMembers (double[][] designConnectivityArray) {
        double[][] allDiagonalMembers = getAllDiagonalMembers();
        ArrayList<ArrayList<Double>> absentDiagonalMembers = new ArrayList<ArrayList<Double>>();
        for (int i = 0; i < allDiagonalMembers.length; i++) {
            boolean presence = checkIfMemberIsPresent(designConnectivityArray, allDiagonalMembers[i]);
            if (!presence) {
                ArrayList<Double> absentDiagonalMember = new ArrayList<>();
                absentDiagonalMember.add(allDiagonalMembers[i][0]);
                absentDiagonalMember.add(allDiagonalMembers[i][1]);
                absentDiagonalMembers.add(absentDiagonalMember);
            }
        }
        return convertArrayListTo2DArray(absentDiagonalMembers);
    }

    private boolean checkIfMemberIsPresent (double[][] designConnectivityArray, double[] member) {
        boolean isPresent = false;
        for (int i = 0; i < designConnectivityArray.length; i++) {
            if (member[0] == designConnectivityArray[i][0]) {
                if (member[1] == designConnectivityArray[i][1]) {
                    isPresent = true;
                    break;
                }
            }
        }
        return isPresent;
    }

    private double[] addRandomRadiusMemberToArray (double[][] fullDesignConnArray, double[] memberToAdd, double[] oldRadiusArray) {
        int addPosition = 0;
        for (int i = 0; i < fullDesignConnArray.length; i++) {
            if (fullDesignConnArray[i][0] == memberToAdd[0]) {
                if (fullDesignConnArray[i][1] > memberToAdd[1])
                    break;
            }
            if (fullDesignConnArray[i][0] > memberToAdd[0])
                break;
            addPosition = i+1;
        }
        addPosition = Math.min(addPosition, (oldRadiusArray.length-1));
        oldRadiusArray[addPosition] = PRNG.nextInt(radii.length);
        return oldRadiusArray;
    }
}

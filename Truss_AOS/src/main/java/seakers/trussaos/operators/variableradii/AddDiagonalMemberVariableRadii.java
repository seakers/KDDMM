package seakers.trussaos.operators.variableradii;

import com.mathworks.engine.MatlabEngine;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import seakers.trussaos.architecture.VariableRadiiRepeatableArchitecture;

import java.util.*;
import java.util.concurrent.ExecutionException;

public class AddDiagonalMemberVariableRadii implements Variation {

    private final boolean keepFeasible;
    private static MatlabEngine engine;
    private final double[][] nodalConnectivityArray;
    private final double sidenum;
    private final double sel;
    private final double[] radiusLowerBounds;
    private final double[] radiusUpperBounds;
    private static double lowestRadiusFactor;

    public AddDiagonalMemberVariableRadii(boolean keepFeasible, double[] radiusLowerBounds, double[] radiusUpperBounds, double smallestRadiusFactor, MatlabEngine eng, double[][] nodalConnArray, double sidenum, double sel) {
        this.keepFeasible = keepFeasible;
        engine = eng;
        this.nodalConnectivityArray = nodalConnArray;
        this.sidenum = sidenum;
        this.sel = sel;
        this.radiusLowerBounds = radiusLowerBounds;
        this.radiusUpperBounds = radiusUpperBounds;
        lowestRadiusFactor = smallestRadiusFactor;
    }


    @Override
    public int getArity() {
        return 1;
    }

    @Override
    public Solution[] evolve(Solution[] sols) {
        VariableRadiiRepeatableArchitecture architecture = new VariableRadiiRepeatableArchitecture(sols[0],sidenum,radiusLowerBounds,radiusUpperBounds,lowestRadiusFactor);
        int[][] connectivityArray = architecture.getFullCAFromSolutionVariableRadii();
        int[][] completeConnectivityArray = architecture.getCompleteConnectivityArrayFromSidenum();
        double[] radiusArray = architecture.getFullRadiusArray();

        VariableRadiiRepeatableArchitecture newArchitecture = null;
        AddMemberVariableRadii addMemberObject = new AddMemberVariableRadii(keepFeasible,radiusLowerBounds,radiusUpperBounds,lowestRadiusFactor,engine,nodalConnectivityArray,sidenum,sel);
        int[][] updatedConnectivityArray;
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
                int[][] diagonalMembersAbsent = findAbsentDiagonalMembers(connectivityArray);
                if (diagonalMembersAbsent.length == 0) {
                    updatedConnectivityArray = connectivityArray.clone();
                    updatedRadiusArray = radiusArray.clone();
                }
                else if (diagonalMembersAbsent.length == 1) {
                    updatedConnectivityArray = addMemberObject.addMemberToConnectivityArray(connectivityArray,diagonalMembersAbsent[0]);
                    updatedRadiusArray = addRandomRadiusMemberToArray(completeConnectivityArray,diagonalMembersAbsent[0],radiusArray);
                }
                else {
                    Random random = new Random();
                    int diagonalMemberChoice = random.nextInt(diagonalMembersAbsent.length);
                    updatedConnectivityArray = addMemberObject.addMemberToConnectivityArray(connectivityArray,diagonalMembersAbsent[diagonalMemberChoice]);
                    updatedRadiusArray = addRandomRadiusMemberToArray(completeConnectivityArray,diagonalMembersAbsent[diagonalMemberChoice],radiusArray);
                }
                try {
                    newDesignTrussIntersections = addMemberObject.findNumberOfIntersectingTrusses(updatedConnectivityArray);
                } catch (ExecutionException | InterruptedException e) {
                    e.printStackTrace();
                }
                newArchitecture = architecture.getVariableRadiiArchitectureFromRadiusArray(updatedRadiusArray);
                if (oldDesignTrussIntersections >= newDesignTrussIntersections) {
                    feasibilitySameOrBetter = true;
                }
                else {
                    numberOfAttempts += 1;
                }
            }
        }
        else {
            int[][] diagonalMembersAbsent = findAbsentDiagonalMembers(connectivityArray);
            if (diagonalMembersAbsent.length == 0) {
                //updatedConnectivityArray = connectivityArray.clone();
                updatedRadiusArray = radiusArray.clone();
            }
            else if (diagonalMembersAbsent.length == 1) {
                //updatedConnectivityArray = addMemberObject.addMemberToConnectivityArray(connectivityArray,diagonalMembersAbsent[0]);
                updatedRadiusArray = addRandomRadiusMemberToArray(completeConnectivityArray,diagonalMembersAbsent[0],radiusArray);
            }
            else {
                Random random = new Random();
                int diagonalMemberChoice = random.nextInt(diagonalMembersAbsent.length);
                //updatedConnectivityArray = addMemberObject.addMemberToConnectivityArray(connectivityArray,diagonalMembersAbsent[diagonalMemberChoice]);
                updatedRadiusArray = addRandomRadiusMemberToArray(completeConnectivityArray,diagonalMembersAbsent[diagonalMemberChoice],radiusArray);
            }
            newArchitecture = architecture.getVariableRadiiArchitectureFromRadiusArray(updatedRadiusArray);
        }
        return new Solution[]{newArchitecture};
    }

    private int[] getCornerNodes () {
        // Array order: {bottom left, top left, bottom right, top right}
        return new int[]{1, (int) sidenum, (int) (sidenum * sidenum - (sidenum - 1)), (int) (sidenum * sidenum)};
    }

    private int[] getCentreNodes () {
        // Array order: {left, down, top, right} (assuming only odd values for sidenum)
        return new int[]{(int) (sidenum - ((sidenum - 1)/2)), (int) ((((sidenum - 1)/2) * sidenum) + 1), (int) (sidenum * (sidenum + 1)/2), (int) ((sidenum * sidenum) - ((sidenum - 1)/2))};

    }

    private int[][] getAllDiagonalMembers () {
        ArrayList<ArrayList<Integer>> diagonalMembers = new ArrayList<ArrayList<Integer>>();

        for (int i = 0; i < ((sidenum * sidenum) - 1); i++) {
            int nodeNumber = i+1;
            int closestTopNode = findClosestTopNode(nodeNumber);
            for (int j = nodeNumber+1; j < ((sidenum * sidenum)+1); j++) {
                if ((j <= closestTopNode) || ((j - nodeNumber)%sidenum == 0)) {
                    continue;
                } else {
                    ArrayList<Integer> currentDiagonalMember = new ArrayList<Integer>();
                    currentDiagonalMember.add(nodeNumber);
                    currentDiagonalMember.add(j);
                    diagonalMembers.add(currentDiagonalMember);
                }
            }
        }
        return convertArrayListTo2DArray(diagonalMembers);
    }

    private int findClosestTopNode (int node) {
        boolean topReached = false;
        int currentNode = node;
        while (!topReached) {
            if (currentNode%sidenum == 0) {
                topReached = true;
            } else {
                currentNode += 1;
            }
        }
        return currentNode;
    }

    private int[][] convertArrayListTo2DArray (ArrayList<ArrayList<Integer>> arrayList2D) {
        int[][] array2D = new int[arrayList2D.size()][2];
        for (int i = 0; i < arrayList2D.size(); i++) {
            ArrayList<Integer> row = arrayList2D.get(i);
            int[] rowArray = row.stream().mapToInt(j->j).toArray();
            array2D[i] = rowArray;
        }
        return array2D;
    }

    private int[][] findAbsentDiagonalMembers (int[][] designConnectivityArray) {
        int[][] allDiagonalMembers = getAllDiagonalMembers();
        ArrayList<ArrayList<Integer>> absentDiagonalMembers = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < allDiagonalMembers.length; i++) {
            boolean presence = checkIfMemberIsPresent(designConnectivityArray, allDiagonalMembers[i]);
            if (!presence) {
                ArrayList<Integer> absentDiagonalMember = new ArrayList<>();
                absentDiagonalMember.add(allDiagonalMembers[i][0]);
                absentDiagonalMember.add(allDiagonalMembers[i][1]);
                absentDiagonalMembers.add(absentDiagonalMember);
            }
        }
        return convertArrayListTo2DArray(absentDiagonalMembers);
    }

    private boolean checkIfMemberIsPresent (int[][] designConnectivityArray, int[] member) {
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

    private double[] addRandomRadiusMemberToArray (int[][] fullDesignConnArray, int[] memberToAdd, double[] oldRadiusArray) {
        synchronized (PRNG.getRandom()) {
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
            oldRadiusArray[addPosition] = PRNG.nextDouble(radiusLowerBounds[addPosition],radiusUpperBounds[addPosition]);
            return oldRadiusArray;
        }
    }

}

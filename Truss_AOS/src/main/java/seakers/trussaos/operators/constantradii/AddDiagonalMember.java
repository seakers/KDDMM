package seakers.trussaos.operators.constantradii;

import com.mathworks.engine.MatlabEngine;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;

import java.util.*;
import java.util.concurrent.ExecutionException;

public class AddDiagonalMember implements Variation {

    private final boolean keepFeasible;
    private static MatlabEngine engine;
    private final double[][] nodalConnectivityArray;
    private final double sidenum;
    private final double sel;
    private final int numHeurObjectives;
    private final int numHeurConstraints;

    public AddDiagonalMember(boolean keepFeasible, MatlabEngine eng, double[][] nodalConnArray, double sidenum, double sel, int numHeuristicObjectives, int numHeuristicConstraints) {
        this.keepFeasible = keepFeasible;
        engine = eng;
        this.nodalConnectivityArray = nodalConnArray;
        this.sidenum = sidenum;
        this.sel = sel;
        this.numHeurObjectives = numHeuristicObjectives;
        this.numHeurConstraints = numHeuristicConstraints;
    }


    @Override
    public int getArity() {
        return 1;
    }

    @Override
    public Solution[] evolve(Solution[] sols) {
        TrussRepeatableArchitecture architecture = new TrussRepeatableArchitecture(sols[0],sidenum,numHeurObjectives,numHeurConstraints);
        double[][] connectivityArray = architecture.getConnectivityArrayFromSolution(sols[0]);
        TrussRepeatableArchitecture newArchitecture = null;
        AddMember addMemberObject = new AddMember(keepFeasible,engine,nodalConnectivityArray,sidenum,sel,numHeurObjectives,numHeurConstraints);
        double[][] updatedConnectivityArray;

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
                }
                else if (diagonalMembersAbsent.length == 1) {
                    updatedConnectivityArray = addMemberObject.addMemberToConnectivityArray(connectivityArray,diagonalMembersAbsent[0]);
                    newArchitecture = architecture.getArchitectureFromConnectivityArray(updatedConnectivityArray);
                }
                else {
                    Random random = new Random();
                    int diagonalMemberChoice = random.nextInt(diagonalMembersAbsent.length);
                    updatedConnectivityArray = addMemberObject.addMemberToConnectivityArray(connectivityArray,diagonalMembersAbsent[diagonalMemberChoice]);
                }
                try {
                    newDesignTrussIntersections = addMemberObject.findNumberOfIntersectingTrusses(updatedConnectivityArray);
                } catch (ExecutionException | InterruptedException e) {
                    e.printStackTrace();
                }
                //int newNumberOfIntersectingTrusses = newDesignTrussIntersections.get(0) == null ? 0 : newDesignTrussIntersections.get(0).size();
                //int oldNumberOfIntersectingTrusses = oldDesignTrussIntersections.get(0) == null ? 0 : oldDesignTrussIntersections.get(0).size();
                if (oldDesignTrussIntersections >= newDesignTrussIntersections) {
                    feasibilitySameOrBetter = true;
                    newArchitecture = architecture.getArchitectureFromConnectivityArray(updatedConnectivityArray);
                }
                else {
                    numberOfAttempts += 1;
                }
            }
        }
        else {
            double[][] diagonalMembersAbsent = findAbsentDiagonalMembers(connectivityArray);
            if (diagonalMembersAbsent.length == 0) {
                updatedConnectivityArray = connectivityArray.clone();
            }
            else if (diagonalMembersAbsent.length == 1) {
                updatedConnectivityArray = addMemberObject.addMemberToConnectivityArray(connectivityArray,diagonalMembersAbsent[0]);
            }
            else {
                Random random = new Random();
                int diagonalMemberChoice = random.nextInt(diagonalMembersAbsent.length);
                updatedConnectivityArray = addMemberObject.addMemberToConnectivityArray(connectivityArray,diagonalMembersAbsent[diagonalMemberChoice]);
            }
            newArchitecture = architecture.getArchitectureFromConnectivityArray(updatedConnectivityArray);
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

    private double[][] getAllDiagonalMembers () {
        ArrayList<ArrayList<Double>> diagonalMembers = new ArrayList<ArrayList<Double>>();

        for (int i = 0; i < ((sidenum * sidenum) - 1); i++) {
            double nodeNumber = i+1;
            double closestTopNode = findClosestTopNode(nodeNumber);
            for (double j = nodeNumber+1; j < ((sidenum * sidenum)+1); j++) {
                if ((j <= closestTopNode) || ((j - nodeNumber)%sidenum == 0)) {
                    continue;
                } else {
                    ArrayList<Double> currentDiagonalMember = new ArrayList<Double>();
                    currentDiagonalMember.add(nodeNumber);
                    currentDiagonalMember.add(j);
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

}

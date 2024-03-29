package seakers.trussaos.operators.constantradii;

import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.core.PRNG;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;

import java.lang.Math;
import java.util.ArrayList;

public class ImproveOrientation implements Variation {

    private final boolean arteryProblem;
    private final double[][] nodalConnectivityArray;
    private final int sidenum;
    private final double targetCRatio;
    private final int numHeurObjectives;
    private final int numHeurConstraints;

    public ImproveOrientation(boolean arteryProblem, double[][] nodalConnArray, double targetCRatio, int sidenum, int numHeuristicObjectives, int numHeuristicConstraints) {
        this.arteryProblem = arteryProblem;
        this.nodalConnectivityArray = nodalConnArray;
        this.targetCRatio = targetCRatio;
        this.sidenum = sidenum;
        this.numHeurObjectives = numHeuristicObjectives;
        this.numHeurConstraints = numHeuristicConstraints;
    }

    @Override
    public int getArity() {
        return 1;
    }

    @Override
    public Solution[] evolve(Solution[] solutions) {
        TrussRepeatableArchitecture architecture = new TrussRepeatableArchitecture(solutions[0],sidenum,numHeurObjectives,numHeurConstraints);
        double[][] connectivityArray = architecture.getConnectivityArrayFromSolution(solutions[0]);
        ArrayList<ArrayList<Double>> connectivityArrayList = new ArrayList<ArrayList<Double>>();
        for (double[] member:connectivityArray) {
            ArrayList<Double> memberList = new ArrayList<Double>();
            for (double node:member) {
                memberList.add(node);
            }
            connectivityArrayList.add(memberList);
        }

        TrussRepeatableArchitecture newArchitecture = null;

        double[] architectureOrientations = findMemberOrientations(connectivityArray);

        double totalOrientation = 0;
        for (int i = 0; i < architectureOrientations.length; i++) {
            totalOrientation += architectureOrientations[i];
        }
        double meanOrientation = totalOrientation/architectureOrientations.length;

        double targetOrientation = findTargetOrientation(targetCRatio);

        int[] horizontalCornerNodes = new int[sidenum];
        int[] verticalCornerNodes = new int[sidenum];

        for (int j = 0; j < sidenum; j++) {
            horizontalCornerNodes[j] = sidenum*j + 1;
            verticalCornerNodes[j] = j+1;
        }

        ArrayList<ArrayList<Integer>> verticalMembers = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> horizontalMembers = new ArrayList<ArrayList<Integer>>();

        for (int k:horizontalCornerNodes) {
            for (int m = 0; m < sidenum-1; m++) {
                ArrayList<Integer> verticalMember = new ArrayList<Integer>();
                verticalMember.add(k+m);
                verticalMember.add(k+m+1);
                verticalMembers.add(verticalMember);
            }
        }

        for (int n:verticalCornerNodes) {
            for (int p = 0; p < sidenum-1; p++) {
                ArrayList<Integer> horizontalMember = new ArrayList<Integer>();
                horizontalMember.add(n+sidenum*p);
                horizontalMember.add(n+sidenum*(p+1));
                horizontalMembers.add(horizontalMember);
            }
        }

        boolean memberAdded = false;
        int numAttempts = 0;
        double[][] newConnArray = null;
        if (meanOrientation > targetOrientation) {
            while (!memberAdded && (numAttempts < 5)) {
                ArrayList<Integer> memberChoice = new ArrayList<Integer>();
                int randomMemberIndex = PRNG.nextInt(horizontalMembers.size());
                memberChoice = horizontalMembers.get(randomMemberIndex);
                if (!connectivityArrayList.contains(memberChoice)) {
                    newConnArray = addMemberToConnectivityArray(connectivityArray,memberChoice.stream().filter(t -> t != null).mapToDouble(t -> t).toArray());
                    memberAdded = true;
                }
                else {
                    numAttempts += 1;
                }
            }
            if (numAttempts == 5) {
                newConnArray = connectivityArray.clone();
            }
        }
        else if (meanOrientation < targetOrientation) {
            while (!memberAdded && (numAttempts < 5)) {
                ArrayList<Integer> memberChoice = new ArrayList<Integer>();
                int randomMemberIndex = PRNG.nextInt(verticalMembers.size());
                memberChoice = verticalMembers.get(randomMemberIndex);
                if (!connectivityArrayList.contains(memberChoice)) {
                    newConnArray = addMemberToConnectivityArray(connectivityArray,memberChoice.stream().filter(t -> t != null).mapToDouble(t -> t).toArray());
                    memberAdded = true;
                }
                else {
                    numAttempts += 1;
                }
            }
            if (numAttempts == 5) {
                newConnArray = connectivityArray.clone();
            }
        }
        else if (meanOrientation == targetOrientation) {
            newConnArray = connectivityArray.clone();
        }
        assert newConnArray != null;
        newArchitecture = architecture.getArchitectureFromConnectivityArray(newConnArray, arteryProblem);
        return new Solution[]{newArchitecture};
    }

    public double[] findMemberOrientations(double[][] designConnArray) {
        double[] memberOrientations = new double[designConnArray.length];
        for (int i = 0; i < designConnArray.length; i++) {
            memberOrientations[i] = findMemberOrientation(designConnArray[i]);
        }
        return memberOrientations;
    }

    private double findMemberOrientation(double[] memberNodes) { // CHANGE BASED ON NEW ORIENTATION HEURISTIC
        double x1 = nodalConnectivityArray[(int) memberNodes[0]-1][0];
        double y1 = nodalConnectivityArray[(int) memberNodes[0]-1][1];
        double x2 = nodalConnectivityArray[(int) memberNodes[1]-1][0];
        double y2 = nodalConnectivityArray[(int) memberNodes[1]-1][1];

        double memberLength = Math.sqrt(((x2-x1)*(x2-x1)) + ((y2-y1)*(y2-y1)));

        return Math.toDegrees(Math.acos((x2-x1)/memberLength));
    }

    public double findTargetOrientation(double targetRatio) { // MUST BE UPDATED
        return Math.toDegrees((0.5d*Math.atan(targetRatio-1)) + (Math.PI/4));
    }

    private double[][] addMemberToConnectivityArray (double[][] oldConnectivityArray, double[] trussToAdd) {
        int addPosition = 0;
        //boolean positionReached = false;
        for (int i = 0; i < oldConnectivityArray.length; i++) {
            if (oldConnectivityArray[i][0] == trussToAdd[0]) {
                //if (oldConnectivityArray[i][1] < trussToAdd[1])
                //addPosition += 1;
                if (oldConnectivityArray[i][1] > trussToAdd[1])
                    //positionReached = true;
                    break;
            }
            if (oldConnectivityArray[i][0] > trussToAdd[0])
                break;
            //if (positionReached)
            //break;
            addPosition = i+1;
        }
        double[][] newConnectivityArray = new double[oldConnectivityArray.length + 1][2];
        int currentPosition = 0;
        for (int i = 0; i < newConnectivityArray.length; i++) {
            if (i == addPosition) {
                for (int j = 0; j < 2; j++) {
                    newConnectivityArray[i][j] = trussToAdd[j];
                }
            }
            else {
                for (int j = 0; j < 2; j++) {
                    newConnectivityArray[i][j] = oldConnectivityArray[currentPosition][j];
                }
                currentPosition += 1;
            }
        }
        return newConnectivityArray;
    }
}

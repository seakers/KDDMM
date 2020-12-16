package seakers.trussaos.operators;

import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;

import java.lang.Math;
import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;

public class ImproveOrientation implements Variation {

    private final double[][] nodalConnectivityArray;

    private final int sidenum;

    private final double targetCRatio;

    public ImproveOrientation(double[][] nodalConnArray, double targetCRatio, int sidenum) {
        this.nodalConnectivityArray = nodalConnArray;
        this.targetCRatio = targetCRatio;
        this.sidenum = sidenum;
    }

    @Override
    public int getArity() {
        return 1;
    }

    @Override
    public Solution[] evolve(Solution[] solutions) {
        TrussRepeatableArchitecture architecture = new TrussRepeatableArchitecture(solutions[0]);
        int[][] connectivityArray = architecture.getConnectivityArrayFromSolution(solutions[0]);
        ArrayList<ArrayList<Integer>> connectivityArrayList = new ArrayList<ArrayList<Integer>>();
        for (int[] member:connectivityArray) {
            ArrayList<Integer> memberList = new ArrayList<Integer>();
            for (int node:member) {
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
        int[][] newConnArray = null;
        if (meanOrientation > targetOrientation) {
            while (!memberAdded && (numAttempts < 5)) {
                ArrayList<Integer> memberChoice = new ArrayList<Integer>();
                int randomMemberIndex = ThreadLocalRandom.current().nextInt(0,horizontalMembers.size());
                memberChoice = horizontalMembers.get(randomMemberIndex);
                if (!connectivityArrayList.contains(memberChoice)) {
                    newConnArray = addMemberToConnectivityArray(connectivityArray,memberChoice.stream().filter(t -> t != null).mapToInt(t -> t).toArray());
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
                int randomMemberIndex = ThreadLocalRandom.current().nextInt(0,verticalMembers.size());
                memberChoice = verticalMembers.get(randomMemberIndex);
                if (!connectivityArrayList.contains(memberChoice)) {
                    newConnArray = addMemberToConnectivityArray(connectivityArray,memberChoice.stream().filter(t -> t != null).mapToInt(t -> t).toArray());
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
        newArchitecture = architecture.getArchitectureFromConnectivityArray(newConnArray);
        return new Solution[]{newArchitecture};
    }

    public double[] findMemberOrientations(int[][] designConnArray) {
        double[] memberOrientations = new double[designConnArray.length];
        for (int i = 0; i < designConnArray.length; i++) {
            memberOrientations[i] = findMemberOrientation(designConnArray[i]);
        }
        return memberOrientations;
    }

    private double findMemberOrientation(int[] memberNodes) {
        double x1 = nodalConnectivityArray[memberNodes[0]-1][0];
        double y1 = nodalConnectivityArray[memberNodes[0]-1][1];
        double x2 = nodalConnectivityArray[memberNodes[1]-1][0];
        double y2 = nodalConnectivityArray[memberNodes[1]-1][1];

        double memberLength = Math.sqrt(((x2-x1)*(x2-x1)) + ((y2-y1)*(y2-y1)));

        return Math.toDegrees(Math.acos((x2-x1)/memberLength));
    }

    public double findTargetOrientation(double targetRatio) {
        return Math.toDegrees((0.5d*Math.atan(targetRatio-1)) + (Math.PI/4));
    }

    private int[][] addMemberToConnectivityArray (int[][] oldConnectivityArray, int[] trussToAdd) {
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
        int[][] newConnectivityArray = new int[oldConnectivityArray.length + 1][2];
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

package seakers.trussaos.operators.integerradii;

import org.moeaframework.core.PRNG;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import seakers.trussaos.architecture.IntegerRepeatableArchitecture;

import java.lang.Math;
import java.util.ArrayList;

public class ImproveOrientationIntegerRadii implements Variation {

    private final double[][] nodalConnectivityArray;
    private final double sidenum;
    private final double targetCRatio;
    private final double[] radii;
    private final double[] youngsModulii;
    private final int numHeurObjectives;
    private final int numHeurConstraints;

    public ImproveOrientationIntegerRadii(double[][] nodalConnArray, double[] radii, double[] youngsModulii, double targetCRatio, double sidenum, int numHeurObjectives, int numHeurConstraints) {
        this.nodalConnectivityArray = nodalConnArray;
        this.targetCRatio = targetCRatio;
        this.sidenum = sidenum;
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
        int youngsModulusChoice = architecture.getYoungsModulusChoice();

        ArrayList<ArrayList<Double>> connectivityArrayList = new ArrayList<ArrayList<Double>>();
        for (double[] member:connectivityArray) {
            ArrayList<Double> memberList = new ArrayList<Double>();
            for (double node:member) {
                memberList.add(node);
            }
            connectivityArrayList.add(memberList);
        }

        IntegerRepeatableArchitecture newArchitecture = null;

        double[] architectureOrientations = findMemberOrientations(connectivityArray);

        double totalOrientation = 0;
        for (int i = 0; i < architectureOrientations.length; i++) {
            totalOrientation += architectureOrientations[i];
        }
        double meanOrientation = totalOrientation/architectureOrientations.length;

        double targetOrientation = findTargetOrientation(targetCRatio);

        int[] horizontalCornerNodes = new int[(int) sidenum];
        int[] verticalCornerNodes = new int[(int) sidenum];

        for (int j = 0; j < sidenum; j++) {
            horizontalCornerNodes[j] = (int) (sidenum*j + 1);
            verticalCornerNodes[j] = j+1;
        }

        ArrayList<ArrayList<Double>> verticalMembers = new ArrayList<ArrayList<Double>>();
        ArrayList<ArrayList<Double>> horizontalMembers = new ArrayList<ArrayList<Double>>();

        for (int k:horizontalCornerNodes) {
            for (int m = 0; m < sidenum-1; m++) {
                ArrayList<Double> verticalMember = new ArrayList<Double>();
                verticalMember.add((double) (k+m));
                verticalMember.add((double) (k+m+1));
                verticalMembers.add(verticalMember);
            }
        }

        for (int n:verticalCornerNodes) {
            for (int p = 0; p < sidenum-1; p++) {
                ArrayList<Double> horizontalMember = new ArrayList<Double>();
                horizontalMember.add(n+sidenum*p);
                horizontalMember.add(n+sidenum*(p+1));
                horizontalMembers.add(horizontalMember);
            }
        }

        boolean memberAdded = false;
        int numAttempts = 0;
        //int[][] newConnArray = null;
        double[] newRadiusArray = null;
        if (meanOrientation > targetOrientation) {
            while (!memberAdded && (numAttempts < 5)) {
                ArrayList<Double> memberChoice = new ArrayList<Double>();
                int randomMemberIndex = PRNG.nextInt(0,horizontalMembers.size());
                memberChoice = horizontalMembers.get(randomMemberIndex);
                if (!connectivityArrayList.contains(memberChoice)) {
                    //newConnArray = addMemberToConnectivityArray(connectivityArray,memberChoice.stream().filter(t -> t != null).mapToInt(t -> t).toArray());
                    newRadiusArray = addRandomRadiusMemberToArray(completeConnectivityArray,memberChoice.stream().filter(t -> t != null).mapToDouble(t -> t).toArray(),radiusArray);
                    memberAdded = true;
                }
                else {
                    numAttempts += 1;
                }
            }
            if (numAttempts == 5) {
                //newConnArray = connectivityArray.clone();
                newRadiusArray = radiusArray.clone();
            }
        }
        else if (meanOrientation < targetOrientation) {
            while (!memberAdded && (numAttempts < 5)) {
                ArrayList<Double> memberChoice = new ArrayList<Double>();
                int randomMemberIndex = PRNG.nextInt(0,verticalMembers.size());
                memberChoice = verticalMembers.get(randomMemberIndex);
                if (!connectivityArrayList.contains(memberChoice)) {
                    //newConnArray = addMemberToConnectivityArray(connectivityArray,memberChoice.stream().filter(t -> t != null).mapToInt(t -> t).toArray());
                    newRadiusArray = addRandomRadiusMemberToArray(completeConnectivityArray,memberChoice.stream().filter(t -> t != null).mapToDouble(t -> t).toArray(),radiusArray);
                    memberAdded = true;
                }
                else {
                    numAttempts += 1;
                }
            }
            if (numAttempts == 5) {
                //newConnArray = connectivityArray.clone();
                newRadiusArray = radiusArray.clone();
            }
        }
        else if (meanOrientation == targetOrientation) {
            //newConnArray = connectivityArray.clone();
            newRadiusArray = radiusArray.clone();
        }
        //assert newConnArray != null;
        newArchitecture = architecture.getIntegerArchitectureFromRadiusArray(newRadiusArray, youngsModulusChoice);
        return new Solution[]{newArchitecture};
    }

    public double[] findMemberOrientations(double[][] designConnArray) {
        double[] memberOrientations = new double[designConnArray.length];
        for (int i = 0; i < designConnArray.length; i++) {
            memberOrientations[i] = findMemberOrientation(designConnArray[i]);
        }
        return memberOrientations;
    }

    private double findMemberOrientation(double[] memberNodes) {
        double x1 = nodalConnectivityArray[(int) (memberNodes[0]-1)][0];
        double y1 = nodalConnectivityArray[(int) (memberNodes[0]-1)][1];
        double x2 = nodalConnectivityArray[(int) (memberNodes[1]-1)][0];
        double y2 = nodalConnectivityArray[(int) (memberNodes[1]-1)][1];

        double memberLength = Math.sqrt(((x2-x1)*(x2-x1)) + ((y2-y1)*(y2-y1)));

        return Math.toDegrees(Math.acos((x2-x1)/memberLength));
    }

    public double findTargetOrientation(double targetRatio) {
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

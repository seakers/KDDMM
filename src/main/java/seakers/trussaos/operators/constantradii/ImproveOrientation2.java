package seakers.trussaos.operators.constantradii;

import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.core.PRNG;
import seakers.trussaos.architecture.TrussRepeatableArchitecture;

import java.lang.Math;
import java.util.*;

/**
 * This version of the orientation operator is greedier than the previous "ImproveOrientation" class. Here, a random member with the
 * angle closest to the difference between the mean design angle and target angle is added to the design.
 */

public class ImproveOrientation2 implements Variation {

    private final boolean arteryProblem;
    private final double[][] nodalConnectivityArray;
    private final int sidenum;
    private final double targetCRatio;
    private final int numHeurObjectives;
    private final int numHeurConstraints;

    public ImproveOrientation2(boolean arteryProblem, double[][] nodalConnArray, double targetCRatio, int sidenum, int numHeuristicObjectives, int numHeuristicConstraints) {
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

        // Find average member angle of design
        double[] memberOrientations = new double[connectivityArray.length];
        for (int i = 0; i < connectivityArray.length; i++) {
            memberOrientations[i] = findMemberOrientation2(connectivityArray[i]);
        }
        double totalOrientation = 0;
        for (int i = 0; i < memberOrientations.length; i++) {
            totalOrientation += memberOrientations[i];
        }
        double meanOrientation = totalOrientation/memberOrientations.length;

        // Find target orientation based on target stiffness ratio
        double targetOrientation = Math.toDegrees(Math.atan(targetCRatio));

        double orientationDifference = meanOrientation - targetOrientation;

        if (orientationDifference == 0) {
            newArchitecture = architecture.getArchitectureFromConnectivityArray(connectivityArray, arteryProblem);
        } else {
            // Store the candidate members with their orientations (candidate members are those that are not present in the current design)
            HashMap<ArrayList<Double>, Double> candidateMemberMeanOrientations = new HashMap<>();
            for (int i = 1; i < sidenum*sidenum; i++) {
                for (int j = i+1; j <= sidenum*sidenum; j++) {
                    ArrayList<Double> currentMember = new ArrayList<>();
                    currentMember.add((double) i);
                    currentMember.add((double) j);
                    if (memberPresentInDesign(currentMember, connectivityArrayList)) {
                        continue;
                    } else {
                        double currentMemberMeanOrientation = findDesignOrientationWithMember(totalOrientation, memberOrientations.length, currentMember.stream().mapToDouble(Double::doubleValue).toArray());
                        candidateMemberMeanOrientations.put(currentMember, currentMemberMeanOrientation);
                    }
                }
            }

            // Choose and add a member based on the difference between the mean and target orientations
            ArrayList<Double> candidateMeanOrientations = new ArrayList<Double>(candidateMemberMeanOrientations.values());
            ArrayList<Double> memberToAdd = new ArrayList<>();
            double[] candidateMeanOrientationDifferences = new double[candidateMeanOrientations.size()];
            for (int i = 0; i < candidateMeanOrientations.size(); i++) {
                candidateMeanOrientationDifferences[i] = Math.abs(candidateMeanOrientations.get(i) - targetOrientation);
            }

            int indexOfMinimumDifference = 0;
            double minimumAbsoluteDifference = candidateMeanOrientationDifferences[0];
            for (int i = 0; i < candidateMeanOrientationDifferences.length; i++) {
                if (candidateMeanOrientationDifferences[i] < minimumAbsoluteDifference) {
                    minimumAbsoluteDifference = candidateMeanOrientationDifferences[i];
                    indexOfMinimumDifference = i;
                }
            }
            ArrayList<ArrayList<Double>> membersWithMinimumDifference =  getMembersForOrientationDifference(candidateMemberMeanOrientations, candidateMeanOrientations.get(indexOfMinimumDifference));
            memberToAdd = membersWithMinimumDifference.get(PRNG.nextInt(membersWithMinimumDifference.size()));
            if (findDesignOrientationWithMember(totalOrientation, memberOrientations.length, memberToAdd.stream().mapToDouble(Double::doubleValue).toArray()) > meanOrientation) {
                newArchitecture = architecture.getArchitectureFromConnectivityArray(connectivityArray, arteryProblem);
            } else {
                double[][] newConnectivityArray = addMemberToConnectivityArray(connectivityArray, memberToAdd.stream().mapToDouble(d->d).toArray());
                newArchitecture = architecture.getArchitectureFromConnectivityArray(newConnectivityArray, arteryProblem);
            }
        }

        // Store the candidate members with their orientations (candidate members are those that are not present in the current design)
        // TWO SEPARATE HASHMAPS FOR AVERAGE ORIENTATIONS < & > THAN TARGET ORIENTATION
        //HashMap<ArrayList<Double>, Double> candidateMemberLowerOrientations = new HashMap<>();
        //HashMap<ArrayList<Double>, Double> candidateMemberHigherOrientations = new HashMap<>();
        //for (int i = 1; i < sidenum*sidenum; i++) {
            //for (int j = i+1; j <= sidenum*sidenum; j++) {
                //ArrayList<Double> currentMember = new ArrayList<>();
                //currentMember.add((double) i);
                //currentMember.add((double) j);
                //if (memberPresentInDesign(currentMember, connectivityArrayList)) {
                    //continue;
                //} else {
                    ////double currentMemberMeanOrientation = findDesignOrientationWithMember(totalOrientation, memberOrientations.length, currentMember.stream().mapToDouble(Double::doubleValue).toArray());
                    //double currentMemberMeanOrientation = findMemberOrientation2(currentMember.stream().mapToDouble(Double::doubleValue).toArray());
                    //if (currentMemberMeanOrientation < targetOrientation) {
                        //candidateMemberLowerOrientations.put(currentMember, currentMemberMeanOrientation);
                    //} else {
                        //candidateMemberHigherOrientations.put(currentMember, currentMemberMeanOrientation);
                    //}
                //}
            //}
        //}

        // Choose and add a member based on the difference between the mean and target orientations
        // IF (ORIENTATIONDIFFERENCE < 0 ) { USE HIGHER ORIENTATION HASHMAP } ELSE { USE LOWER ORIENTATION HASHMAP }
        //ArrayList<ArrayList<Double>> candidatesLowerOrientations = new ArrayList<ArrayList<Double>>(candidateMemberLowerOrientations.keySet());
        //ArrayList<ArrayList<Double>> candidatesHigherOrientations = new ArrayList<ArrayList<Double>>(candidateMemberHigherOrientations.keySet());
        //ArrayList<Double> memberToAdd = new ArrayList<>();

        //if (orientationDifference < 0) {
            //double[] candidateOrientationDifferences = new double[candidatesHigherOrientations.size()];
            //for (int i = 0; i < candidatesHigherOrientations.size(); i++) {
                //double currentCandidateMeanOrientation = findDesignOrientationWithMember(totalOrientation, memberOrientations.length, candidatesHigherOrientations.get(i).stream().mapToDouble(Double::doubleValue).toArray());
                //candidateOrientationDifferences[i] = Math.abs(currentCandidateMeanOrientation - targetOrientation);
            //}

            //int indexOfMinimumDifference = 0;
            //double minimumAbsoluteDifference = candidateOrientationDifferences[0];
            //for (int i = 0; i < candidateOrientationDifferences.length; i++) {
                //if (candidateOrientationDifferences[i] < minimumAbsoluteDifference) {
                    //minimumAbsoluteDifference = candidateOrientationDifferences[i];
                    //indexOfMinimumDifference = i;
                //}
            //}
            //// CHANGE LINE BELOW USE CANDIDATELOWERORIENTATIONS/CANDIDATEHIGHERORIENTATIONS AND CANDIDATEORIENTATIONDIFFERENCES
            //ArrayList<ArrayList<Double>> membersWithMinimumDifference =  getMembersForOrientationDifference(candidateMemberHigherOrientations, minimumAbsoluteDifference);
            //memberToAdd = membersWithMinimumDifference.get(PRNG.nextInt(membersWithMinimumDifference.size()));
        //} else {
            //double[] candidateOrientationDifferences = new double[candidatesLowerOrientations.size()];
            //for (int i = 0; i < candidatesLowerOrientations.size(); i++) {
                //double currentCandidateMeanOrientation = findDesignOrientationWithMember(totalOrientation, memberOrientations.length, candidatesLowerOrientations.get(i).stream().mapToDouble(Double::doubleValue).toArray());
                //candidateOrientationDifferences[i] = Math.abs(currentCandidateMeanOrientation - targetOrientation);
            //}

            //int indexOfMinimumDifference = 0;
            //double minimumAbsoluteDifference = candidateOrientationDifferences[0];
            //for (int i = 0; i < candidateOrientationDifferences.length; i++) {
                //if (candidateOrientationDifferences[i] < minimumAbsoluteDifference) {
                    //minimumAbsoluteDifference = candidateOrientationDifferences[i];
                    //indexOfMinimumDifference = i;
                //}
            //}
            //// CHANGE LINE BELOW USE CANDIDATELOWERORIENTATIONS/CANDIDATEHIGHERORIENTATIONS AND CANDIDATEORIENTATIONDIFFERENCES
            //ArrayList<ArrayList<Double>> membersWithMinimumDifference =  getMembersForOrientationDifference(candidateMemberLowerOrientations, minimumAbsoluteDifference);
            //memberToAdd = membersWithMinimumDifference.get(PRNG.nextInt(membersWithMinimumDifference.size()));
        //}

        return new Solution[]{newArchitecture};
    }

    private double findMemberOrientation2(double[] memberNodes) {
        double x1 = nodalConnectivityArray[(int) memberNodes[0]-1][0];
        double y1 = nodalConnectivityArray[(int) memberNodes[0]-1][1];
        double x2 = nodalConnectivityArray[(int) memberNodes[1]-1][0];
        double y2 = nodalConnectivityArray[(int) memberNodes[1]-1][1];

        return Math.abs(Math.toDegrees(Math.atan((y2 - y1)/(x2 - x1))));
    }

    private double findDesignOrientationWithMember(double totalConnectivityArrayOrientation, double numberOfMembers, double[] addedMember) {
        double memberOrientation = findMemberOrientation2(addedMember);
        return (totalConnectivityArrayOrientation + memberOrientation)/(numberOfMembers + 1);
    }

    private boolean memberPresentInDesign(ArrayList<Double> member, ArrayList<ArrayList<Double>> designConnectivityArrayList) {
        boolean isPresent = false;
        for (ArrayList<Double> designMember : designConnectivityArrayList) {
            if (member.equals(designMember)) {
                isPresent = true;
                break;
            }
        }
        return isPresent;
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

    private ArrayList<ArrayList<Double>> getMembersForOrientationDifference(HashMap<ArrayList<Double>, Double> candidateMembersMap, double orientationValue) {
        ArrayList<ArrayList<Double>> candidateMembers = new ArrayList<>();
        if (candidateMembersMap.containsValue(orientationValue)) {
            for (Map.Entry<ArrayList<Double>, Double> entry : candidateMembersMap.entrySet()) {
                if (Objects.equals(entry.getValue(), orientationValue)) {
                    candidateMembers.add(entry.getKey());
                }
            }
        }
        return candidateMembers;
    }
}

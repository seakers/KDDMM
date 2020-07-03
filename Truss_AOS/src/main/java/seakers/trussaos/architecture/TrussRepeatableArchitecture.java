package seakers.trussaos.architecture;

import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.EncodingUtils;

/**
 * Design class for the binary design vector in the 2D 3x3 nodal grid case
 *
 * @author roshan94
 */

public class TrussRepeatableArchitecture extends Solution{
    private static final long serialVersionUID = -3246610489443538032L;

    private int NumberOfTrusses;

    private boolean[] CompleteBooleanDesign;

    private int[][] FullConnectivityArray = {{1,2}, {1,3}, {1,4}, {1,5}, {1,6}, {1,7}, {1,8}, {1,9},
                                            {2,3}, {2,4}, {2,5}, {2,6}, {2,7}, {2,8}, {2,9},
                                            {3,4}, {3,5}, {3,6}, {3,7}, {3,8}, {3,9},
                                            {4,5}, {4,6}, {4,7}, {4,8}, {4,9},
                                            {5,6}, {5,7}, {5,8}, {5,9},
                                            {6,7}, {6,8}, {6,9},
                                            {7,8}, {7,9},
                                            {8,9}};

    private int[][] DesignConnectivityArray;

    public TrussRepeatableArchitecture(Solution solution) {
        super(solution);
        // Solution sol = (TrussRepeatableArchitecture) solution;

        //// Extract the design as a boolean array
        boolean[] BooleanDesign = getBooleanDesignArray(solution);
        //design = EncodingUtils.getBinary(solution.getVariable());

        //// Obtain corresponding Connectivity Array
        DesignConnectivityArray= ConvertToFullConnectivityArray(BooleanDesign);

        //// Compute number of trusses in the design
        CompleteBooleanDesign = getCompleteBooleanDesignArray(BooleanDesign);
        NumberOfTrusses = getTrussCount(CompleteBooleanDesign);

    }

    private int[][] ConvertToFullConnectivityArray(boolean[] CurrentDesign){
        boolean[] FullDesign = getCompleteBooleanDesignArray(CurrentDesign);
        int TrussCount = 0;
        int[][] ConnArray = new int[FullDesign.length][2];
        //design = EncodingUtils.getBinary(solution.getVariable());
        for (int index = 0; index < FullDesign.length; index++){
            boolean decision = FullDesign[index];
            if (decision){
                ConnArray[TrussCount] = FullConnectivityArray[index];
                TrussCount += 1;
            }
        }
        int[][] designConnArray;
        designConnArray = sliceFullConnectivityArray(ConnArray, TrussCount);
        return designConnArray;
    }

    private int getTrussCount(boolean[] CompleteDesign){
        int numVars = CompleteDesign.length;
        int TrussCount = 0;
        for (int index = 0; index < numVars; index++){
            boolean decision =  CompleteDesign[index];
            if (decision){
                TrussCount += 1;
            }
        }
        return TrussCount;
    }

    private boolean[] getBooleanDesignArray (Solution soln){
        int numVars = soln.getNumberOfVariables();
        boolean[] design = new boolean[numVars];
        for (int index = 0; index < numVars; index++){
            boolean decision =  EncodingUtils.getBoolean(soln.getVariable(index));
            design[index] = decision;
        }
        return design;
    }

    private boolean[] getCompleteBooleanDesignArray (boolean[] CurrentDesign){
        return new boolean[]{CurrentDesign[0], CurrentDesign[1], CurrentDesign[2], CurrentDesign[3], CurrentDesign[4],
                CurrentDesign[5], CurrentDesign[6], CurrentDesign[7], CurrentDesign[8], CurrentDesign[9],
                CurrentDesign[10], CurrentDesign[11], CurrentDesign[12], CurrentDesign[13], CurrentDesign[14],
                CurrentDesign[15], CurrentDesign[16], CurrentDesign[2], CurrentDesign[17], CurrentDesign[18],
                CurrentDesign[19], CurrentDesign[20], CurrentDesign[21], CurrentDesign[22], CurrentDesign[23],
                CurrentDesign[24], CurrentDesign[25], CurrentDesign[26], CurrentDesign[27], CurrentDesign[28],
                CurrentDesign[29], CurrentDesign[30], CurrentDesign[22], CurrentDesign[0], CurrentDesign[31],
                CurrentDesign[8]};
    }

    public boolean[] getBooleanDesignFromSolution (Solution solution) {
        //// Extract the design as a boolean array
        boolean[] BooleanDesign = getBooleanDesignArray(solution);
        //design = EncodingUtils.getBinary(solution.getVariable());

        //// Compute number of trusses in the design
        return getCompleteBooleanDesignArray(BooleanDesign);
    }

    public int[][] getConnectivityArrayFromSolution (Solution solution) {
        //// Extract the design as a boolean array
        boolean[] BooleanDesign = getBooleanDesignArray(solution);

        //// Obtain corresponding Connectivity Array
        return ConvertToFullConnectivityArray(BooleanDesign);
    }

    private int[][] sliceFullConnectivityArray (int[][] fullConnectivityArray, int trussCount) {
        int[][] connectivityArray = new int[trussCount][2];
        for (int i = 0; i < trussCount; i++) {
            for (int j = 0; j < 2; j++) {
                connectivityArray[i][j] = fullConnectivityArray[i][j];
            }
        }
        return connectivityArray;
    }

    public int getTrussCountFromSolution () {
        return NumberOfTrusses;
    }


}

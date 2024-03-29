package seakers.trussaos;

import seakers.aos.aos.AOS;
import seakers.aos.history.AOSHistoryIO;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Population;
import org.moeaframework.core.Solution;
import org.moeaframework.util.TypedProperties;

public class IntegerEvolutionarySearch implements Callable<Algorithm> {

    private final String savePath;
    private final String name;
    private final Algorithm alg;
    private final TypedProperties properties;
    private static double sidenum;
    private final double[] radii;
    private final double[] youngsModulii;
    private static int numHeurObjectives;
    private static int numHeurConstraints;

    public IntegerEvolutionarySearch(Algorithm alg, TypedProperties properties, String savePath, String name, double sideNodeNumber, double[] radii, double[] youngsModulii, int numHeuristicObjectives, int numHeuristicConstraints) {
        this.alg = alg;
        this.properties = properties;
        this.savePath = savePath;
        this.name = name;
        sidenum = sideNodeNumber;
        this.radii = radii;
        this.youngsModulii = youngsModulii;
        numHeurObjectives = numHeuristicObjectives;
        numHeurConstraints = numHeuristicConstraints;
    }

    @Override
    public Algorithm call() throws Exception {
        int populationSize = (int) properties.getDouble("populationSize", 100);
        int maxEvaluations = (int) properties.getDouble("maxEvaluations", 10000);

        // run the executor using the listener to collect results
        System.out.println("Starting " + alg.getClass().getSimpleName() + " on " + alg.getProblem().getName() + " with pop size: " + populationSize);
        alg.step();
        long startTime = System.currentTimeMillis();

        HashSet<Solution> allSolutions = new HashSet<>();
        Population initPop = ((AbstractEvolutionaryAlgorithm) alg).getPopulation();
        //int numberOfVariables = initPop.get(0).getNumberOfVariables();
        for (int i = 0; i < initPop.size(); i++) {
            initPop.get(i).setAttribute("NFE", 0);
            allSolutions.add( initPop.get(i));
        }

        while (!alg.isTerminated() && (alg.getNumberOfEvaluations() < maxEvaluations)) {
            if (alg.getNumberOfEvaluations() % 250 == 0) {
                System.out.println("NFE: " + alg.getNumberOfEvaluations());
                System.out.print("Popsize: " + ((AbstractEvolutionaryAlgorithm) alg).getPopulation().size());
                System.out.println("  Archivesize: " + ((AbstractEvolutionaryAlgorithm) alg).getArchive().size());
            }
            alg.step();
            Population pop = ((AbstractEvolutionaryAlgorithm) alg).getPopulation();
            for(int i=1; i<3; i++){
                Solution s = pop.get(pop.size() - i);
                s.setAttribute("NFE", alg.getNumberOfEvaluations());
                allSolutions.add(s);
            }
        }

        alg.terminate();
        long finishTime = System.currentTimeMillis();
        System.out.println("Done with optimization. Execution time: " + ((finishTime - startTime) / 1000) + "s");

        IntegerResultIO resultIOInteger = new IntegerResultIO(alg.getProblem().getName(),sidenum,radii,youngsModulii,numHeurObjectives,numHeurConstraints);;

        //String filename = savePath + File.separator + alg.getClass().getSimpleName() + "_" + name;
        //ResultIO.savePopulation(((AbstractEvolutionaryAlgorithm) alg).getPopulation(), filename);
        //ResultIO.savePopulation(new Population(allSolutions), filename + "_all");
        //ResultIO.saveObjectives(alg.getResult(), filename);

        String popCsvFilename;
        popCsvFilename = savePath + File.separator + alg.getClass().getSimpleName() + "_" + name + "_fullpop" + ".csv";
        resultIOInteger.savePopulationHistoryToCsvInteger(allSolutions, popCsvFilename);

        String csvFilename;
        csvFilename = savePath + File.separator + alg.getClass().getSimpleName() + "_" + name + ".csv";
        resultIOInteger.saveFinalResultToCsvInteger(((AbstractEvolutionaryAlgorithm) alg).getPopulation(), csvFilename);

        if (alg instanceof AOS) {
            AOS algAOS = (AOS) alg;
            if (properties.getBoolean("saveQuality", false)) {
                AOSHistoryIO.saveQualityHistory(algAOS.getQualityHistory(), new File(savePath + File.separator + name + "_qual" + ".csv"), ",");
            }
            if (properties.getBoolean("saveCredits", false)) {
                AOSHistoryIO.saveCreditHistory(algAOS.getCreditHistory(), new File(savePath + File.separator + name + "_credit" + ".csv"), ",");
            }
            if (properties.getBoolean("saveSelection", false)) {
                AOSHistoryIO.saveSelectionHistory(algAOS.getSelectionHistory(), new File(savePath + File.separator + name + "_hist" + ".csv"), ",");
            }
        }
        return alg;
    }
}

package seakers.trussaos;

import seakers.ahs.AHSMOEA;
import seakers.ahs.history.AHSHistoryIO;
import seakers.aos.aos.AOS;
import seakers.aos.history.AOSHistoryIO;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.atomic.AtomicInteger;
import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Population;
import org.moeaframework.core.Solution;
import org.moeaframework.util.TypedProperties;

public class EvolutionarySearch implements Callable<Algorithm> {

    private final String savePath;
    private final String name;
    private final Algorithm alg;
    private AtomicInteger nfe;
    private final TypedProperties properties;
    private static double sidenum;
    private HashSet<Solution> currentSolutions;
    private static int numHeurObjectives;
    private static int numHeurConstraints;

    public EvolutionarySearch(Algorithm alg, AtomicInteger nfe, TypedProperties properties, String savePath, String name, double sideNodeNumber, int numHeuristicObjectives, int numHeuristicConstraints) {
        this.alg = alg;
        this.nfe = nfe;
        this.properties = properties;
        this.savePath = savePath;
        this.name = name;
        sidenum = sideNodeNumber;
        this.currentSolutions = new HashSet<>();
        numHeurObjectives = numHeuristicObjectives;
        numHeurConstraints = numHeuristicConstraints;
    }

    public EvolutionarySearch(Algorithm alg, AtomicInteger nfe, TypedProperties properties, String savePath, String name, double sideNodeNumber, HashSet<Solution> currentSolutions, int numHeuristicObjectives, int numHeuristicConstraints) {
        this.alg = alg;
        this.nfe = nfe;
        this.properties = properties;
        this.savePath = savePath;
        this.name = name;
        sidenum = sideNodeNumber;
        this.currentSolutions = currentSolutions;
        numHeurObjectives = numHeuristicObjectives;
        numHeurConstraints = numHeuristicConstraints;
    }

    @Override
    public Algorithm call() throws IOException, ExecutionException, InterruptedException {

        int populationSize = (int) properties.getDouble("populationSize", 100);
        int maxEvaluations = (int) properties.getDouble("maxEvaluations", 10000);

        // run the executor using the listener to collect results
        System.out.println("Starting " + alg.getClass().getSimpleName() + " on " + alg.getProblem().getName() + " with pop size: " + populationSize);
        alg.step();
        long startTime = System.currentTimeMillis();

        //HashSet<Solution> allSolutions = new HashSet<>();
        Population initPop = ((AbstractEvolutionaryAlgorithm) alg).getPopulation();
        for (int i = 0; i < initPop.size(); i++) {
            initPop.get(i).setAttribute("NFE", 0);
            currentSolutions.add(initPop.get(i));
            //for (int j = 0; j < populationSize; j++) {
                //currentSolutions.add(initPop.get(j));
            //}
        }

        while (!alg.isTerminated() && (alg.getNumberOfEvaluations() < maxEvaluations)) {
            if (alg.getNumberOfEvaluations() % 250 == 0) {
                System.out.println("NFE: " + alg.getNumberOfEvaluations());
                System.out.print("Popsize: " + ((AbstractEvolutionaryAlgorithm) alg).getPopulation().size());
                System.out.println("  Archivesize: " + ((AbstractEvolutionaryAlgorithm) alg).getArchive().size());
            }
            nfe.set(alg.getNumberOfEvaluations());
            alg.step();
            Population pop = ((AbstractEvolutionaryAlgorithm) alg).getPopulation();
            for (int j = 0; j < populationSize; j++) {
                currentSolutions.add(pop.get(j));
            }
            for(int i=1; i<3; i++){
                Solution s = pop.get(pop.size() - i);
                s.setAttribute("NFE", alg.getNumberOfEvaluations());
                //currentSolutions.add(s);
            }
        }

        alg.terminate();
        long finishTime = System.currentTimeMillis();
        System.out.println("Done with optimization. Execution time: " + ((finishTime - startTime) / 1000) + "s");

        ResultIO resultIO = new ResultIO(alg.getProblem().getName(),sidenum,numHeurObjectives,numHeurConstraints);;

        //String filename = savePath + File.separator + alg.getClass().getSimpleName() + "_" + name;
        //ResultIO.savePopulation(((AbstractEvolutionaryAlgorithm) alg).getPopulation(), filename);
        //ResultIO.savePopulation(new Population(allSolutions), filename + "_all");
        //ResultIO.saveObjectives(alg.getResult(), filename);

        String popCsvFilename;
        popCsvFilename = savePath + File.separator + alg.getClass().getSimpleName() + "_" + name + "_fullpop" + ".csv";
        resultIO.savePopulationHistoryToCsv(currentSolutions, popCsvFilename);

        String csvFilename;
        csvFilename = savePath + File.separator + alg.getClass().getSimpleName() + "_" + name + ".csv";
        resultIO.saveFinalResultToCsv(((AbstractEvolutionaryAlgorithm) alg).getPopulation(), csvFilename);

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

        if (alg instanceof AHSMOEA) {
            if (properties.getBoolean("saveQuality", false)) {
                AHSHistoryIO.saveQualityHistory(((AHSMOEA) alg).getAdaptiveHeuristicSelection().getQualityHistory(), new File(savePath + File.separator + name + "_qual" + ".csv"), ",");
            }
            if (properties.getBoolean("saveCredits", false)) {
                AHSHistoryIO.saveCreditHistory(((AHSMOEA) alg).getAdaptiveHeuristicSelection().getCreditHistory(), new File(savePath + File.separator + name + "_credit" + ".csv"), ",");
            }
        }
        return alg;
    }
}

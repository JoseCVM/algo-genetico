package problems.qbf.solvers;

import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import metaheuristics.ga.AbstractGA;
import metaheuristics.ga.AbstractGA.Chromosome;
import metaheuristics.ga.AbstractGA.Population;
import problems.qbf.QBF;
import solutions.Solution;

/**
 * Metaheuristic GA (Genetic Algorithm) for
 * obtaining an optimal solution to a QBF (Quadractive Binary Function --
 * {@link #QuadracticBinaryFunction}). 
 * 
 * @author ccavellucci, fusberti
 */
public class GA_QBF extends AbstractGA<Integer, Integer> {

	/**
	 * Constructor for the GA_QBF class. The QBF objective function is passed as
	 * argument for the superclass constructor.
	 * 
	 * @param generations
	 *            Maximum number of generations.
	 * @param popSize
	 *            Size of the population.
	 * @param mutationRate
	 *            The mutation rate.
	 * @param filename
	 *            Name of the file for which the objective function parameters
	 *            should be read.
	 * @throws IOException
	 *             Necessary for I/O operations.
	 */
	public GA_QBF(Integer generations, Integer popSize, Double mutationRate, String filename) throws IOException {
		super(new QBF(filename), generations, popSize, mutationRate,0,false);
	}

	/**
	 * {@inheritDoc}
	 * 
	 * This createEmptySol instantiates an empty solution and it attributes a
	 * zero cost, since it is known that a QBF solution with all variables set
	 * to zero has also zero cost.
	 */
	@Override
	public Solution<Integer> createEmptySol() {
		Solution<Integer> sol = new Solution<Integer>();
		sol.cost = 0.0;
		return sol;
	}
	public Population fixChildren(Population children) {
		return children;
	}
	public Population latinHypercubeGet(int popMult) {
		int alleleNumber = 2;
		int popSize = popMult*alleleNumber;
		Set<List<Integer>> columns = new HashSet<List<Integer>>();
		while(columns.size() < this.ObjFunction.getDomainSize()) {
			List<Integer> range = IntStream.rangeClosed(0,popSize-1).boxed().collect(Collectors.toList());
			Collections.shuffle(range);
			columns.add(range);
		}
		int latinSquare[][] = new int[popSize][this.ObjFunction.getDomainSize()];
		int k = 0;
		for(List<Integer> l : columns) {
			int j = 0;
			for(int i : l) {
				latinSquare[k][j] = i%alleleNumber;
			}
		}
		Population p = new Population();
		for(int i = 0;i<popSize;i++) {
			Chromosome c = new Chromosome();
			for(int j = 0;j< this.ObjFunction.getDomainSize();j++) {
				c.add(latinSquare[i][j]);
			}
			p.add(c);
		}
		return p;	
	}
	/*
	 * (non-Javadoc)
	 * 
	 * @see metaheuristics.ga.AbstractGA#decode(metaheuristics.ga.AbstractGA.
	 * Chromosome)
	 */
	@Override
	protected Solution<Integer> decode(Chromosome chromosome) {

		Solution<Integer> solution = createEmptySol();
		for (int locus = 0; locus < chromosome.size(); locus++) {
			if (chromosome.get(locus) == 1) {
				solution.add(new Integer(locus));
			}
		}

		ObjFunction.evaluate(solution);
		return solution;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see metaheuristics.ga.AbstractGA#generateRandomChromosome()
	 */
	@Override
	protected Chromosome generateRandomChromosome() {

		Chromosome chromosome = new Chromosome();
		for (int i = 0; i < chromosomeSize; i++) {
			chromosome.add(rng.nextInt(2));
		}

		return chromosome;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see metaheuristics.ga.AbstractGA#fitness(metaheuristics.ga.AbstractGA.
	 * Chromosome)
	 */
	@Override
	protected Double fitness(Chromosome chromosome) {

		return decode(chromosome).cost;

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * metaheuristics.ga.AbstractGA#mutateGene(metaheuristics.ga.AbstractGA.
	 * Chromosome, java.lang.Integer)
	 */
	@Override
	protected void mutateGene(Chromosome chromosome, Integer locus) {

		chromosome.set(locus, 1 - chromosome.get(locus));

	}

	/**
	 * A main method used for testing the GA metaheuristic.
	 * 
	 */
	public static void main(String[] args) throws IOException {

		long startTime = System.currentTimeMillis();
		GA_QBF ga = new GA_QBF(1000, 100, 1.0 / 100.0, "instances/qbf400");
		Solution<Integer> bestSol = ga.solve();
		System.out.println("maxVal = " + bestSol);
		long endTime = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("Time = " + (double) totalTime / (double) 1000 + " seg");

	}

}

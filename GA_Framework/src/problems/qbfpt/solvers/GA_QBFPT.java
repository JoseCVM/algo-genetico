package problems.qbfpt.solvers;

import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import metaheuristics.ga.AbstractGA;
import metaheuristics.ga.AbstractGA.Chromosome;
import metaheuristics.ga.AbstractGA.Population;
import problems.qbf.QBF;
import problems.qbfpt.QBFPT;
import solutions.Solution;

/**
 * Metaheuristic GA (Genetic Algorithm) for
 * obtaining an optimal solution to a QBF (Quadractive Binary Function --
 * {@link #QuadracticBinaryFunction}). 
 * 
 * @author ccavellucci, fusberti
 */
public class GA_QBFPT extends AbstractGA<Integer, Integer> {

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
	public GA_QBFPT(Integer generations, Integer popSize, Double mutationRate, String filename,int latinMult,boolean steadyState) throws IOException {
		super(new QBFPT(filename), generations, popSize, mutationRate, latinMult,steadyState);
		this.prohibited_triples = this.mountProhibitedList();
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
				latinSquare[j][k] = i%alleleNumber;
				j++;
			}
			k++;
		}
		Population p = new Population();
		for(int i = 0;i<popSize;i++) {
			Chromosome c = new Chromosome();
			for(int j = 0;j< this.ObjFunction.getDomainSize();j++) {
				c.add(latinSquare[i][j]);
				//System.out.print(latinSquare[i][j]);
			}
			//System.out.println(" ");
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
	public Integer[][] prohibited_triples;

	public Integer[][] mountProhibitedList() {
		int size = this.ObjFunction.getDomainSize();
		Integer[][] triples = new Integer[size][3];
		for (int i = 0; i < size; i++) {
			triples[i][0] = i+1;

			if (lFunction(i, 131, 1031) != i+1) {
				triples[i][1] = lFunction(i, 131, 1031);
			} else {
				triples[i][1] = 1 + (lFunction(i, 131, 1031) % size);
			}

			Integer x = 1 + (lFunction(i, 193, 1093) % size);
			if (lFunction(i, 193, 1093) != i+1 && lFunction(i, 193, 1093) != triples[i][1]) {
				triples[i][2] = lFunction(i, 193, 1093);
			} else if (x != i+1 && x != triples[i][1]) {
				triples[i][2] = x;
			} else {
				triples[i][2] = 1 + ((lFunction(i, 193, 1093) + 1) % size);
			}
			Integer maxi = Math.max(triples[i][0], Math.max(triples[i][1], triples[i][2]));
			Integer mini = Math.min(triples[i][0], Math.min(triples[i][1], triples[i][2]));
			Integer middle = triples[i][0] + triples[i][1] + triples[i][2] - maxi - mini;
			triples[i][0] = mini-1;
			triples[i][1] = middle-1;
			triples[i][2] = maxi-1;
		}
		return triples;
	}
	
	public void printProhibitedList() {
		int size = this.ObjFunction.getDomainSize();
		for (int i = 0; i < size; i++) {
			System.out.println(prohibited_triples[i][0] + " " + prohibited_triples[i][1] + " " +  prohibited_triples[i][2]);
		}
	}

	private Integer lFunction(Integer u, Integer pi_1, Integer pi_2) {
		int size = this.ObjFunction.getDomainSize();
		return 1 + ((pi_1 * u + pi_2) % size);
	}
	
	/**
	 * A main method used for testing the GA metaheuristic.
	 * 
	 */
	public static void main(String[] args) throws IOException {
		long startTime = System.currentTimeMillis();
			GA_QBFPT ga = new GA_QBFPT(1000, 100, 1.0 / 100.0, "instances/qbf400",0,false);
		Solution<Integer> bestSol = ga.solve();
		System.out.println("maxVal = " + bestSol);
		long endTime = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("Time = " + (double) totalTime / (double) 1000 + " seg");

	}
	@Override
	public Population fixChildren(Population children) {
		Random r = new Random();
		int size = this.ObjFunction.getDomainSize();
		Population fixedChildren = new Population();
		for(Chromosome c : children) {
			for(int i = 0;i<size;i++) {
				if(c.get(this.prohibited_triples[i][0]) == 1 && c.get(this.prohibited_triples[i][1]) == 1 && c.get(this.prohibited_triples[i][2]) == 1) {
					int sai = r.nextInt(3);
					c.set(this.prohibited_triples[i][sai],0);
				}
			}			
			fixedChildren.add(c);
		}
		return fixedChildren;
	}

}

package problems.qbfpt.solvers;

import java.io.IOException;
import java.util.concurrent.TimeUnit;

import solutions.Solution;

public class GA_QBFPT_multithread extends GA_QBFPT implements Runnable{
    public GA_QBFPT_multithread(Integer generations, Integer popSize, Double mutationRate, String filename,int latinMult,boolean steadyState) throws IOException {
        super(generations,popSize,mutationRate,filename,latinMult,steadyState);
        this.done = false;
    }
    public Solution<Integer> finalSolution;
    public long totalTime;
    public Boolean done;
    public void run() {
        try {
            long startTime = System.currentTimeMillis();
            this.finalSolution = this.solve();
            this.totalTime = System.currentTimeMillis() - startTime;
            this.done = true;
        }catch(Exception e){
            e.printStackTrace();
        }
    }
    public void displaySolution() {
        System.out.println("maxVal = " + bestSol);
        System.out.println("Time = "+(double)totalTime/(double)1000+" seg");
    }
    public static void main(String args[]) throws IOException {
        String instances[] = {"instances/qbf020","instances/qbf040","instances/qbf060","instances/qbf080","instances/qbf100","instances/qbf200","instances/qbf400"};
        GA_QBFPT_multithread.verbose = false;
        int maxIter = 60 * 15 * 1000;
        int popSize = 100;
        int latinMult = 50;
        Double mutRate = 0.5/100.0, mutRateDiff = 0.15/100.0;
        for(int i = 6;i<7;i++) {
            System.out.println("Instancia: "+instances[i]);

            GA_QBFPT_multithread gaPadrao = new GA_QBFPT_multithread(maxIter,popSize, mutRate,instances[i], 0, false);
            GA_QBFPT_multithread gaDiffPopSize = new GA_QBFPT_multithread(maxIter,popSize*2, mutRate,instances[i], 0, false);
            GA_QBFPT_multithread gaDiffMutRate = new GA_QBFPT_multithread(maxIter,popSize,mutRateDiff, instances[i], 0, false);
            GA_QBFPT_multithread gaLatinCube = new GA_QBFPT_multithread(maxIter,popSize, mutRate,instances[i], latinMult, false);
            GA_QBFPT_multithread gaSteadyState = new GA_QBFPT_multithread(maxIter,popSize, mutRate,instances[i], 0, true);


            Thread gaPadraoThread = new Thread(gaPadrao);
            Thread gaDiffPopSizeThread = new Thread(gaDiffPopSize);
            Thread gaDiffMutRateThread = new Thread(gaDiffMutRate);
            Thread gaLatinCubeThread = new Thread(gaLatinCube);
            Thread gaSteadyStateThread = new Thread(gaSteadyState);


            gaPadraoThread.start();
            gaDiffPopSizeThread.start();
            gaDiffMutRateThread.start();
            gaLatinCubeThread.start();
            gaSteadyStateThread.start();


            try {
                TimeUnit.MILLISECONDS.sleep((long) (maxIter*1.05));
            } catch (InterruptedException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }

            System.out.println("GA padrao: ");
            gaPadrao.displaySolution();
            System.out.println("");
            System.out.println("GA popsize mudado: ");
            gaDiffPopSize.displaySolution();
            System.out.println("");
            System.out.println("GA mut rate mudado: ");
            gaDiffMutRate.displaySolution();
            System.out.println("");
            System.out.println("GA latin hypercube: ");
            gaLatinCube.displaySolution();
            System.out.println("");
            System.out.println("GA steady state: ");
            gaSteadyState.displaySolution();
            System.out.println("");
            System.out.println("");
            System.out.println("");
        }
    }
}
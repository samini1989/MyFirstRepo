using DNmodel.Common;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace DNmodel.GA
{
    public class GeneticAlgorithm
    {
        public int n_DGS = 2;///N of DG
        public int nBranch; // Nº of ramos
        public int nG = 200;//N of generations//200 
        public int nP = 20; // Nº of Individuos in population//20
        public individual[] Population;
        public individual best;
        public int n_switch = 10; // nº of  switices of ramos
        public bool defaultSwitchPosition = false; //if true, use default case position from xml
        public double mutation_factor;
        public double recombine_factor;
        public double p_crossover = 0.8;//  probability of crossover//.8
        public double p_mutation = 0.2;//  probability of mutation//.2
        public double cp = 1e6;
        public Network energetic;

        public void metodoAg(string _path)
        {
            // generate the initial  population
            Population = new individual[nP];
            best = new individual();

            best.fitness = 0.0;
            double fitness_pop = 0.0;
            for (int i = 0; i < nP; i++)
            {
                Thread.Sleep(1);
                Population[i] = new individual();
                int N_chor = energetic.buses.Length + n_DGS * Population[i].powerActive_k;//number of columns of matrix chromosome

                bool[][] n_cromossomo = new bool[2][];//row 0=position of switches, row 1=position of DGs   100+2*8 116     sw 100    
                n_cromossomo[0] = new bool[energetic.branches.Length];
                n_cromossomo[1] = new bool[N_chor];
                Population[i].chromosome = n_cromossomo;
                // ramos with switches
                Random aleatorio = new Random();




                for (int j = 0; j < n_switch; j++)
                {
                    int position = 0;
                    do
                    {
                        Thread.Sleep(1);
                        position = aleatorio.Next(1, energetic.branches.Length);

                    } while (Population[i].chromosome[0][position]);

                    Population[i].chromosome[0][position] = true;
                }

                Random aleatorio3 = new Random();
                for (int j = 0; j < n_DGS; j++)
                {
                    int position = 0;
                    do
                    {
                        Thread.Sleep(1);
                        position = aleatorio3.Next(1, energetic.buses.Length);

                    } while (Population[i].chromosome[1][position]);

                    Population[i].chromosome[1][position] = true;
                }


                //valor bin do Mse (valor  codificado)
                for (int j = energetic.buses.Length; j < N_chor; j++) // 
                {
                    Thread.Sleep(1);
                    if (aleatorio3.Next(0, 100) >= 50)
                    {
                        Population[i].chromosome[1][j] = true;
                    }

                }

                //cost of ens for each indivual configurations
                Population[i].decodification(energetic, cp, n_DGS, defaultSwitchPosition);
                fitness_pop += Population[i].fitness; // accumulates the fitness of the entire population
                if (Population[i].fitness > best.fitness) best = Population[i];



            }

            // generation counter here below
            Random aleatorio1 = new Random();  // WE STOP HERE, WE WILL START THE RECOMBINATION AND MUTATION PROCESS FOR THE Individuos
            individual[] pop_prole = new individual[nP];
            using (StreamWriter sw = new StreamWriter(_path))
            {
                for (int i = 0; i < nG; i++)
                {
                    for (int j = 0; j < nP; j++)
                    {
                        double arrow_rolet = aleatorio1.NextDouble() * fitness_pop;//?
                        // choose of parents
                        int indice_parent1 = 0;
                        double stop_rolet = 0.0;
                        while (stop_rolet < arrow_rolet && indice_parent1 < nP - 1)
                        {
                            stop_rolet += Population[indice_parent1].fitness;
                            indice_parent1++;
                        }
                        int indice_parent2;

                        // to guarantee the parent ratings are not the same
                        do
                        {
                            arrow_rolet = aleatorio1.NextDouble() * fitness_pop;
                            // choose two prantes
                            indice_parent2 = 0;
                            stop_rolet = 0;
                            while (stop_rolet < arrow_rolet && indice_parent2 < nP - 1)
                            {
                                stop_rolet += Population[indice_parent2].fitness;
                                indice_parent2++;
                            }

                        }
                        while (indice_parent1 == indice_parent2);

                        // calculate a probability de cross over
                        double prob = aleatorio1.NextDouble();
                        individual child = new individual();
                        if (prob < p_crossover)
                        {

                            // aplication o cross over of 1 ponto
                            child = recombine(Population[indice_parent1], Population[indice_parent2]);
                        }
                        else
                        {
                            child = Population[indice_parent1];
                            if (Population[indice_parent2].fitness > Population[indice_parent1].fitness)
                            {
                                child = Population[indice_parent2];
                            }

                        }
                        // calculate a probability de mutation
                        prob = aleatorio1.NextDouble();
                        if (prob < p_mutation)
                        {
                            muta(ref child);
                            child.decodification(energetic, cp, n_DGS, defaultSwitchPosition);
                        }
                        pop_prole[j] = child;
                        if (child.fitness > best.fitness)
                        {
                            best = new individual();
                            //best solution
                            best.fitness = child.fitness;

                            best.costofENS = child.costofENS;

                            best.chromosome = new bool[child.chromosome.Length][];

                            for (int k = 0; k < child.chromosome.Length; k++)
                            {
                                best.chromosome[k] = new bool[child.chromosome[k].Length];

                                for (int r = 0; r < child.chromosome[k].Length; r++)
                                {
                                    best.chromosome[k][r] = child.chromosome[k][r];
                                }
                            }

                        }
                    }
                    fitness_pop = 0.0;
                    for (int j = 0; j < Population.Length; j++)
                    {
                        Population[j] = pop_prole[j];
                        fitness_pop += Population[j].fitness;
                    }
                    sw.WriteLine("{0}: {1}", i.ToString(), best.fitness.ToString()); ///generate the best fitness file??
                    //DG[] list_dg = new DG[1];

                    //Switch[] list_sw = new Switch[1];
                    // best.getcoded(energetic, ref list_dg, ref list_sw);


                    bool printList = true;
                    int count = 1; //Counter of output data types
                    if (printList)
                    {
                        DG[] list_dg = new DG[1];

                        Switch[] list_sw = new Switch[1];
                        best.getcoded(energetic, ref list_dg, ref list_sw);


                        //Information DGs
                        sw.WriteLine("<<<" + count.ToString() + ". Information DGs>>>-----------------------------------------------------------------------------------------");
                        for (int j = 0; j < list_dg.Length; j++)
                            sw.WriteLine(" bFrombus:{0}   capacity:{1}  minpowerActive:{2}  maxpowerActive:{3}",
                                list_dg[j].bFrombus.ToString(),
                                 list_dg[j].capacity.ToString(),
                                 list_dg[j].minpowerActive.ToString(),
                                list_dg[j].maxpowerActive.ToString());
                        count++;

                        //Information SWitch
                        sw.WriteLine("<<<" + count.ToString() + ". Information SWitch>>>----------------------------------------------------------------------------------------");
                        for (int k = 0; k < list_sw.Length; k++)
                            sw.WriteLine("bFrom_switch:{0} bTo_switch:{1}  isclosed:{2}  ",
                                list_sw[k].bFrom_switch.ToString(),
                                list_sw[k].bTo_switch.ToString(),
                                 list_sw[k].isclosed.ToString()
                                );

                        count++;
                    }


                }


            }
        }

        public individual recombine(individual parent1, individual parent2)
        {
            // Here the recombination process takes place
            individual child1 = new individual();
            individual child2 = new individual();


            int conta_switch = 0;

            do
            {

                child1.chromosome = new bool[parent1.chromosome.Length][];
                child2.chromosome = new bool[parent2.chromosome.Length][];

                for (int k = 0; k < parent1.chromosome.Length; k++)
                {
                    child1.chromosome[k] = new bool[parent1.chromosome[k].Length];

                }


                for (int k = 0; k < parent2.chromosome.Length; k++)
                {
                    child2.chromosome[k] = new bool[parent2.chromosome[k].Length];

                }


                Random aleatorio = new Random();

                int ponto = aleatorio.Next(1, n_switch);
                int i = 0;

                int nsw11 = 0;
                int nsw12 = 0;
                while (i < energetic.branches.Length)
                {
                    if (parent1.chromosome[0][i])
                    {
                        if (nsw11 <= ponto) child1.chromosome[0][i] = parent1.chromosome[0][i];
                        else child2.chromosome[0][i] = parent1.chromosome[0][i];
                        nsw11++;
                    }
                    if (parent2.chromosome[0][i])
                    {
                        if (nsw12 <= ponto) child2.chromosome[0][i] = parent2.chromosome[0][i];
                        else child1.chromosome[0][i] = parent2.chromosome[0][i];
                        nsw12++;
                    }
                    i++;
                }
                // switches amount adjustments (CHILDREN CONSIDER THE SAME AMOUNT AS THEIR PARENTS)
                conta_switch = 0;
                for (int j = 0; j < energetic.branches.Length; j++)
                {
                    conta_switch = child1.chromosome[0][j] ? conta_switch + 1 : conta_switch;
                    conta_switch = child2.chromosome[0][j] ? conta_switch + 1 : conta_switch;
                }
            }
            while (conta_switch != (2 * n_switch));

            int conta_DG = 0;

            do
            {

                Random aleatorio = new Random();

                int ponto = aleatorio.Next(energetic.buses.Length, parent1.chromosome[1].Length);
                for (int i = 0; i <= ponto; i++)
                {
                    child1.chromosome[1][i] = parent1.chromosome[1][i];
                    child2.chromosome[1][i] = parent2.chromosome[1][i];
                }
                for (int i = ponto + 1; i < parent1.chromosome[1].Length; i++)
                {
                    child1.chromosome[1][i] = parent2.chromosome[1][i];
                    child2.chromosome[1][i] = parent1.chromosome[1][i];
                }
                // ajustes da  quantidade de FACTS  (FILHOS CONSIDERAREM A MESMA QUANTIDADE DOS PAIS)
                conta_DG = 0;
                for (int j = 0; j < energetic.buses.Length; j++)
                {
                    conta_DG = child1.chromosome[1][j] ? conta_DG + 1 : conta_DG;
                    conta_DG = child2.chromosome[1][j] ? conta_DG + 1 : conta_DG;
                }
            }
            while (conta_DG != (2 * n_DGS));




            child1.decodification(energetic, cp, n_DGS, defaultSwitchPosition);
            child2.decodification(energetic, cp, n_DGS, defaultSwitchPosition);

            if (child1.fitness > child2.fitness)
            {
                return child1;
            }
            else
            {
                return child2;
            }
        }

        public void muta(ref individual prole)
        {   //here the mutation process takes place


            Random aleatorio = new Random();
            int ponto = aleatorio.Next(1, energetic.branches.Length);

            if (ponto < energetic.branches.Length)
            {
                if (prole.chromosome[0][ponto])
                {
                    int ponto2 = 1;
                    do
                    {
                        ponto2 = aleatorio.Next(1, energetic.branches.Length);
                    }
                    while (prole.chromosome[0][ponto2]);
                    prole.chromosome[0][ponto2] = !prole.chromosome[0][ponto2];

                }
                else
                {
                    int selfacts1 = aleatorio.Next(1, n_switch + 1);
                    int i = 0;

                    while (selfacts1 > 0)
                    {
                        i++;

                        if (prole.chromosome[0][i])
                        {
                            selfacts1--;

                        }


                    }
                    prole.chromosome[0][i] = !prole.chromosome[0][i];

                }

            }

            prole.chromosome[0][ponto] = !prole.chromosome[0][ponto];



            Random aleatorio_dg = new Random();
            int ponto_dg = aleatorio_dg.Next(1, energetic.buses.Length);

            if (ponto_dg < energetic.buses.Length)
            {
                if (prole.chromosome[1][ponto_dg])
                {
                    int ponto2 = 1;
                    do
                    {
                        ponto2 = aleatorio_dg.Next(1, energetic.buses.Length);
                    }
                    while (prole.chromosome[1][ponto2]);
                    prole.chromosome[1][ponto2] = !prole.chromosome[1][ponto2];

                }
                else
                {
                    int selfacts = aleatorio_dg.Next(1, n_DGS + 1);
                    int j = 0;
                    while (selfacts > 0)
                    {
                        j++;

                        if (prole.chromosome[1][j])
                        {
                            selfacts--;
                        }
                    }
                    prole.chromosome[1][j] = !prole.chromosome[1][j];
                }
            }
            prole.chromosome[1][ponto_dg] = !prole.chromosome[1][ponto_dg];





        }
    }
}

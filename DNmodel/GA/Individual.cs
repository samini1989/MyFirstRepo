using DNmodel.Common;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace DNmodel.GA
{
    public class individual
    {
        // public double prob_rolet;
        //public Network energetic;

        public double fitness;
        public bool[][] chromosome;
        public int powerActive_k = 8;
        public double costofENSRestor;
        public double costofENSSwitching;
        public double costofENS;//cost of ENS of indivuial
                                // public Branch branches;
        public void getcoded(Network energetic, ref DG[] list_dgs, ref Switch[] list_switch)
        {
            //list of switches
            int cont_switch = 0;
            list_switch = new Switch[chromosome[0].Count(b => b == true)];
            for (int i = 0; i < chromosome[0].Length; i++)
            {
                if (chromosome[0][i])
                {
                    list_switch[cont_switch] = new Switch();
                    list_switch[cont_switch].bFrom_switch = energetic.branches[i].From;
                    list_switch[cont_switch].bTo_switch = energetic.branches[i].To;
                    list_switch[cont_switch].isclosed = true;
                    cont_switch++;

                }
            }


            int expoenent = 0;
            int position_inicial = energetic.buses.Length - 1;
            int conta_DG = 0;

            list_dgs = new DG[chromosome[1].TakeWhile((b, index) => index < energetic.buses.Length).Count(b => b == true)];

            for (int i = 0; i < energetic.buses.Length; i++)
            {

                if (chromosome[1][i])
                {
                    list_dgs[conta_DG] = new DG();
                    list_dgs[conta_DG].bFrombus = i;
                    conta_DG++;
                }
            }

            //list of DGs
            for (int i = 0; i < list_dgs.Length; i++)
            {

                list_dgs[i].capacity = 0.0;
                for (int j = position_inicial + powerActive_k * i; j < position_inicial + powerActive_k * (i + 1); j++)
                {
                    list_dgs[i].capacity += this.chromosome[1][j] ? Math.Pow(2, expoenent) : 0;
                    expoenent++;
                }
                // expoenent capacity
                list_dgs[i].capacity = list_dgs[i].minpowerActive + list_dgs[i].capacity * list_dgs[i].maxpowerActive / (Math.Pow(2, powerActive_k) - 1);
                expoenent = 0;
            }
        }

        public void decodification(Network energetic, double c1, int dg3, bool defaultSwitch)
        {


            if (defaultSwitch == false)
            {
                for (int i = 0; i < energetic.branches.Length; i++)
                {
                    energetic.branches[i].switch_ = chromosome[0][i];
                }

            }



            void obtainSections()//obtainSection 
            {
                //defind the sectionalizing branch, lenght section, demand load           

                void Sectionlizing_branch(Bus advanced, Section newsection)
                {
                    //Associate a tie switch
                    for (int j = 0; j < energetic.tieSwitches.Length; j++)
                    {
                        if (advanced.Number == energetic.tieSwitches[j].bTo_tieswitch)
                        {
                            newsection.tie_switch_section = energetic.tieSwitches[j];
                            break;
                        }
                    }

                    for (int i = 0; i < advanced.DownBranches.Count; i++)
                    {
                        if (advanced.DownBranches[i].switch_ == true)
                        {
                            newsection.SectionalizingBranches.Add(advanced.DownBranches[i]);
                        }
                        else
                        {
                            newsection.length_section += advanced.DownBranches[i].length_bran;
                            Sectionlizing_branch(advanced.DownBranches[i].bTo, newsection);
                        }

                    }



                }

                //define the number of sections
                energetic.root.DownBranches[0].switch_ = true;
                for (int i = 0; i < energetic.branches.Length; i++)
                {
                    if (energetic.branches[i].switch_)
                    {
                        Section newSection = new Section();
                        newSection.rootBranch = energetic.branches[i];
                        newSection.residentiallLoad = new double[3];
                        newSection.industrialLoad = new double[3];
                        newSection.commercialLoad = new double[3];
                        newSection.totalLoad = new double[3];
                        newSection.length_section = 0;
                        energetic.sections.Add(newSection);
                        newSection.statusOfFault = false;
                        newSection.statusOfSection = false;
                        //newSection.failureRate = 0.025;
                        newSection.failureRate = 0.0;
                        //newSection.CENS = new double[3];
                    }
                }
                energetic.totalLengthOffeeder = 0;
                //Make the  Sectionalizing Branches Association
                for (int j = 0; j < energetic.sections.Count; j++)
                {
                    Sectionlizing_branch(energetic.sections[j].rootBranch.bTo, energetic.sections[j]);
                    // total_load(sections[j].rootBranch.bTo, sections[j]);
                    energetic.totalLengthOffeeder += energetic.sections[j].length_section;


                }

                //Make the upstream section association
                for (int i = 0; i < energetic.sections.Count; i++)
                {
                    //failure rate of section j
                    energetic.sections[i].failureRate = (0.025) * energetic.sections[i].length_section * energetic.totalLengthOffeeder;


                    for (int j = 0; j < energetic.sections.Count; j++)
                    {
                        for (int k = 0; k < energetic.sections[j].SectionalizingBranches.Count; k++)
                        {
                            if (energetic.sections[j].SectionalizingBranches[k] == energetic.sections[i].rootBranch)
                            {
                                energetic.sections[i].UpSection = energetic.sections[j];
                                break;
                            }
                        }
                        if (energetic.sections[i].UpSection != null)
                        {
                            break;
                        }
                    }
                }

            }//obtain sections

            //calculate cost of 
            obtainSections();

            void statusofcustomer()
            {
                //Define the status of sections witout neighbour feeder
                for (int i = 0; i < energetic.sections.Count; i++)
                {
                    bool changedSS = false; //flag indicating the status of any section was changed
                    for (int j = 0; j < energetic.sections.Count; j++)
                    {
                        if (energetic.sections[j].statusOfFault == true)
                        {
                            energetic.sections[j].statusOfSection = true;
                        }
                        if (energetic.sections[j].UpSection != null && energetic.sections[j].UpSection.statusOfSection == true && energetic.sections[j].statusOfSection == false)
                        {
                            energetic.sections[j].statusOfSection = true;

                            changedSS = energetic.sections[j].statusOfSection;
                        }
                    }
                    if (!changedSS) break;
                }

                //Check the neighbour feeder to finilize the definition of section state
                for (int j = 0; j < energetic.sections.Count; j++)
                {
                    double pdg = energetic.sections[j].DGs_section == null ? 0.0 : energetic.sections[j].DGs_section.capacity;
                    if (!energetic.sections[j].statusOfFault && energetic.useTS && energetic.sections[j].tie_switch_section != null && energetic.sections[j].tie_switch_section.capacity >= (energetic.sections[j].totalLoad.Sum() - pdg))
                    {
                        energetic.sections[j].statusOfSection = false;

                    }
                }

                //Check the DGs to finilize the definition of section state
                for (int j = 0; j < energetic.sections.Count; j++)
                {
                    if (!energetic.sections[j].statusOfFault && energetic.sections[j].DGs_section != null && energetic.sections[j].DGs_section.capacity >= energetic.sections[j].totalLoad.Sum())
                    {
                        energetic.sections[j].statusOfSection = false;

                    }
                }

                //Define the status of custumer
                void getSC(Bus advanced, Section newsection)
                {
                    advanced.statusofthejthcustomer = newsection.statusOfSection;

                    for (int i = 0; i < advanced.DownBranches.Count; i++)
                    {
                        if (advanced.DownBranches[i].switch_ == false)
                        {
                            getSC(advanced.DownBranches[i].bTo, newsection);
                        }

                    }
                } //Recurssive function for each section

                //Call the recurssive function to every section in the network
                //for (int j = 0; j < energetic.sections.Count; j++)
                // {
                //    getSC(energetic.sections[j].rootBranch.bTo, energetic.sections[j]);
                //  }
            }


            //calculo do fiteness
            DG[] DGs = new DG[dg3];

            int conta_DG = 0;

            for (int i = 0; i < energetic.buses.Length; i++)
            {
                energetic.buses[i].DGs = null;

                if (chromosome[1][i])
                {
                    energetic.buses[i].DGs = DGs[conta_DG];
                    DGs[conta_DG] = new DG();
                    DGs[conta_DG].bFrombus = i;
                    conta_DG++;

                }
            }


            // for (int k = 0; k < conta_DG; k++)
            // {
            //     DGs[k].maxpowerActive = 0.0;

            // determination poweractive
            //    int expoenent1 = 0;

            //for (int i = 0; i < DGs[k].maxpowerActive; i++)
            //{
            //  for (int j = position_inicial + powerActive_k * i; j < position_inicial + powerActive_k * (i + 1); j++)
            // {
            //       DGs[k].maxpowerActive += this.chromosome[1][j] ? Math.Pow(2, expoenent1) : 0;
            //          expoenent1++;
            //    }
            // expoenent poweractive
            //     DGs[k].maxpowerActive *= 0.5 / (Math.Pow(2, powerActive_k) - 1);
            //     expoenent1 = 0;
            //  }
            //determination capacity
            // position_inicial = energetic.buses.Length + ((k + 1) * mse_k + k * teta_k) * FACTS[k].teta_se.Length;
            int expoenent = 0;
            int position_inicial = energetic.buses.Length - 1;

            for (int i = 0; i < DGs.Length; i++)
            {
                DGs[i].capacity = 0.0;
                for (int j = position_inicial + powerActive_k * i; j < position_inicial + powerActive_k * (i + 1); j++)
                {
                    DGs[i].capacity += this.chromosome[1][j] ? Math.Pow(2, expoenent) : 0;
                    expoenent++;
                }
                // expoenent capacity
                DGs[i].capacity = DGs[i].minpowerActive + DGs[i].capacity * DGs[i].maxpowerActive / (Math.Pow(2, powerActive_k) - 1);
                expoenent = 0;
            }




            for (int j = 0; j < energetic.buses.Length; j++)
            {
                energetic.buses[j].DGs = null;
                for (int k = 0; k < conta_DG; k++)
                {
                    if (energetic.buses[j].Number == DGs[k].bFrombus)
                    {
                        energetic.buses[j].DGs = DGs[k];
                        break;
                    }
                }

            }


            // energetic.calculatePowerFlow(6, 7);



            energetic.list_CENSnetwork_Switching.Clear();
            energetic.list_CENSnetwork_restor.Clear();


            Random rnd = new Random(); //Randon variable 504
                                       // int max_omega1 = 4;//set of time series comprised in the interruption  time interval;
                                       //int x = rnd.Next(0, energetic.Curves[0].dataCurve.Count - max_omega1 - 3);
            int x = 6; //for x at min or max load x=201(max) and x=6(min)
            statusofcustomer();
            int omega_max = 4;

            int omega1 = 0;


            while (omega1 < omega_max)
            {
                //Solve the power flow and get  the percentage of energy consumption
                energetic.calculatePowerFlow(x, x + 1);

                //Calculate the total for each feeder section
                void total_load(Bus advanced, Section newsection)
                {
                    double PDG = 0;
                    if (advanced.DGs != null)
                    {
                        newsection.DGs_section = advanced.DGs;
                        PDG = advanced.DGs.capacity / 3;

                    }
                    for (int j = 0; j < 3; j++)
                    {
                        if (advanced.type == "residential")
                        {
                            newsection.residentiallLoad[j] += (Complex.Conjugate(advanced.I[j]) * advanced.V[j]).Real + PDG;
                        }
                        if (advanced.type == "commercial")
                        {
                            newsection.commercialLoad[j] += (Complex.Conjugate(advanced.I[j]) * advanced.V[j]).Real + PDG;
                        }
                        if (advanced.type == "industrial")
                        {
                            newsection.industrialLoad[j] += (Complex.Conjugate(advanced.I[j]) * advanced.V[j]).Real + PDG;
                        }
                        newsection.totalLoad[j] += (Complex.Conjugate(advanced.I[j]) * advanced.V[j]).Real;



                        for (int k = 0; k < 3; k++)
                        {
                            newsection.totalLoad[j] += Math.Pow(Complex.Abs(advanced.UpBranch.J[j]), 2) * advanced.UpBranch.Z[k][j].Real;
                        }
                    }
                    for (int i = 0; i < advanced.DownBranches.Count; i++)
                    {
                        if (advanced.DownBranches[i].switch_ == false)
                        {
                            total_load(advanced.DownBranches[i].bTo, newsection);
                        }

                    }

                }

                energetic.costENSofNetwork_restor = new double[3];
                energetic.costENSofNetwork_Switching = new double[3];

                //Calculate the network C_ENS at time x
                for (int i = 0; i < energetic.sections.Count; i++)
                {
                    total_load(energetic.sections[i].rootBranch.bTo, energetic.sections[i]);

                    energetic.sections[i].statusOfFault = true;
                    if (energetic.USEFLISER == false)
                    {
                        statusofcustomer(); //Default condition where downstream feeder section is interrupted
                    }
                    else
                    {
                        energetic.sections[i].statusOfSection = true; //Just faulted section is interrupted
                    }


                    //Put the code to calculate the CostRestoration of ENS and the CostSwitching of ENS

                    //cost of energy not supplied for jth section in the ith time step($);furmula(6)
                    double[] costENSofResidential_restor = new double[3];
                    double[] costENSOfCommercial_restor = new double[3];
                    double[] costENSOfIndustrial_restor = new double[3];
                    energetic.costENS_restor = new double[3];

                    double[] costENSofResidential_Switching = new double[3];
                    double[] costENSOfCommercial_Switching = new double[3];
                    double[] costENSOfIndustrial_Switching = new double[3];
                    energetic.costENS_Switching = new double[3];

                    if (omega1 != 0 || energetic.USEFLISER == true)
                    {
                        energetic.electricityRatesResidential_restor = 7.0;//residentialS  restorCost
                        energetic.electricityRatesCommercialS_restor = 200.0;//commercial  restorCost
                        energetic.electricityRatesIndustrial_restor = 30.0;//industrial    restorCost

                        //Put the code to calculate the CostRestoration of ENS 

                        for (int j = 0; j < energetic.sections.Count; j++)
                        {
                            if (energetic.sections[j].statusOfSection)
                            {
                                for (int k = 0; k < 3; k++)
                                {
                                    costENSofResidential_restor[k] += energetic.electricityRatesResidential_restor * energetic.sections[j].residentiallLoad[k] * 1e-3;
                                    costENSOfCommercial_restor[k] += energetic.electricityRatesCommercialS_restor * energetic.sections[j].commercialLoad[k] * 1e-3;
                                    costENSOfIndustrial_restor[k] += energetic.electricityRatesIndustrial_restor * energetic.sections[j].industrialLoad[k] * 1e-3;

                                    energetic.costENS_restor[k] = costENSofResidential_restor[k] + costENSOfCommercial_restor[k] + costENSOfIndustrial_restor[k];
                                    energetic.costENSofNetwork_restor[k] += energetic.costENS_restor[k] * energetic.sections[i].failureRate * 8.76e-6;
                                }
                            }
                        }

                    }
                    else
                    {
                        energetic.electricityRatesResidential_Switching = 0.55;//residential  SwitchingCost
                        energetic.electricityRatesCommercialS_Switching = 4.0;//commercial   SwitchingCost
                        energetic.electricityRatesIndustrial_Switching = 5.0;//industrial     SwitchingCost
                        for (int j = 0; j < energetic.sections.Count; j++)
                            energetic.sections[j].statusOfSection = true;


                        //Put the code to calculate  and the CostSwitching of ENS


                        for (int j = 0; j < energetic.sections.Count; j++)
                        {
                            if (energetic.sections[j].statusOfSection)
                            {
                                for (int k = 0; k < 3; k++)
                                {
                                    costENSofResidential_Switching[k] += energetic.electricityRatesResidential_Switching * energetic.sections[j].residentiallLoad[k] * 1e-3;
                                    costENSOfCommercial_Switching[k] += energetic.electricityRatesCommercialS_Switching * energetic.sections[j].commercialLoad[k] * 1e-3;
                                    costENSOfIndustrial_Switching[k] += energetic.electricityRatesIndustrial_Switching * energetic.sections[j].industrialLoad[k] * 1e-3;

                                    energetic.costENS_Switching[k] = costENSofResidential_Switching[k] + costENSOfCommercial_Switching[k] + costENSOfIndustrial_Switching[k];
                                    energetic.costENSofNetwork_Switching[k] += energetic.costENS_Switching[k] * energetic.sections[i].failureRate * 8.76e-6;
                                }
                            }
                        }

                    }

                    //Clear the sections atributes associated to fault
                    energetic.sections[i].statusOfFault = false;
                    for (int j = 0; j < energetic.sections.Count; j++) energetic.sections[j].statusOfSection = false;
                }

                energetic.list_CENSnetwork_Switching.Add(energetic.costENSofNetwork_Switching);
                energetic.list_CENSnetwork_restor.Add(energetic.costENSofNetwork_restor);


                //Clear the load of section

                for (int i = 0; i < energetic.sections.Count; i++)
                {
                    energetic.sections[i].residentiallLoad = new double[3];
                    energetic.sections[i].commercialLoad = new double[3];
                    energetic.sections[i].industrialLoad = new double[3];
                    energetic.sections[i].totalLoad = new double[3];

                }


                omega1++;
                x += omega1;
            }

            //calculate the total cost_restor of ENS of network  for 3 phases


            for (int i = 0; i < energetic.list_CENSnetwork_restor.Count; i++)
            {
                for (int k = 0; k < 3; k++)
                {

                    costofENSRestor += energetic.list_CENSnetwork_restor[i][k];
                }
            }

            //calculate the total cost_switch of ENS of network  for 3 phases

            for (int i = 0; i < energetic.list_CENSnetwork_Switching.Count; i++)
            {
                for (int k = 0; k < 3; k++)
                {

                    costofENSSwitching += energetic.list_CENSnetwork_Switching[i][k];
                }
            }


            costofENS = costofENSSwitching + costofENSRestor;

            //calculate of fiteness
            this.fitness = Double.IsNaN(costofENS) ? 0 : c1 / (costofENS);
            energetic.sections.Clear();

        }
    }

}

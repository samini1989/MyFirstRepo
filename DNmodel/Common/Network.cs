using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Xml;

namespace DNmodel.Common
{

    public class Network
    {
        //Atributes
        public Bus[] buses; //Network buses
        public Branch[] branches; //Network branches
        public TieSwitch[] tieSwitches; //Network tie switches
        public DG[] generators;
        public List<Section> sections = new List<Section>(); //Network sections
        public Bus root; //Root bus of the distribution network
        public Complex[] Vse; //Voltage reference Substation
        public List<double> absolutePower = new List<double>(); //Abs of the difference between calculated and specified power injection

        public List<Curve> Curves = new List<Curve>();

      
        public double lambda_0 = 0.0;
        public double lambda_1 = 0.5;
        public double lambda_2 = 0.5;
        public Complex[] SZIP = new Complex[3];//ZIP load model
        public List<double[]> R_TH = new List<double[]>();
        public List<double[]> F_TH = new List<double[]>();
        //ZIP load model
        public List<Complex[]> percentageEnergy = new List<Complex[]>();//The function of the percentage of energy consumption

        public List<double[]> normalizedOfPowerofConsumed = new List<double[]>();//public double[] normalizedOfPowerofConsumed = new double[3];//the normalized magnitude of power Consumed
        public List<double[]> powerFactorofPowerofConsumed = new List<double[]>();//the normalized magnitude of of power Consumed

        public double[] alphap = new double[3];//the random real variables, from 0 to 1, and
        public double[] alphapf = new double[3];//the random real variables, from 0 to 1, and
        public List<double[]> projectionp = new List<double[]>();
        public List<double[]> projectionpf = new List<double[]>();//the projections of the points NR, and NF over the abscissa indicate the values of rpi (ts ) and f pi (ts ),
        public double w = 0.9;//Defaull = 0.7   the variable ωlimits the fluctuation range between the maximum and minimum projection.

        public int sumLenghtsOfbranches;//summation of total lengths of branches.
        public List<double> probabilityOfselectionofBranch_fault = new List<double>();//a probabilidade de seleção do ramo onde a falt
        public List<double> frequenyMin = new List<double>();//
        public List<double> frequenyMax = new List<double>();//frequenyMin and frequenyMax limit the range of the accumulated frequency of the chosen path j.
        public List<double> selectOfrandomVariable = new List<double>();//the discrete selection variable, Sel(γ), determines a path as a function of the random variable γ.
        public double randomVariable;//a function of the random variable γ
        public double[] probabilityOfEachTypeOfFault = new double[4];
        public List<double> probabilityOfselectofEachhourwithFault = new List<double>();        //probabilityOfselectofEachhourwithFault
        public double deltaTime_switch;
        public double deltaTime_fault;
        public Complex eneryNotSupplied_BS;
        public Complex eneryNotSupplied_AS;
        public Complex eneryNotSupplied_total;
        public int numberOfSections;
        public double lossepoweractive;

        public double[] listOfCENSNetOFTime;


        public double costofenergynotsupplied;//cost of energy not supplied
        public double[] LoadPercentegeDemandrese;//the load percentage demand hour-by-hour
        public List<double[]> LoadPercentegeDemandcom = new List<double[]>();//the load percentage demand hour-by-hour
        public double[] LoadPercentegeDemandind;//the load percentage demand hour-by-hour

        public bool USEFLISER = false; // if true defult setuation ideal condition
        public bool useTS = true;//use trasfration tie sw

        public double probabilityOffaultNetwork;
        public double electricityRatesResidential;//residential  
        public double electricityRatesCommercialS;//commercial  
        public double electricityRatesIndustrial;//Industrial 

        public double electricityRatesResidential_restor;//residentialSwitchingCost
        public double electricityRatesCommercialS_restor;//commercialSwitchingCost
        public double electricityRatesIndustrial_restor;//industrialSwitchingCost
        public List<double[]> list_CENSnetwork_restor = new List<double[]>();
        public double[] costENS_restor; //Interruption cost of energy not supplied
        public double[] costENSofNetwork_restor; //Interruption cost of energy not supplied of network


        public double electricityRatesResidential_Switching;//residentialSwitchingCost
        public double electricityRatesCommercialS_Switching;//commercialSwitchingCost
        public double electricityRatesIndustrial_Switching;//industrialSwitchingCost
        public List<double[]> list_CENSnetwork_Switching = new List<double[]>();
        public double[] costENS_Switching; //Interruption cost of energy not supplied
        public double[] costENSofNetwork_Switching; //Interruption cost of energy not supplied of network


        public List<double[]> list_CENSnetwork = new List<double[]>();
        double[] costENS; //Interruption cost of energy not supplied
        double[] costENSofNetwork; //Interruption cost of energy not supplied of network
        double[] costIC; //Interruption cost of custumer
        double[] costICofNetwork; //Interruption cost of custumer
        public List<double[]> list_ICnetwork = new List<double[]>();
        public double totalLengthOffeeder;
        public double[] installedPower;
        public bool FultHappenNetwork = false;
        
        //  public int sumLenghtsOfbranches;//summation of total lengths of branches.
        public List<double> probabilityOfselectionofsection_fault = new List<double>();//a probabilidade de seleção do ramo onde a falt
        public List<double> frequenyMin_section = new List<double>();//
        public List<double> frequenyMax_section = new List<double>();//frequenyMin and frequenyMax limit the range of the accumulated frequency of the chosen path j.
        public List<double> selectOfrandomVariable_section = new List<double>();//the discrete selection variable, Sel(γ), determines a path as a function of the random variable γ.
        public double randomVariable_section;//a function of the random variable γ
        public double[] probabilityOfEachTypeOfFault_section = new double[4];
        public List<double> probabilityOfselectofEachhourwithFault_section = new List<double>();        //probabilityOfselectofEachhourwithFault

        //Functions
        static double NextDouble(Random rand, double minValue, double maxValue)
        {
            return rand.NextDouble() * (maxValue - minValue) + minValue;
        }

        public double variance(double[] sources)
        {
            double meanSquares = 0.0, mean = 0.0;
            for (int i = 0; i < sources.Length; i++)
            {
                meanSquares += Math.Pow(sources[i], 2);
                mean += sources[i];
            }
            double d = meanSquares / (sources.Length - 1);
            double d1 = Math.Pow(mean, 2) / (sources.Length * (sources.Length - 1));
            return d - d1;

        }

        public void readNetwork(XmlDocument inputDN)
        {
            //Define Substation Voltage
            Vse = new Complex[3];
            Vse[0] = new Complex(13800 / Math.Sqrt(3), 0.0);
            Vse[1] = new Complex(Vse[0].Magnitude * Math.Cos(4 * Math.PI / 3), Vse[0].Magnitude * Math.Sin(4 * Math.PI / 3));
            Vse[2] = new Complex(Vse[0].Magnitude * Math.Cos(2 * Math.PI / 3), Vse[0].Magnitude * Math.Sin(2 * Math.PI / 3));

            XmlNode branchXML = inputDN.DocumentElement.SelectSingleNode("Branch");
            XmlNode busXML = inputDN.DocumentElement.SelectSingleNode("Bus");
            XmlNode LoadCurveXML = inputDN.DocumentElement.SelectSingleNode("LoadCurve");
            XmlNode TieSwitchXML = inputDN.DocumentElement.SelectSingleNode("TieSwitch");
            XmlNode GeneratorXML = inputDN.DocumentElement.SelectSingleNode("Generator");
            //Reading the input data TieSwitch
            tieSwitches = new TieSwitch[TieSwitchXML.ChildNodes.Count];

            for (int i = 0; i < tieSwitches.Length; i++)
            {
                tieSwitches[i] = new TieSwitch();
                tieSwitches[i].capacity = Convert.ToDouble(TieSwitchXML.ChildNodes[i].Attributes["Capacity"].Value) * 1e3; //Default = 1e3
                tieSwitches[i].bFrom_tieswitch = Convert.ToInt32(TieSwitchXML.ChildNodes[i].Attributes["from"].Value);
                tieSwitches[i].bTo_tieswitch = Convert.ToInt32(TieSwitchXML.ChildNodes[i].Attributes["to"].Value);
            }


            //Reading the input data Genarator
            generators = new DG[GeneratorXML.ChildNodes.Count];

            for (int i = 0; i < generators.Length; i++)
            {
                generators[i] = new DG();
                generators[i].capacity = Convert.ToDouble(GeneratorXML.ChildNodes[i].Attributes["p"].Value);
                generators[i].bFrombus = Convert.ToInt32(GeneratorXML.ChildNodes[i].Attributes["bus"].Value);
            }

            //Reading the input data curves

            Curve LC = new Curve();
            for (int i = 0; i < LoadCurveXML.ChildNodes.Count; i++)
            {
                DataCurve data = new DataCurve();
                data.year = Convert.ToInt32(LoadCurveXML.ChildNodes[i].Attributes["year"].Value);
                data.month = Convert.ToInt32(LoadCurveXML.ChildNodes[i].Attributes["month"].Value);
                data.day = Convert.ToInt32(LoadCurveXML.ChildNodes[i].Attributes["day"].Value);
                data.hour = Convert.ToInt32(LoadCurveXML.ChildNodes[i].Attributes["hour"].Value);
                data.p = Convert.ToDouble(LoadCurveXML.ChildNodes[i].Attributes["p"].Value);
                data.pf = Convert.ToDouble(LoadCurveXML.ChildNodes[i].Attributes["pf"].Value);
                LC.dataCurve.Add(data);
            }
            Curves.Add(LC);
            //List<double[]> R_TH = new List<double[]>();
            //List<double[]> F_TH = new List<double[]>();
            for (int i = 1; i < 25; i++)
            {
                List<double> store_R = new List<double>();
                List<double> store_F = new List<double>();
                for (int j = 0; j < Curves[0].dataCurve.Count; j++)
                {
                    if (Curves[0].dataCurve[j].hour == i)
                    {
                        store_R.Add(Curves[0].dataCurve[j].p);
                        store_F.Add(Curves[0].dataCurve[j].pf);

                    }
                }
                R_TH.Add(store_R.ToArray());
                F_TH.Add(store_F.ToArray());
                store_R.Clear();
                store_F.Clear();
            }

            //the variance of hourly Normalized 
            Curves[0].varianceOfHourlyNormalized = new double[R_TH.Count];
            Curves[0].deviationOfHourlyNormalized = new double[R_TH.Count];
            for (int i = 0; i < R_TH.Count; i++)
            {
                Curves[0].varianceOfHourlyNormalized[i] = variance(R_TH[i]);
                Curves[0].deviationOfHourlyNormalized[i] = Math.Sqrt(Curves[0].varianceOfHourlyNormalized[i]);
            }
            //the variance of hourly power factor 
            Curves[0].varianceOfHourlyOfPowerFactor = new double[F_TH.Count];
            Curves[0].deviationOfHourlyOfPowerFactor = new double[F_TH.Count];

            for (int i = 0; i < F_TH.Count; i++)
            {
                Curves[0].varianceOfHourlyOfPowerFactor[i] = variance(F_TH[i]);

                Curves[0].deviationOfHourlyOfPowerFactor[i] = Math.Sqrt(variance(F_TH[i]));
            }

            //the projections of the points NR, and NF over the abscissan
            Random rand = new Random();
            for (int i = 0; i < R_TH.Count; i++)
            {
                double[] projectP = new double[3];

                for (int k = 0; k < 3; k++)
                {
                    alphap[k] = NextDouble(rand, 0, 1);

                    projectP[k] = (w + alphap[k] * (1 - w)) / (Curves[0].deviationOfHourlyNormalized[i]
                     * Math.Sqrt(2 * Math.PI));
                }
                projectionp.Add(projectP);
            }
            for (int i = 0; i < F_TH.Count; i++)
            {
                double[] projectPf = new double[3];
                for (int k = 0; k < 3; k++)
                {
                    alphapf[k] = NextDouble(rand, 0, 1);
                    projectPf[k] = (w + alphapf[k] * (1 - w)) / (Curves[0].deviationOfHourlyOfPowerFactor[i]
                        * Math.Sqrt(2 * Math.PI));
                }
                projectionpf.Add(projectPf);
            }
            //the probability density function strongly depends on the calculated variance
            for (int i = 0; i < R_TH.Count; i++)
            {
                double[] nPc = new double[3];
                // double[] fPc = new double[3];

                for (int k = 0; k < 3; k++)
                {
                    // double  beta = NextDouble(rand, 0, 1);
                    nPc[k] = Math.Log(1 / (projectionp[i][k] * Curves[0].deviationOfHourlyNormalized[i] *
                        Math.Sqrt(2 * Math.PI)));

                    nPc[k] = Math.Sqrt(2 * Curves[0].varianceOfHourlyNormalized[i] * nPc[k]);

                    // nPc[k] =Math.Sign(beta-0.5)* nPc[k] + Curves[0].dataCurve[i].p;
                }
                normalizedOfPowerofConsumed.Add(nPc);
                //powerFactorofPowerofConsumed.Add(fPc);
            }
            for (int i = 0; i < F_TH.Count; i++)
            {
                //  double[] nPc = new double[3];
                double[] fPc = new double[3];

                for (int k = 0; k < 3; k++)
                {
                    //double beta = NextDouble(rand, 0, 1);
                    fPc[k] = Math.Log(1 / (projectionp[i][k] * Curves[0].deviationOfHourlyOfPowerFactor[i] *
                        Math.Sqrt(2 * Math.PI)));

                    fPc[k] = Math.Sqrt(2 * Curves[0].varianceOfHourlyOfPowerFactor[i] * fPc[k]);

                    //fPc[k] = Math.Sign(beta - 0.5) * fPc[k] + Curves[0].dataCurve[i].pf;
                }
                //normalizedOfPowerofConsumed.Add(nPc);
                powerFactorofPowerofConsumed.Add(fPc);
            }

            //Reading the input data of buses
            buses = new Bus[busXML.ChildNodes.Count];
            for (int i = 0; i < buses.Length; i++)
            {
                buses[i] = new Bus();
                buses[i].Number = Convert.ToInt32(busXML.ChildNodes[i].Attributes["id"].Value);
                buses[i].S = new Complex[3];
                buses[i].S[0] = new Complex(Convert.ToDouble(busXML.ChildNodes[i].Attributes["pa"].Value) * 1e3, Convert.ToDouble(busXML.ChildNodes[i].Attributes["qa"].Value) * 1e3);
                buses[i].S[1] = new Complex(Convert.ToDouble(busXML.ChildNodes[i].Attributes["pb"].Value) * 1e3, Convert.ToDouble(busXML.ChildNodes[i].Attributes["qb"].Value) * 1e3);
                buses[i].S[2] = new Complex(Convert.ToDouble(busXML.ChildNodes[i].Attributes["pc"].Value) * 1e3, Convert.ToDouble(busXML.ChildNodes[i].Attributes["qc"].Value) * 1e3);
                buses[i].type = Convert.ToString(busXML.ChildNodes[i].Attributes["type"].Value);
            }

            //Reading the input data of branches
            branches = new Branch[branchXML.ChildNodes.Count];
            for (int i = 0; i < branches.Length; i++)
            {
                branches[i] = new Branch();
                branches[i].length_bran = Convert.ToInt32(branchXML.ChildNodes[i].Attributes["length"].Value);
                branches[i].switch_ = Convert.ToBoolean(Convert.ToInt32(branchXML.ChildNodes[i].Attributes["switch"].Value));
                branches[i].From = Convert.ToInt32(branchXML.ChildNodes[i].Attributes["from"].Value);
                branches[i].To = Convert.ToInt32(branchXML.ChildNodes[i].Attributes["to"].Value);

                branches[i].Z = new Complex[3][];
                for (int j = 0; j < 3; j++) branches[i].Z[j] = new Complex[3];

                //Phase A line impedance
                branches[i].Z[0][0] = new Complex(Convert.ToDouble(branchXML.ChildNodes[i].Attributes["raa"].Value), Convert.ToDouble(branchXML.ChildNodes[i].Attributes["xaa"].Value));
                branches[i].Z[0][1] = new Complex(Convert.ToDouble(branchXML.ChildNodes[i].Attributes["rab"].Value), Convert.ToDouble(branchXML.ChildNodes[i].Attributes["xab"].Value));
                branches[i].Z[0][2] = new Complex(Convert.ToDouble(branchXML.ChildNodes[i].Attributes["rac"].Value), Convert.ToDouble(branchXML.ChildNodes[i].Attributes["xac"].Value));
                //Phase B line impedance
                branches[i].Z[1][0] = branches[i].Z[0][1];
                branches[i].Z[1][1] = new Complex(Convert.ToDouble(branchXML.ChildNodes[i].Attributes["rbb"].Value), Convert.ToDouble(branchXML.ChildNodes[i].Attributes["xbb"].Value));
                branches[i].Z[1][2] = new Complex(Convert.ToDouble(branchXML.ChildNodes[i].Attributes["rbc"].Value), Convert.ToDouble(branchXML.ChildNodes[i].Attributes["xbc"].Value));
                //Phase C line impedance
                branches[i].Z[2][0] = branches[i].Z[0][2];
                branches[i].Z[2][1] = branches[i].Z[1][2];
                branches[i].Z[2][2] = new Complex(Convert.ToDouble(branchXML.ChildNodes[i].Attributes["rcc"].Value), Convert.ToDouble(branchXML.ChildNodes[i].Attributes["xcc"].Value));

                //Make the origin bus association 
                for (int j = 0; j < buses.Length; j++)
                {
                    if (branches[i].From == buses[j].Number)
                    {
                        branches[i].bFrom = buses[j];
                        break;
                    }
                }

                //Make the destination bus association
                for (int j = 0; j < buses.Length; j++)
                {
                    if (branches[i].To == buses[j].Number)
                    {
                        branches[i].bTo = buses[j];
                        break;
                    }
                }
            }

            //the probability of selection of the branch where the fault

            List<double> setOfLength = new List<double>();
            for (int i = 0; i < branches.Length; i++)
            {
                setOfLength.Add(branches[i].length_bran);
            }
            for (int i = 0; i < branches.Length; i++)
            {
                probabilityOfselectionofBranch_fault.Add(setOfLength[i] / setOfLength.Sum());
            }
            //frequenyMin and frequenyMax limit the range of the accumulated frequency of the chosen path j.
            double sumFrequency = 0.0;
            for (int i = 0; i < branches.Length; i++)
            {
                frequenyMin.Add(sumFrequency);
                sumFrequency += probabilityOfselectionofBranch_fault[i];
                frequenyMax.Add(sumFrequency);
            }
            //the discrete selection variable, Sel(γ), determines a path as a function of the random variable γ.

            for (int i = 0; i < branches.Length; i++)
            {

                if (frequenyMin[i] < randomVariable && randomVariable < frequenyMax[i])
                {
                    selectOfrandomVariable[i] = i;
                }
            }

            //select the type of fault: single phase, two phase or three phase.

            double probabilityOfselectionofsinglePhase_fault = 0.79;
            double probabilityOfselectionoftwoPhase_fault = 0.11;
            double probabilityOfselectionofThreePhase_fault = 0.02;
            double probabilityOfselectionofothers_fault = 0.08;

            probabilityOfEachTypeOfFault[0] = probabilityOfselectionofsinglePhase_fault;
            probabilityOfEachTypeOfFault[1] = probabilityOfselectionoftwoPhase_fault;
            probabilityOfEachTypeOfFault[2] = probabilityOfselectionofThreePhase_fault;
            probabilityOfEachTypeOfFault[3] = probabilityOfselectionofothers_fault;

            // probabilityOfselectofEachhourwithFault

            for (int i = 0; i < branches.Length; i++)
            {
                probabilityOfselectofEachhourwithFault.Add(1 / 24);
            }
            eneryNotSupplied_BS = 0;
            //energy not supplied
            for (int i = 0; i < buses.Length; i++)
            {
                for (int k = 0; k < 3; k++)
                {
                    eneryNotSupplied_BS += buses[i].S[k] * deltaTime_fault;
                }
            }
            eneryNotSupplied_AS = 0;
            for (int i = 0; i < branches.Length; i++)
            {
                if (branches[i].switch_ == false)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        eneryNotSupplied_AS += buses[i].S[k] * (deltaTime_fault - deltaTime_switch);
                    }
                }
            }

            eneryNotSupplied_total = eneryNotSupplied_AS + eneryNotSupplied_BS;
            //Functions

            //Make the upstream branch association
            for (int i = 0; i < branches.Length; i++)
            {
                for (int j = 0; j < branches.Length; j++)
                {
                    if (branches[i].From == branches[j].To)
                    {
                        branches[i].UpBranch = branches[j];
                        break;
                    }
                }
            }

            //Make the downstream branch association
            for (int i = 0; i < buses.Length; i++)
            {
                for (int j = 0; j < branches.Length; j++)
                {
                    if (buses[i].Number == branches[j].To)
                    {
                        buses[i].UpBranch = branches[j];
                    }
                    if (buses[i].Number == branches[j].From)
                    {
                        buses[i].DownBranches.Add(branches[j]);
                    }
                }
            }

            //Make the root bus association
            for (int i = 0; i < branches.Length; i++)
            {
                if (branches[i].UpBranch == null)
                {
                    root = branches[i].bFrom;
                    break;
                }
            }

            //creat the code of the feeder and thier assocations
            //obtainSections();
            //calculateCIC();







            for (int i = 0; i < generators.Length; i++)
            {
                generators[i].capacity = 0.0;

                // expoenent capacity
                generators[i].capacity = generators[i].minpowerActive + generators[i].capacity * generators[i].maxpowerActive;

            }




            for (int j = 0; j < buses.Length; j++)
            {
                buses[j].DGs = null;
                for (int k = 0; k < generators.Length; k++)
                {
                    if (buses[j].Number == generators[k].bFrombus)
                    {
                        buses[j].DGs = generators[k];
                        break;
                    }
                }

            }

        } //Read the network data from xml file 

        void DFS(Bus v)
        {
            //update the bus voltage 

            for (int i = 0; i < 3; i++)
            {
                v.V[i] = v.UpBranch.bFrom.V[i]; //Vi=Vu


                for (int j = 0; j < 3; j++)
                {
                    v.V[i] -= v.UpBranch.Z[i][j] * v.UpBranch.J[j];  //vi= vu-z.j
                }
            }

            // using sts for probabilistic insted of s(determistec)
            for (int i = 0; i < 3; i++)
            {
                SZIP[i] = v.S_ts[i] * (lambda_0 +
                            lambda_1 * Complex.Abs(v.V[i]) / Complex.Abs(Vse[i]) +
                            lambda_2 * Math.Pow(Complex.Abs(v.V[i]) / Complex.Abs(Vse[i]), 2));

                Complex s_DG = new Complex();
                if (v.DGs != null)
                {
                    s_DG = new Complex(v.DGs.capacity / 3, 0);


                }
                absolutePower.Add(Complex.Abs(v.V[i] * Complex.Conjugate(v.I[i]) - SZIP[i] + s_DG));
                v.I[i] = Complex.Conjugate((SZIP[i] - s_DG) / v.V[i]);

                // v.I[i] = Complex.Conjugate(v.S[i] / v.V[i]); //update Bus current injection Ii=Si/Vi
                v.UpBranch.J[i] = v.I[i]; //ji=Ii
            }

            for (int i = 0; i < v.DownBranches.Count; i++)
                DFS(v.DownBranches[i].bTo);

            //update the upstream curent flow
            if (v.UpBranch.UpBranch != null)
            {
                for (int i = 0; i < 3; i++)
                {
                    //v.UpBranch.UpBranch.J[i] += v.I[i]; //Ju=Ju+Ji=Ju+Ii 
                    v.UpBranch.UpBranch.J[i] += v.UpBranch.J[i];
                }
            }
        }

        public void calculatePowerFlow(int x, int max_x)
        {
            //int x = 0; //counter of simulation time
            //int max_x = 25; // number of simulation //6time
            while (x < max_x)
            {
                x++;
                //construct the 3 phase voltage and current vector
                for (int i = 0; i < buses.Length; i++)
                {
                    buses[i].V = new Complex[3];
                    buses[i].I = new Complex[3];
                    buses[i].S_ts = new Complex[3];

                }

                //construct the 3 phase current vector
                for (int j = 0; j < branches.Length; j++)
                {
                    branches[j].J = new Complex[3];
                }

                //The function of the percentage of energy consumption
                Random rand = new Random();
                double[] beta = new double[3];
                double[] beta_1 = new double[3];
                for (int i = 0; i < buses.Length; i++)
                {
                    //buses[i].loadpercentagedemand = new double[3];

                    double[] nPc = new double[3];
                    double[] fPc = new double[3];
                    Complex[] pEnergy = new Complex[3];
                    double[] pEnergyReal = new double[3];
                    int t_h;
                    t_h = x - (int)(x / 24) * 24;

                    for (int k = 0; k < 3; k++)
                    {
                        // LoadPercentegeDemandrese[k] = 0.0;
                        beta[k] = NextDouble(rand, 0, 1);

                        nPc[k] = Math.Sign(beta[k] - 0.5) * normalizedOfPowerofConsumed[t_h][k] + Curves[0].dataCurve[x].p;

                        fPc[k] = Math.Sign(beta_1[k] - 0.5) * powerFactorofPowerofConsumed[t_h][k] + Curves[0].dataCurve[x].pf;

                        pEnergy[k] = nPc[k] * new Complex(fPc[k], Math.Sin(Math.Acos(fPc[k])));

                        buses[i].S_ts[k] = buses[i].S[k].Magnitude * pEnergy[k];//  percentageEnergy[i][k];

                    }
                    percentageEnergy.Add(pEnergy);

                }

                //for (int i = 0; i < buses.Length; i++)
                //{
                //    buses[i].S_ts = new Complex[3];

                //    for (int k = 0; k < 3; k++)
                //    {
                //        buses[i].S_ts[k] = buses[i].S[k].Magnitude * percentageEnergy[i][k];

                //    }
                //}

                root.V = Vse;
                double epsilon = 1e-4;
                double mismatch = 1;
                while (mismatch > epsilon)
                {
                    //Update The Mismatch using eq4 

                    for (int i = 0; i < root.DownBranches.Count; i++)
                    {
                        DFS(root.DownBranches[i].bTo);
                    }

                    mismatch = absolutePower.Max();
                    absolutePower.Clear();

                }

                percentageEnergy.Clear();
                lossepoweractive = 0;
                for (int i = 0; i < branches.Length; i++)
                {
                    lossepoweractive += Math.Pow(branches[i].J[0].Magnitude, 2) * branches[i].Z[0][0].Real;
                }

                //print value of voltage and current by depend on period time x
                //printNetwork(System.IO.Directory.GetParent(@"../../").FullName + " \\tDN_135b_results" + x.ToString() + ".txt");

            }

        } //Calculate the power flow

        public void calculateCIC(XmlDocument xdoc)//calculate cic
        {
            void statusofcustomer()
            {
                //Define the status of sections witout neighbour feeder
                for (int i = 0; i < sections.Count; i++)
                {
                    bool changedSS = false; //flag indicating the status of any section was changed
                    for (int j = 0; j < sections.Count; j++)
                    {
                        if (sections[j].statusOfFault == true)
                        {
                            sections[j].statusOfSection = true;
                        }
                        if (sections[j].UpSection != null && sections[j].UpSection.statusOfSection == true && sections[j].statusOfSection == false)
                        {
                            sections[j].statusOfSection = true;

                            changedSS = sections[j].statusOfSection;
                        }
                    }
                    if (!changedSS) break;
                }

                //Check the neighbour feeder to finilize the definition of section state
                for (int j = 0; j < sections.Count; j++)
                {
                    double pdg = sections[j].DGs_section == null ? 0.0 : sections[j].DGs_section.capacity;

                    if (!sections[j].statusOfFault && useTS && sections[j].tie_switch_section != null && sections[j].tie_switch_section.capacity >= (sections[j].totalLoad.Sum() - pdg))
                    {
                        sections[j].statusOfSection = false;
                    }
                }

                //Check the DG to finilize the definition of section state
                for (int j = 0; j < sections.Count; j++)
                {
                    if (!sections[j].statusOfFault && sections[j].DGs_section != null && sections[j].DGs_section.capacity >= sections[j].totalLoad.Sum())
                    {
                        sections[j].statusOfSection = false;
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
                //for (int j = 0; j < sections.Count; j++)
                // {
                //  getSC(sections[j].rootBranch.bTo, sections[j]);
                // }
            }

            Random rnd = new Random(); //Randon variable
            int max_omega1 = 4;//set of time series comprised in the interruption  time interval;
            int max_run = 1000;//1000years    
            int run = 0;
            Random rand1 = new Random();
            Random rand2 = new Random();

            //saif = 6
            probabilityOffaultNetwork = 0.025 * totalLengthOffeeder * (1e-3);//probablity of fualt in network by year 

            for (int i = 0; i < sections.Count; i++)//probablity of selection section with fualt
            {
                sections[i].probabilityOfselectionofsection = sections[i].length_section / totalLengthOffeeder;
            }
            //frequenyMin and frequenyMax limit the range of the accumulated frequency of the chosen section i.
            double sumFrequency_section = 0.0;
            for (int i = 0; i < sections.Count; i++)
            {
                frequenyMin_section.Add(sumFrequency_section);
                sumFrequency_section += sections[i].probabilityOfselectionofsection;
                frequenyMax_section.Add(sumFrequency_section);
            }
            while (run < max_run)
            {
                Thread.Sleep(1);
                double randomVariable_net = NextDouble(rand1, 0, 1);
                if (probabilityOffaultNetwork > randomVariable_net) //compare probabilityOffaultNetwork > randomVariable_net
                {
                    costENSofNetwork = new double[3];
                    costICofNetwork = new double[3];
                    for (int r = 0; r < 6; r++)
                    {
                        FultHappenNetwork = true;
                        int omega1 = 0;
                        double randomVariable_Section1 = NextDouble(rand2, 0, 1);
                        int x = rnd.Next(0, Curves[0].dataCurve.Count - max_omega1 - 3);
                        //int x = 201; //for x at min or max load x=201(max) and x=6(min)
                        while (omega1 < max_omega1)
                        {
                            //Solve the power flow and get  the percentage of energy consumption
                            calculatePowerFlow(x, x + 1);

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

                            //costENSofNetwork = new double[3];

                            //Calculate the network C_ENS at time x
                            //the discrete selection variable, Sel(γ), determines a path as a function of the random variable γ.

                            for (int i = 0; i < sections.Count; i++)
                            {
                                if (frequenyMin_section[i] < randomVariable_Section1 && randomVariable_Section1 < frequenyMax_section[i])
                                {
                                    sections[i].statusOfFault = true;
                                    //total_load(sections[i].rootBranch.bTo, sections[i]);
                                    if (USEFLISER == false)
                                    {
                                        statusofcustomer(); //Default condition where downstream feeder section is interrupted
                                    }
                                    else
                                    {
                                        sections[i].statusOfSection = true; //Just faulted section is interrupted
                                    }
                                    //Put the code to calculate the Cost of ENS

                                    double customerDamageFunctionResidential;//for 4 hours
                                    double customerDamageFunctionCommercial;//for 4 hours
                                    double customerDamageFunctionIndustrial;//for 4 hours

                                    if (omega1 != 0 || USEFLISER == true)
                                    {
                                        //  electricityRatesResidential = 7.0;//residentialS  restorCost
                                        //electricityRatesCommercialS = 200.0;//commercial  restorCost
                                        //  electricityRatesIndustrial = 30.0;//industrial    restorCost
                                        customerDamageFunctionResidential = 1.8126 / (3);//for 4 hours
                                        customerDamageFunctionCommercial = 166.2123 / (3);//for 4 hours
                                        customerDamageFunctionIndustrial = 46.3678 / (3);//for 4 hours
                                    }
                                    else
                                    {
                                        //  electricityRatesResidential = 0.55;//residential  SwitchingCost
                                        //  electricityRatesCommercialS = 40.0;//commercial   SwitchingCost
                                        //  electricityRatesIndustrial = 5.0;//industrial     SwitchingCost
                                        customerDamageFunctionResidential = 0.55;//for fisr hours
                                        customerDamageFunctionCommercial = 40;//for first hours
                                        customerDamageFunctionIndustrial = 5.0;//for first hours
                                        for (int j = 0; j < sections.Count; j++)
                                            sections[j].statusOfSection = true;
                                    }


                                    double[] costICofResidential = new double[3];
                                    double[] costICOfCommercial = new double[3];
                                    double[] costICOfIndustrial = new double[3];
                                    //costIC = new double[3];

                                    //cost of energy not supplied for jth section in the ith time step($);furmula(6)
                                    double[] costENSofResidential = new double[3];
                                    double[] costENSOfCommercial = new double[3];
                                    double[] costENSOfIndustrial = new double[3];
                                    costENS = new double[3];
                                    for (int j = 0; j < sections.Count; j++)
                                    {
                                        if (sections[j].statusOfSection)
                                        {
                                            total_load(sections[j].rootBranch.bTo, sections[j]);

                                            for (int k = 0; k < 3; k++)
                                            {
                                                //CostIC=calculation of interruption cost of the j  th customer the formula 8 and 

                                                costICofResidential[k] += customerDamageFunctionResidential * sections[j].residentiallLoad[k] * 1e-3;
                                                costICOfCommercial[k] += customerDamageFunctionCommercial * sections[j].commercialLoad[k] * 1e-3;
                                                costICOfIndustrial[k] += customerDamageFunctionIndustrial * sections[j].industrialLoad[k] * 1e-3;
                                                //  costENSofResidential[k] += costICofResidential[k] * sections[j].residentiallLoad[k] * 1e-3;
                                                // costENSOfCommercial[k] += costICOfCommercial[k] * sections[j].commercialLoad[k] * 1e-3;
                                                // costENSOfIndustrial[k] += costICOfIndustrial[k] * sections[j].industrialLoad[k] * 1e-3;
                                                costENS[k] = costICofResidential[k] + costICOfCommercial[k] + costICOfIndustrial[k];
                                                costENSofNetwork[k] += costENS[k];
                                            }
                                        }
                                    }
                                    //Clear the sections atributes associated to fault
                                    sections[i].statusOfFault = false;
                                    for (int j = 0; j < sections.Count; j++) sections[j].statusOfSection = false;
                                    break;
                                }
                            }
                            //Clear the load of section
                            for (int i = 0; i < sections.Count; i++)
                            {
                                sections[i].residentiallLoad = new double[3];
                                sections[i].commercialLoad = new double[3];
                                sections[i].industrialLoad = new double[3];
                                sections[i].totalLoad = new double[3];
                            }
                            omega1++;
                            x += omega1;
                        }
                    }
                    list_CENSnetwork.Add(costENSofNetwork);
                }
                run++;
            }
        }
        //sectionlising

        public void obtainSections(XmlDocument xdoc)//obtainSection
        {
            //defind the sectionalizing branch, lenght section, demand load



            void Sectionlizing_branch(Bus advanced, Section newsection)
            {
                //Associate a tie switch
                for (int j = 0; j < tieSwitches.Length; j++)
                {
                    if (advanced.Number == tieSwitches[j].bTo_tieswitch)
                    {
                        newsection.tie_switch_section = tieSwitches[j];
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
            root.DownBranches[0].switch_ = true;
            //define the number of sections
            for (int i = 0; i < branches.Length; i++)
            {
                if (branches[i].switch_ == true)
                {
                    Section newSection = new Section();
                    newSection.rootBranch = branches[i];
                    newSection.residentiallLoad = new double[3];
                    newSection.industrialLoad = new double[3];
                    newSection.commercialLoad = new double[3];
                    newSection.totalLoad = new double[3];
                    newSection.length_section = 0;
                    sections.Add(newSection);
                    newSection.statusOfFault = false;
                    newSection.statusOfSection = false;
                    // newSection.failureRate = 0.025;
                    newSection.failureRate = 0.0;

                }
            }

            //Make the  Sectionalizing Branches Association
            for (int j = 0; j < sections.Count; j++)
            {
                Sectionlizing_branch(sections[j].rootBranch.bTo, sections[j]);



            }
            totalLengthOffeeder = 0;
            //Make the  Sectionalizing Branches Association
            for (int j = 0; j < sections.Count; j++)
            {

                totalLengthOffeeder += sections[j].length_section;


            }
            //Make the upstream section association
            for (int i = 0; i < sections.Count; i++)
            {
                //failure rate of section j
                sections[i].failureRate = sections[i].length_section / totalLengthOffeeder;

                for (int j = 0; j < sections.Count; j++)
                {
                    for (int k = 0; k < sections[j].SectionalizingBranches.Count; k++)
                    {
                        if (sections[j].SectionalizingBranches[k] == sections[i].rootBranch)
                        {
                            sections[i].UpSection = sections[j];
                            break;
                        }
                    }
                    if (sections[i].UpSection != null)
                    {
                        break;
                    }
                }
            }

        }

        public void printNetwork(string _path)
        {
            //Create a file to write to.
            bool printPF = false; //print voltages and currents
            bool printHD = false; //print historical data
            bool printCENS = true; //print C_ENS
            bool printIC = false; //print C_ENS
            using (StreamWriter sw = new StreamWriter(_path))
            {
                int count = 1; //Counter of output data types
                if (printPF)
                {
                    double Vse = 13800 / Math.Sqrt(3);//base voltage

                    //Voltage values
                    sw.WriteLine("<<<" + count.ToString() + ". Voltage(pu)>>>-----------------------------------------------------------------------------------------");
                    for (int i = 0; i < buses.Length - 1; i++)
                        sw.WriteLine("{0}: {1}[{2}°] {3}[{4}°] {5}[{6}°]",
                            buses[i].Number.ToString(),
                            (buses[i].V[0].Magnitude / Vse).ToString(),
                            (buses[i].V[0].Phase * 180 / Math.PI).ToString(),
                            (buses[i].V[1].Magnitude / Vse).ToString(),
                            (buses[i].V[1].Phase * 180 / Math.PI).ToString(),
                            (buses[i].V[2].Magnitude / Vse).ToString(),
                            (buses[i].V[2].Phase * 180 / Math.PI).ToString());
                    count++;

                    //Current Values
                    sw.WriteLine("<<<" + count.ToString() + ". Current(A)>>>----------------------------------------------------------------------------------------");
                    for (int i = 0; i < branches.Length; i++)
                        sw.WriteLine("{0}: {1}[{2}°] {3}[{4}°] {5}[{6}°]",
                            branches[i].bTo.Number.ToString(),
                            branches[i].J[0].Magnitude.ToString(),
                            (branches[i].J[0].Phase * 180 / Math.PI).ToString(),
                            branches[i].J[1].Magnitude.ToString(),
                            (branches[i].J[1].Phase * 180 / Math.PI).ToString(),
                            branches[i].J[2].Magnitude.ToString(),
                            (branches[i].J[2].Phase * 180 / Math.PI).ToString());
                    count++;
                }
                if (printHD)
                {
                    // R_TH Values
                    sw.WriteLine("<<<" + count.ToString() + ". R_TH (hourly measurements of the normalized magnitude)>>>----------------------------------------------------------------------------------------");
                    for (int j = 0; j < Curves[0].dataCurve.Count / 24; j++)


                        sw.WriteLine("{0}: 1`:{1}  2`:{1}  3`:{2}  4`:{3}  5`:{4}  6`:{5}  7`:{6}  8`:{7}  9`:{8}  10`:{9}  11`:{10}  12`:{11}  13`:{12}  14`:{13}  15`:{14}  16`:{15}  17`:{16}  18`:{17}  19`:{18}  20`:{19}  21`:{20}  22`:{21} 23`:{22} 24`:{23}",
                        j + 1,
                         R_TH[0][j].ToString(), R_TH[1][j].ToString(), R_TH[2][j].ToString(), R_TH[3][j].ToString(), R_TH[4][j].ToString(),
                        R_TH[5][j].ToString(), R_TH[6][j].ToString(), R_TH[7][j].ToString(), R_TH[8][j].ToString(), R_TH[9][j].ToString(),
                        R_TH[10][j].ToString(), R_TH[11][j].ToString(), R_TH[12][j].ToString(), R_TH[13][j].ToString(), R_TH[14][j].ToString(),
                        R_TH[15][j].ToString(), R_TH[16][j].ToString(), R_TH[17][j].ToString(), R_TH[18][j].ToString(), R_TH[16][j].ToString(),
                        R_TH[20][j].ToString(), R_TH[21][j].ToString(), R_TH[22][j].ToString(), R_TH[23][j].ToString()
                        );
                    count++;

                    // F_TH Values
                    sw.WriteLine("<<<" + count.ToString() + ". F_TH (hourly measurements of power factor)>>>----------------------------------------------------------------------------------------");
                    for (int j = 0; j < Curves[0].dataCurve.Count / 24; j++)


                        sw.WriteLine("{0}: 1`:{1}  2`:{1}  3`:{2}  4`:{3}  5`:{4}  6`:{5}  7`:{6}  8`:{7}  9`:{8}  10`:{9}  11`:{10}  12`:{11}  13`:{12}  14`:{13}  15`:{14}  16`:{15}  17`:{16}  18`:{17}  19`:{18}  20`:{19}  21`:{20}  22`:{21} 23`:{22} 24`:{23}",
                       j + 1,
                      F_TH[0][j].ToString(), F_TH[1][j].ToString(), F_TH[2][j].ToString(), F_TH[3][j].ToString(), F_TH[4][j].ToString(),
                      F_TH[5][j].ToString(), F_TH[6][j].ToString(), F_TH[7][j].ToString(), F_TH[8][j].ToString(), F_TH[9][j].ToString(),
                       F_TH[10][j].ToString(), F_TH[11][j].ToString(), F_TH[12][j].ToString(), F_TH[13][j].ToString(), F_TH[14][j].ToString(),
                      F_TH[15][j].ToString(), F_TH[16][j].ToString(), F_TH[17][j].ToString(), F_TH[18][j].ToString(), F_TH[16][j].ToString(),
                        F_TH[20][j].ToString(), F_TH[21][j].ToString(), F_TH[22][j].ToString(), F_TH[23][j].ToString()
                        );
                    count++;
                }
                if (printCENS)
                {
                    //CENS values
                    sw.WriteLine("<<<" + count.ToString() + ". printCENS>>>-----------------------------------------------------------------------------------------");
                    for (int i = 0; i < list_CENSnetwork.Count; i++)
                        sw.WriteLine("{0} {1} {2}               {3},",

                            list_CENSnetwork[i][0].ToString(),
                            list_CENSnetwork[i][1].ToString(),
                            list_CENSnetwork[i][2].ToString(),
                            list_CENSnetwork[i].Sum()
                            );
                    count++;
                }
                if (printIC)
                {
                    //CENS values
                    sw.WriteLine("<<<" + count.ToString() + ". printCENS>>>-----------------------------------------------------------------------------------------");
                    for (int i = 0; i < list_ICnetwork.Count; i++)
                        sw.WriteLine("{0} {1} {2}               {3}",

                            list_ICnetwork[i][0].ToString(),
                            list_ICnetwork[i][1].ToString(),
                            list_ICnetwork[i][2].ToString(),
                            list_ICnetwork[i].Sum()
                            );
                    count++;
                }

            }
        }//Save distribution network data into a text file

    }

}

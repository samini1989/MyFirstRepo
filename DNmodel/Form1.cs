using DNmodel.Common;
using DNmodel.GA;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Xml;

namespace DNmodel
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void Runbtn_Click(object sender, EventArgs e)
        {

            //Object for input data
            XmlDocument Xdoc = new XmlDocument();
            Xdoc.Load(System.IO.Directory.GetParent(@"../../").FullName + "\\tDN_135b_V2.xml");


            //Object of the distribuion network (DN)
            Network DN = new Network();

            //Reading the input data of distribution network
            DN.readNetwork(Xdoc);

            //Reading  data of section
            DN.obtainSections(Xdoc);

            //Reading  data of caculating cost
            DN.calculateCIC(Xdoc);

            //Culculate Power Flow
            //DN.calculatePowerFlow(8,12); //8 is the starting time and 12 is the stopping time

            //print the power flow results
            DN.printNetwork(System.IO.Directory.GetParent(@"../../").FullName + "\\tDN_135b_results.txt");
            //DN.printNetwork(@"C:\Users\RVargas\Documents\Visual Studio 2012\Projects\ConsoleAppAli\ConsoleAppAli\LoadCurve.xml");

            MessageBox.Show("Run Finished.");
        }

        private void RunGAbtn_Click(object sender, EventArgs e)
        {
            //Object for input data
            XmlDocument Xdoc = new XmlDocument();
            Xdoc.Load(System.IO.Directory.GetParent(@"../../").FullName + "\\tDN_135b_V2.xml");


            //Object of the distribuion network (DN)

            Network DN = new Network();

            //Reading the input data of distribution network
            DN.readNetwork(Xdoc);
            //GA
            GeneticAlgorithm geneticAlguritm = new GeneticAlgorithm();

            geneticAlguritm.energetic = DN;


            geneticAlguritm.metodoAg(System.IO.Directory.GetParent(@"../../").FullName + "\\tDN_135b_resultsAG3.txt");
            //  AG.melhor.decodifica(aG.energisa, aG.cp, aG.cp2, aG.n_facts, aG.custofacts_A, System.IO.Directory.GetParent(@"../../").FullName);
            //   geneticAlguritm.best.decodification(geneticAlguritm.energetic, geneticAlguritm.cp, geneticAlguritm.CostofEns, System.IO.Directory.GetParent(@"../../").FullName);

            //getin the best solution variable
            DG[] list_dg = new DG[1];

            Switch[] list_sw = new Switch[1];
            geneticAlguritm.best.getcoded(DN, ref list_dg, ref list_sw);

            //Reading  data of section
            // DN.obtainSections(Xdoc);

            //Reading  data of caculating costs
            // DN.calculateCIC(Xdoc);

            //Culculate Power Flow
            //DN.calculatePowerFlow(8,12); //8 is the starting time and 12 is the stopping time

            //print the power flow results
            DN.printNetwork(System.IO.Directory.GetParent(@"../../").FullName + "\\tDN_135b_results.txt");

            //DN.printNetwork(@"C:\Users\RVargas\Documents\Visual Studio 2012\Projects\ConsoleAppAli\ConsoleAppAli\LoadCurve.xml");

            MessageBox.Show("Run GA Finished");
        }
    }
}

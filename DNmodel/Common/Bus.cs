using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace DNmodel.Common
{
    public class Bus
    {
        //Input data
        public int Number; //Bus identification (ID)
        public Complex[] S; // Complex power vector

        public string consumptionType;

        //Class association
        public List<Branch> DownBranches = new List<Branch>(); //Downstream branches
        public Branch UpBranch; //Upstream Branch
                                //Power flow variables
        public Complex[] V; //Voltage vector
        public Complex[] I; //Bus current injection
        public Complex[] S_ts; //Load power by time

        public string type;//type of customer of bus
        public bool statusofthejthcustomer;// = new Boolean[3];// binary variable that indicates the status of the jthcustomer.If equal to 1, the energy supply is interrupted
        public DG DGs;
    }
}

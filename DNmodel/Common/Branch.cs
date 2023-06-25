using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace DNmodel.Common
{
    public class Branch
    {
        //Input data
        public int From; //Origin bus (From)
        public int To; // Destination bus (to)
        public Complex[][] Z; //Branch impedance matrix
        public int length_bran;//lenght of each branch
        public Boolean switch_;
        //Class associations
        public Bus bFrom; //Relation with origin bus
        public Bus bTo; //Relation with destination bus
        public Branch UpBranch; //Relation with upstream branch 

        public List<TieSwitch> tie_switch = new List<TieSwitch>();

        //Power Flow variables
        public Complex[] J; //Current flow through the branch

    }
}

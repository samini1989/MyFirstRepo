using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DNmodel.Common
{
    public class TieSwitch
    {
        public double capacity;
        //Class associations
        public int bFrom_tieswitch; //Relation with origin bus
        public int bTo_tieswitch; //Relation with destination bus

    }
}

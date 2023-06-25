using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DNmodel.Common
{
    public class DG
    {
        public double capacity;
        //Class associations
        public int bFrombus; //Relation with origin bus
        public int bTo_tieswitch; //Relation with destination bus

        public double maxpowerActive = 0;// 3e6;
        public double minpowerActive = 0;// 100e3;

    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DNmodel.Common
{
    public class Section
    {
        //Input data
        public double[] commercialLoad;
        public double[] industrialLoad;
        public double[] residentiallLoad;
        public double failureRate;
        public int length_section;//lenght of section
        public double[] totalLoad;

        public Boolean statusOfFault;    //status of the fault of the section 
        public Boolean statusOfSection; //the status(energized or de-energized) of the other sections

        //Class association
        public Section UpSection; //Relation with upstream section
        public List<Branch> SectionalizingBranches = new List<Branch>(); //SectioningBranches
        public Branch rootBranch; //Root branch of the distribution network

        //public TieSwitch
        public TieSwitch tie_switch_section;
        public double[] powerLossSection;

        public double probabilityOfselectionofsection;
        public double probabilityOffaultofsection;

        public DG DGs_section;
    }

}

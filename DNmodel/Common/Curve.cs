using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DNmodel.Common
{
    public class Curve
    {
        //Class assocation
        public List<DataCurve> dataCurve = new List<DataCurve>(); //
        public double[] varianceOfHourlyNormalized;//the variance of hourly normalized magnitude,
        public double[] varianceOfHourlyOfPowerFactor;// and the variance of hourly power factor
        public double[] deviationOfHourlyOfPowerFactor;
        public double[] deviationOfHourlyNormalized;
    }

}

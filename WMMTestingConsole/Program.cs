using System;
using System.Collections.Generic;
using System.Linq;
using WMM_Csharp;

namespace WMMTestingConsole
{
    public class Program
    {
        public static void Main(string[] args)
        {
            Console.Write("Latitude: ");
            double lat = double.Parse(Console.ReadLine());
            Console.Write("Longitude: ");
            double longi = double.Parse(Console.ReadLine());
            Console.Write("Elevation above MSL (km): ");
            double elevation = double.Parse(Console.ReadLine());
            Console.Write("Time (decimal years): ");
            double time = double.Parse(Console.ReadLine());
            var me = MagneticFieldCalculator.CalculateMagneticElements(lat, longi, elevation, time);


            if(time < MagneticFieldCalculator.epoch || time > MagneticFieldCalculator.epoch + 5)
            {
                Console.WriteLine($"Warning: Date Provided is outside of the available epoch.\nCalculation Errors may occur.\nFor Best Results enter dates between {MagneticFieldCalculator.epoch} and {MagneticFieldCalculator.epoch + 5}");
                Console.WriteLine("An updated Coefficient list may be available from the NCEI: https://www.ncei.noaa.gov/products/world-magnetic-model");
            }

            if(me.H < 2000)
            {
                Console.WriteLine("Warning: location is in the blackout zone around the magnetic pole as defined by the WMM military specification (https://www.ngdc.noaa.gov/geomag/WMM/data/MIL-PRF-89500B.pdf). Compass accuracy is highly degraded in this region.");
            } else if (me.H >= 2000 && me.H < 6000)
            {
                Console.WriteLine("“Caution: location is approaching the blackout zone around the magnetic pole as defined by the WMM military specification (https://www.ngdc.noaa.gov/geomag/WMM/data/MIL-PRF-89500B.pdf). Compass accuracy may be degraded in this region.");
            }

            Console.WriteLine("Results:");

            Console.WriteLine($@"Field   
F    = {Math.Round(me.F, 3)} nT
H    = {Math.Round(me.H, 3)} nT
X    = {Math.Round(me.X, 3)} nT
Y    = {Math.Round(me.Y, 3)} nT
Z    = {Math.Round(me.Z, 3)} nT
Decl = {Math.Round(me.D, 3)}°
Incl = {Math.Round(me.I, 3)}°
");

            Console.WriteLine("\nPress Enter to close program.");

            Console.ReadLine();
        }
    }
}

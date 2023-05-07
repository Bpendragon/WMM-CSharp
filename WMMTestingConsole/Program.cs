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
            MagneticFieldCalculator.CalculateMagneticElements(-80.0, 240.0, 100.0, 2022.5);
            Console.ReadLine();
        }
    }
}

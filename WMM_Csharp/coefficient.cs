using System;

namespace WMM_Csharp
{
    public class Coefficient
    {
        //All of these are at t_0 (or more specifically year 2020.0)
        public int N { get; private set; } //Index
        public int M { get; private set; } //Order
        public double G { get; private set; } //Main Field Coefficient
        public double H { get; private set; } //Main Field Coefficient
        public double Gdot { get; private set; } //Secular Variation Coefficient 
        public double Hdot { get; private set; } //Secular Variation Coefficient 

        public Coefficient (string InputLine)
        {
            var vals = InputLine.Split(' ', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries);
            if (vals.Length < 6) throw new FormatException($"{nameof(InputLine)} is not a valid WMM Coefficient, too few values provided.\nInput as given:\n{InputLine}");
            if (vals.Length > 6) throw new FormatException($"{nameof(InputLine)} is not a valid WMM Coefficient, too many values provided.\nInput as given:\n{InputLine}");
            N = int.Parse(vals[0]);   
            M = int.Parse(vals[1]);
            G = double.Parse(vals[2]);
            H = double.Parse(vals[3]);
            Gdot = double.Parse(vals[4]);
            Hdot = double.Parse(vals[5]);
        }

        public override string ToString()
        {
            return $"n: {N}, m: {M}, g: {G}, h: {H}, g_dot: {Gdot}, h_dot: {Hdot}";
        }
    }
}
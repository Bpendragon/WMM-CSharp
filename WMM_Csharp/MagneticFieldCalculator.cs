using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Linq.Expressions;
using System.Net.NetworkInformation;
using System.Security.Cryptography.X509Certificates;
using System.Runtime.CompilerServices;

namespace WMM_Csharp
{
    public static class MagneticFieldCalculator
    {
        private static Dictionary<(int n, int m), Coefficient> WMMCoefficients;
        private const double RADIUS = 6371200.0; //Earth Radius, meters, value taken directly from WMM annual report
        private const double WGS84_RADIUS = 6378137.0;
        private const double WGS84_RF = 298.257223563;
        private static double epoch;

        static MagneticFieldCalculator()
        {
            if(WMMCoefficients == null)
            {
                //Create Coefficient list from file.
                WMMCoefficients = new();
                var lines = File.ReadAllLines("./WMM_COF/WMM.COF");

                epoch = double.Parse(lines[0].Split(' ', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)[0]);

                foreach (var l in lines.Skip(1).Where(a => !a.StartsWith("99999")))
                {
                    var c = new Coefficient(l);
                    WMMCoefficients[(c.N, c.M)] =  c;
                }

            }
        }

        /// <summary>
        /// Calculates the Magnetic Elements for a given location and altitude
        /// </summary>
        /// <param name="latitude">Geodetic Latitude in Decimal Degrees</param>
        /// <param name="longitude">Longitude in Decimal Degrees</param>
        /// <param name="altitude">Height above Mean Sea Level, in km</param>
        /// <param name="time">Time to calculate for in Decimal Years (i.e. January 1, 2023 = 2023.0, while October 31st, 2023 ~= 2023.83)</param>
        /// <returns></returns>
        public static MagneticElements CalculateMagneticElements(double latitude, double longitude, double altitude, double time)
        {
            // Use names used in the report for ease of use/finding
            // Note: report uses actual Greek letters, here transcribed into ascii for obvious reasons
            // Throughout this method variables preceded with an underscore are only needed in the current step, those without are important in the next step(s).
            double lambda = DegToRad(longitude);
            double phi = DegToRad(latitude);
            double h_MSL = altitude * 1000;
            double t = time;

            MagneticElements me = new();
            
            
            //Step 1: Convert Geodetic Latitude to Spherical.
            var _Rc = Rc(phi);
            var _p = p(_Rc, h_MSL, phi);
            var _z = z(_Rc, h_MSL, phi);

            var r = r0(_p, _z);
            var phi_prime = phiP(_z, r);

            //Step 2: Calculate Gauss Coefficients for the given Time

            Dictionary<(int n, int m), (double g, double h)> gCof = new();
            Dictionary<(int n, int m), double> Schmidts = new(); //pre-compute all the schmidt semi-normalized legendre polynomials to make future steps easier

            foreach (var kvp in WMMCoefficients)
            {
                double _g = kvp.Value.G + (t - epoch) * kvp.Value.Gdot;
                double _h = kvp.Value.H + (t - epoch) * kvp.Value.Hdot;

                gCof[kvp.Key] = (_g, _h);
                Schmidts[kvp.Key] = SchmidtNorm(kvp.Key.n, kvp.Key.m, phi_prime);
            }

            for(int m = 0; m <= 13; m++)
            {
                Schmidts[(13, m)] = SchmidtNorm(13, m, phi_prime); //Needed for calculating the derivitives
            }
            
            //Step 3: Calculate Field Vector Components in Geocentric Coordinates
            double X_Prime = me.X_prime = XPrime(lambda, phi_prime, r, Schmidts, gCof);
            double Y_Prime = me.Y_prime = YPrime(lambda, phi_prime, r, Schmidts, gCof);
            double Z_Prime = me.Z_prime = ZPrime(lambda, phi_prime, r, Schmidts, gCof);

            //Step 3.5: Calculate Secular Variation of the field components in Geocentric Coordinates.
            double X_Prime_Dot = me.X_prime_dot = XPrimeDot(lambda, phi_prime, r, Schmidts);
            double Y_Prime_Dot = me.Y_prime_dot = YPrimeDot(lambda, phi_prime, r, Schmidts);
            double Z_Prime_Dot = me.Z_prime_dot = ZPrimeDot(lambda, phi_prime, r, Schmidts);

            //Step 4: Rotate X', Y', and Z' into the ellipsoidal reference frame.
            double X = me.X = (X_Prime * Math.Cos(phi_prime - phi)) - (Z_Prime * Math.Sin(phi_prime - phi));
            double Y = me.Y = Y_Prime;
            double Z = me.Z = (X_Prime * Math.Sin(phi_prime - phi)) + (Z_Prime * Math.Cos(phi_prime - phi));

            //Step 4.5: Rotate Prime-Dots into the ellipsoidal reference frame
            double X_dot = me.X_dot = (X_Prime_Dot * Math.Cos(phi_prime - phi)) - (Z_Prime_Dot * Math.Sin(phi_prime - phi));
            double Y_dot = me.Y_dot = Y_Prime_Dot;
            double Z_dot = me.Z_dot = (X_Prime_Dot * Math.Sin(phi_prime - phi)) + (Z_Prime_Dot * Math.Cos(phi_prime - phi));

            //Step 5: Calculate Magnetic Field Elements from the orthogonal components
            double H = me.H = Math.Sqrt(Math.Pow(X, 2) + Math.Pow(Y, 2));
            double F = me.F = Math.Sqrt(Math.Pow(H, 2) + Math.Pow(Z, 2));
            double I = Math.Atan2(Z, H); 
            double D = Math.Atan2(Y, X);

            //Step 5.5 Calculate Secular Variation for these
            double H_dot = me.H_dot = ((X * X_dot) + (Y * Y_dot)) / H;
            me.F_dot = ((X * X_dot) + (Y * Y_dot) + (Z * Z_dot)) / F;
            double I_dot = ((H * Z_dot) - (Z * H_dot)) / Math.Pow(F, 2);
            double D_dot = ((X * Y_dot) - (Y * X_dot)) / Math.Pow(H, 2);

            //Step 6: Convert Appropriate units from Radians to Degrees (or Radians/Year -> Degrees/Year)
            me.I = RadToDeg(I);
            me.D = RadToDeg(D);
            me.I_dot = RadToDeg(I_dot);
            me.D_dot = RadToDeg(D_dot);

            return me;
        }

        static double f => 1 / WGS84_RF;
        static double e_2 => f * (2 - f);

        static double Rc   (double phi) => WGS84_RADIUS / (Math.Sqrt(1 - (e_2 * Math.Pow(Math.Sin(phi), 2))));
        static double p    (double Rc, double h, double phi) => (Rc + h) * Math.Cos(phi);
        static double z    (double Rc, double h, double phi) => (Rc * (1 - e_2) + h) * Math.Sin(phi);
        static double r0   (double p, double z) => Math.Sqrt(Math.Pow(p, 2) + Math.Pow(z, 2));
        static double phiP (double z, double r) => Math.Asin(z / r); 

        static public double DegToRad (double theta) => (theta * Math.PI) / 180;
        static public double RadToDeg (double theta) => (theta * 180) / Math.PI;

        // All of the following functions are written very explicitly. The vast majority of the arithmetic can be crunched down into a series of one-liners.
        // THAT SAID: This is a lot easier to read and debug. Plus the Compiler optimizes them down already, so no reason not to.
        private static double XPrime(double lambda, double phi_prime, double r, Dictionary<(int n, int m), double> Schmidts, Dictionary<(int n, int m), (double g, double h)> gCof)
        {
            double sum = 0;
            for(int n = 1; n <= 12; n++)
            {
                double inner = 0;

                for(int m = 0; m <= n; m++)
                {
                    double tmp = gCof[(n, m)].g * Math.Cos(m * lambda);
                    tmp += gCof[(n, m)].h * Math.Sin(m * lambda);
                    tmp *= SchmidtDeriv(n, m, phi_prime, Schmidts);
                    inner += tmp;
                }

                inner *= Math.Pow(RADIUS / r, n + 2);

                sum += inner;
            }

            return -sum;
        }

        private static double YPrime(double lambda, double phi_prime, double r, Dictionary<(int n, int m), double> Schmidts, Dictionary<(int n, int m), (double g, double h)> gCof)
        {

            double sum = 0;
            for (int n = 1; n <= 12; n++)
            {
                double inner = 0;

                for (int m = 0; m <= n; m++)
                {
                    inner += ((gCof[(n, m)].g * Math.Sin(m * lambda)) - (gCof[(n, m)].h * Math.Cos(m * lambda))) * Schmidts[(n, m)] * m;
                }

                inner *= Math.Pow(RADIUS / r, n + 2);

                sum += inner;
            }

            return sum / Math.Cos(phi_prime);
        }

        private static double ZPrime(double lambda, double phi_prime, double r, Dictionary<(int n, int m), double> schmidts, Dictionary<(int n, int m), (double g, double h)> gCof)
        {
            double sum = 0;

            for(int n = 1; n<= 12; n++)
            {
                double inner = 0;
                for (int m = 0; m <= n; m++)
                {
                    double tmp = (gCof[(n, m)].g * Math.Cos(m * lambda));
                    tmp += gCof[(n, m)].h * Math.Sin(m * lambda);
                    tmp *= schmidts[(n, m)];
                    inner += tmp;
                }
                inner *= (n + 1) * Math.Pow(RADIUS / r, n + 2);
                sum += inner;
            }

            return -sum;
        }

        private static double XPrimeDot(double lambda, double phi_prime, double r, Dictionary<(int n, int m), double> Schmidts)
        {
            double sum = 0;
            for (int n = 1; n <= 12; n++)
            {
                double inner = 0;

                for (int m = 0; m <= n; m++)
                {
                    double tmp = WMMCoefficients[(n, m)].Gdot * Math.Cos(m * lambda);
                    tmp += WMMCoefficients[(n, m)].Hdot * Math.Sin(m * lambda);
                    tmp *= SchmidtDeriv(n, m, phi_prime, Schmidts);
                    inner += tmp;
                }

                inner *= Math.Pow(RADIUS / r, n + 2);

                sum += inner;
            }

            return -sum;
        }

        private static double YPrimeDot(double lambda, double phi_prime, double r, Dictionary<(int n, int m), double> Schmidts)
        {
            double sum = 0;
            for (int n = 1; n <= 12; n++)
            {
                double inner = 0;

                for (int m = 0; m <= n; m++)
                {
                    inner += ((WMMCoefficients[(n, m)].Gdot * Math.Sin(m * lambda)) - (WMMCoefficients[(n, m)].Hdot * Math.Cos(m * lambda))) * Schmidts[(n, m)] * m;
                }

                inner *= Math.Pow(RADIUS / r, n + 2);

                sum += inner;
            }

            return sum / Math.Cos(phi_prime);
        }

        private static double ZPrimeDot(double lambda, double phi_prime, double r, Dictionary<(int n, int m), double> schmidts)
        {
            double sum = 0;

            for (int n = 1; n <= 12; n++)
            {
                double inner = 0;
                for (int m = 0; m <= n; m++)
                {
                    double tmp = (WMMCoefficients[(n, m)].Gdot * Math.Cos(m * lambda));
                    tmp += WMMCoefficients[(n, m)].Hdot * Math.Sin(m * lambda);
                    tmp *= schmidts[(n, m)];
                    inner += tmp;
                }
                inner *= (n + 1) * Math.Pow(RADIUS / r, n + 2);
                sum += inner;
            }

            return -sum;
        }

        static double AssocLegendrePolynomial(int n, int m, double phi_prime)
        {
            var t = Math.Sin(phi_prime);

            //Function taken from Heiskanen Moritz 1967 Physical Geodesy Equation 1-62

            //Integer Math Abuse
            int r = (n - m) / 2;

            double sum = 0;
            for(int k = 0; k <= r; k++)
            {
                double tmp = factorial((2 * n) - (2 * k));
                tmp /= factorial(k) * factorial(n - k) * factorial(n - m - (2 * k));
                tmp *= Math.Pow(-1, k);
                tmp *= Math.Pow(t, n - m - (2 * k));
                sum += tmp;
            }

            sum *= Math.Pow(2, -n) * Math.Pow(1 - Math.Pow(t, 2), m / 2.0); //Remember to force doubles for math.pow if needed.

            return sum;
        }

        static double SchmidtNorm (int n, int m, double phi_prime)
        {
            var lp = AssocLegendrePolynomial(n, m, phi_prime);
            if (m == 0) return lp;


            double res = lp * Math.Sqrt(2 * (factorial(n - m) / factorial(n + m)));
            return res;
        }

        static double SchmidtDeriv (int n, int m, double phi_prime, Dictionary<(int n, int m), double> Schmidts)
        {
            double res = n + 1;
            res *= Math.Tan(phi_prime);
            res *= Schmidts[(n, m)];

            double res2 = 1 / Math.Cos(phi_prime);
            res2 *= Schmidts[(n + 1, m)];
            res2 *= Math.Sqrt(Math.Pow(n + 1, 2) - Math.Pow(m, 2));

            return res - res2;
        }

        private static double factorial(int k)
        {
            if (k == 0) return 1.0; //Base Case
            if (k < 0) throw new ArgumentOutOfRangeException();
            return k * factorial(k - 1);
        }

    }
}

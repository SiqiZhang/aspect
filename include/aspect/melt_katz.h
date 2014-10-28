#ifndef MELT_KATZ_H
#define MELT_KATZ_H
#include "aspect/bisection.h"

namespace Melt_Katz
{
    struct Parameters
    {
        double A1;
        double A2;
        double A3;

        double B1;
        double B2;
        double B3;

        double C1;
        double C2;
        double C3;

        double r1;
        double r2;

        double beta1;
        double beta2;

        double K;
        double gamma;
        double D_H2O;

        double X1;
        double X2;
        double lambda;

        double Cp;
        double deltaS;
    };
    class Melt_Katz//: public Bisection
    {
        friend double Bisection::bisecion_solve<Melt_Katz>(Melt_Katz a,double min_value,double max_value,double max_error,bool screen_out);

private:
        double T;
        double P;
        double Mcpx;
        double X;
        double error_fun(double T1)
        {
            double F=melt_fraction.get_melt_fraction(T1,P,Mcpx,X);
            double T2=T*parameters.Cp/(parameters.Cp+F*parameters.deltaS);
            return T2-T1;
        };
        const struct Parameters &parameters;
public:
        Melt_Katz(const struct Parameters &parameters):
        parameters(parameters),
        melt_fraction(parameters)
        {
        };
        class Melt_fraction
        {
            friend double Bisection::bisecion_solve<Melt_fraction>(Melt_fraction a,double min_value,double max_value,double max_error,bool screen_out);

public:
            Melt_fraction(const struct Parameters &parameters):
                parameters(parameters)
            {
            };

            /*
             * Get the melt fraction at given temperature, pressure, cpx weight fraction,
             * and water weight fraction within the melt. The water weight fraction within
             * the melt is a function of melt fraction itself, that is why it has to be solved iteraivly.
             */
            double melt_fraction_x_melt(double T, double P, double Mcpx, double X_H2O_melt)const;

            /*
             * Get the melt fraction at given temperature, pressure, cpx weight fraction,
             * and bulk water weight fraction. Using bisction method.
             * X_H2O in unit wt%
             */
            double get_melt_fraction(double T, double P, double Mcpx, double X_H2O);

            /*
             * Solidut temperature
             */
            double T_solidus(double P) const;

            /*
             * Liquidus temperature
             */
            double T_liquidus(double P) const;

            /*
             * Lherz liquidus temperature
             */
            double T_lherz(double P) const;

            /*
             * Solidus depression due to fraction of water in the melt.
             * XH2O is the weight fraction of water in the melt.
             */
            double get_deltaT(double XH2O, double P) const;

            double get_Xs(double P) const;

            double get_X(double X_bulk,double melt_fraction) const;

        private:
            double T;
            double P;
            double Mcpx;
            double X_H2O;
            /*
             *     Parameters for melting of peridotite after Katz, 2003
             */
            const struct Parameters &parameters;

            //void read_parameters(char* filename);

            double error_fun(double F)
            {
                return melt_fraction_x_melt(T,P,Mcpx,get_X(X_H2O,F))-F;
            };

        } melt_fraction;

        /*
         * Get the modified temperature after melting. Using bisction method and
         * calling get_melt_fraction() function.
         */
        double get_modified_temperature(double T0,double P, double Mcpx, double X_H2O,
                                        double &melt);


    };

}

#endif // MELT_KATZ_H

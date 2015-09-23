#ifndef __aspect__melting_h
#define __aspect__melting_h

#include <string>
#include <vector>

namespace aspect
{
  namespace melting
  {
    class Melting_curve
    {
      public:
        Melting_curve();
        void read(const std::string &filename);
        double T(const double p, const double radius) const;
        bool is_radius;
        unsigned int Num_points;
      private:
        //unsigned int Num_points;
        std::vector<double> T_array;
        std::vector<double> P_array;
        //bool is_radius;
    };
    class Melting_data
    {
      public:
        Melting_data();
        void read(const std::string solidus_file, const std::string liquidus_file);
        Melting_curve Solidus;
        Melting_curve Liquidus;
        double dT;
        double dP;
        double Melting_fraction(const double T, const double p, const double radius,
                                const double water, const double depletion) const;
        double get_solidus(const double p, const double radius,const double water, const double depletion) const;
        double get_liquidus(const double p, const double radius,const double water, const double depletion) const;
        double get_melt_fraction_derivative_temperature(double T, double P, double water, double depletion) const;
        double get_melt_fraction_derivative_pressure(double T, double P, double water, double depletion) const;
    };
  }

}

#endif

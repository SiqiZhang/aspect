/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/
/*  $Id: table.h 1624 2013-04-28 21:15:28Z heister $  */


#ifndef __aspect__model_melt_h
#define __aspect__model_melt_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/melting.h>
#include <aspect/melt_katz.h>

#include <aspect/material_model/steinberger.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class Melt: public MaterialModel::Steinberger<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const;

        virtual double viscosity_ratio (const double temperature,
                                        const double pressure,
                                        const std::vector<double>    &compositional_fields,
                                        const SymmetricTensor<2,dim> &strain_rate,
                                        const Point<dim> &position) const;

        virtual double density (const double temperature,
                                const double pressure,
                                const std::vector<double> &compositional_fields,
                                const Point<dim> &position) const;

        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
        * Return true if the viscosity() function returns something that
        * may depend on the variable identifies by the argument.
        */
        virtual bool
        viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
        * Return true if the density() function returns something that
        * may depend on the variable identifies by the argument.
        */
        virtual bool
        density_depends_on (const NonlinearDependence::Dependence dependence) const;


        /**
         * @name Auxiliary material properties used for postprocessing
         * @{
         */

        virtual void revise_composition (const double temperature,
                                         const double pressure,
                                         const std::vector<double> &compositional_fields,
                                         const Point<dim> &position,
                                         std::vector<double> &new_compositional_fields) const;

      /**
       * Return the product of the change in entropy across phase
       * transitions, the pressure derivative of the phase function
       * (if this is the pressure derivative) or the product of the
       * former two and the Clapeyron slope (if this is the
       * temperature derivative).
       * The entropy change across a phase transition can be calculated
       * as $\frac{\gamma \Delta\rho}{\rho_\text{light} \rho_\text{heavy}}$.
       * $\gamma$ is the Clapeyron slope of the phase transition,
       * $\Delta\rho$ is the density jump across the phase transition,
       * $\rho_\text{light}$ is the density of the light material
       * (above the phase transition) and $\rho_\text{heavy}$ the
       * density of the heavy material (below the phase transition).
       * The phase function hat values ranging from 0 to 1 indicating
       * which percentage of the material has already undergone the phase
       * transition. Its argument is usually the excess pressure
       * $\pi = p - p_0 - \gamma T$, where $p_0$ is the zero-degree
       * transition pressure.
       *
      * This function has a default implementation that sets
       * the entropy gradient to zero (assuming no phase changes).
       */
      virtual double entropy_derivative (const double      temperature,
                                         const double      pressure,
                                         const std::vector<double> &compositional_fields,
                                         const Point<dim> &position,
                                         const NonlinearDependence::Dependence dependence) const;
      /*
       * Deal with depletion with reaction_term of compositional field.
       */
      virtual double reaction_term (const double temperature,
                                    const double pressure,
                                    const std::vector<double> &compositional_fields,
                                    const Point<dim> &position,
                                    const unsigned int compositional_variable) const;


      /*
       * Get the changed temperature after melt.
       * model given by Katz et al. (2003) and
       * influnces of Cpx and water are included.
       */
      double new_temperature_after_melt(const double temperature,
                                        const double pressure,
                                        const std::vector<double> &compositional_fields,
                                        const Point<dim> &position) const;

      /*
       * Get the melt fraction iterat through melting relation, 
       * model given by Katz et al. (2003) and
       * influnces of Cpx and water are included.
       */
      double melt_fraction(const double temperature,
                           const double pressure,
                           const std::vector<double> &compositional_fields,
                           const Point<dim> &position) const;

      /**
       * @}
       */

      /**
       * @name Functions used in dealing with run-time parameters
       * @{
       */
      /**
       * Declare the parameters this class takes through input files.
       */
      static
      void
      declare_parameters (ParameterHandler &prm);

      /**
       * Read the parameters this class declares from the parameter
       * file.
       */
      virtual
      void
      parse_parameters (ParameterHandler &prm);
      /**
       * @}
       */

      virtual
      void
      initialize();

      aspect::melting::Melting_data Data_Melt;
    
    double get_viscosity_diffusion_inverse   (double temperature, double pressure, double depth) const;
    double get_viscosity_dislocation_inverse (double temperature, double pressure, double strain_rate_II, double depth) const;
    double get_viscosity_yield_inverse       (double pressure, const std::vector<double> &compositional_fields,double strain_rate_II, double depth) const;
    double get_viscosity_peierls_inverse     (double temperature, double pressure, double strain_rate_II, double depth) const;

      private:
    static const double R_gas = 8.341; //Gas constant.
    struct Melt_Katz::Parameters melting_parameters;
    int i_composition_Cpx;
    int i_composition_H2O;
    double default_Cpx;

    double exponential_melt;
    double depth_lower;
    double viscosity_cutoff_low;
    double viscosity_cutoff_high;
    std::vector<double> density_difference;
    std::vector<double> viscosity_difference;
    
    bool                is_yield_enable;
    bool                is_yield_dependent_on_composition;
    std::vector<double> yield_stress_surface;
    std::vector<double> yield_friction;
    std::vector<double> yield_composition_factor;
    std::vector<double> yield_stress_max;

    bool   is_diffusion_enable;
    double activation_energy_diffusion_um;
    double activation_volume_diffusion_um;
    double prefactor_diffusion_um;
    double activation_energy_diffusion_lm;
    double activation_volume_diffusion_lm;
    double prefactor_diffusion_lm;        
    
    bool   is_dislocation_enable;
    double activation_energy_dislocation;
    double activation_volume_dislocation;
    double prefactor_dislocation;
    double stress_exponent_dislocation;

    bool   is_peierls_enable;
    double activation_energy_peierls;
    double activation_volume_peierls;
    double prefactor_peierls;
    double stress_exponent_peierls;


    double k_value;
    std::string solidus_filename;
    std::string liquidus_filename;
    double Lh;
        
    double get_viscosity_arrhenius (double temperature, double pressure, double strain_rate_II, double A, double E, double V, double n) const;

    };
  }
}

#endif

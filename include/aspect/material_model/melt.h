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

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that reads the essential values of coefficients from
     * tables in input files that describe their dependence as a function of
     * pressure and temperature.
     *
     * @ingroup MaterialModels
     */
	
	namespace internal
	{
		class MaterialLookup;
	};
	


	template <int dim>
    class Melt: public MaterialModel::InterfaceCompatibility<dim>,
	            public SimulatorAccess<dim>
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

        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const std::vector<double> &compositional_fields,
                                        const Point<dim> &position) const;

        virtual double specific_heat (const double temperature,
                                      const double pressure,
                                      const std::vector<double> &compositional_fields,
                                      const Point<dim> &position) const;

        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const std::vector<double> &compositional_fields,
                                             const Point<dim> &position) const;
        /**
         * @}
         */
/*
       virtual double reaction_term (const double temperature,
                                     const double pressure,
                                     const std::vector<double> &compositional_fields,
                                     const Point<dim> &position,
                                     const unsigned int compositional_variable) const;
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
        * Return true if the compressibility() function returns something that
        * may depend on the variable identifies by the argument.
        *
        * This function must return false for all possible arguments if the
        * is_compressible() function returns false.
        */
        virtual bool
        compressibility_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
        * Return true if the specific_heat() function returns something that
        * may depend on the variable identifies by the argument.
        */
        virtual bool
        specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
        * Return true if the thermal_conductivity() function returns something that
        * may depend on the variable identifies by the argument.
        */
        virtual bool
        thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure but is only interpreted as such in the rhs of the Stokes equation.
         * This is consistent with the so called Bousinesq formulation.
         * In the current context, compressibility means whether we should solve the contuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */

        /**
        * A reference viscosity
        *
        * The value here is not used in the computation of things but only in
        * postprocessing the solution when we want dimension-less
        * quantities.
        */
        virtual double reference_viscosity () const;

        /**
         * A reference density
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less
         * quantities.
         */
        virtual double reference_density () const;

        /**
         * A reference thermal diffusivity $\kappa$. $\kappa$ is related to the thermal
         * conductivity $k$ as $\kappa = k/(rho c_p)$.
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less
         * quantities.
         */
        double reference_thermal_diffusivity () const;

        /**
         * A reference thermal expansion coefficient $\alpha$.
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less
         * quantities.
         */
        double reference_thermal_expansion_coefficient () const;

        /**
         * Return the thermal expansion coefficient $\alpha$ of the model,
         * possibly as a function of depth.
         */
        virtual double thermal_expansion_coefficient (const double temperature,
                                                      const double pressure,
                                                      const std::vector<double> &compositional_fields,
                                                      const Point<dim> &position) const;

        /**
        * A reference thermal specific heat $c_p$.
        *
        * The value here is not used in the computation of things but only in
        * postprocessing the solution when we want dimension-less
        * quantities.
        */
        double reference_cp () const;
        /**
         * @}
         */

        /**
         * @name Auxiliary material properties used for postprocessing
         * @{
         */

        /**
         * the seismic pressure wave speed
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less
         * quantities.
         */
        virtual double seismic_Vp (const double temperature,
                                   const double pressure,
                                   const std::vector<double> &compositional_fields,
								   const Point<dim> &position) const;

        /**
         * the seismic shear wave speed
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less
         * quantities.
         */
        virtual double seismic_Vs (const double temperature,
                                   const double pressure,
                                   const std::vector<double> &compositional_fields,
								   const Point<dim> &position) const;

        /**
         * the phase of the composition at given pressure and temperature
         * this returns an integer value that is associated with a specific phase
         * for instance Majorite of PostPerovskite
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less
         * quantities.
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

      private:
    struct Melt_Katz::Parameters melting_parameters;
    int i_composition_Cpx;
    int i_composition_H2O;
    double default_Cpx;

    double reference_rho;
    double reference_T;
    double reference_kappa;
    double reference_specific_heat;
    double reference_alpha;
    std::string composition;
    std::string data_directory;
    bool compute_phases;
    bool model_is_compressible;

    std::string viscosity_model;
    double reference_eta;
    double exponential_T;
    double exponential_P;
		double exponential_melt;
		double viscosity_factor_litho;
		double viscosity_factor_trans;
		double viscosity_factor_lower;
    double depth_litho;
		double depth_trans;
		double depth_lower;
		double viscosity_cutoff_low;
		double viscosity_cutoff_high;
		//double yield_stress;
    //double yield_stress_increase;
    double reference_dT;
		double mantle_thickness;
		double earth_radius;
    std::vector<double> yield_stress;
    std::vector<double> yield_stress_increase;
    std::vector<double> density_difference;
    std::vector<double> viscosity_difference;
    std::vector<double> yield_factor;

		double increase_lower_mantle;
    double activation_energy_diffusion_um;
    double activation_volume_diffusion_um;
    double prefactor_diffusion_um;
    double activation_energy_diffusion_lm;
    double activation_volume_diffusion_lm;
    double prefactor_diffusion_lm;        
    double activation_energy_dislocation;
    double activation_volume_dislocation;
    double prefactor_dislocation;
    double stress_exponent;

    double k_value;
		std::string solidus_filename;
		std::string liquidus_filename;
		double Lh;
        
		bool interpolation;
		double latent_heat;
		std::string datadirectory;
    std::vector<std::string> material_file_names;
    unsigned int n_material_data;
    std::vector<std_cxx1x::shared_ptr<internal::MaterialLookup> > material_lookup;
		double get_deltat (const Point<dim> &position) const;
		
    };
  }
}

#endif

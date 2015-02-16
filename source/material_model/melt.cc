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
/*  $Id: table.cc 1790 2013-06-16 22:14:04Z bangerth $  */


#include <aspect/material_model/melt.h>
#include <aspect/melting.h>
#include <aspect/simulator_access.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/table.h>
#include <deal.II/base/symmetric_tensor.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include <algorithm>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
  
  
    template <int dim>
    void
    Melt<dim>::
    initialize()
    {
      Steinberger<dim>::initialize();
      Data_Melt.read(solidus_filename,liquidus_filename);
    }

    template <int dim>
    double
    Melt<dim>::
    viscosity (const double temperature,
               const double pressure,
               const std::vector<double> &composition,
               const SymmetricTensor<2,dim> &strain_rate,
               const Point<dim> &position) const
    {
      double viscosity;
      const double depth = (this->get_geometry_model()).depth(position);

      const double strain_rate_II=std::max(std::sqrt(std::fabs(second_invariant(strain_rate))),1e-17);

      const double viscosity_diffusion_inverse   = get_viscosity_diffusion_inverse(temperature,pressure,depth);
      const double viscosity_dislocation_inverse = get_viscosity_dislocation_inverse(temperature,pressure,strain_rate_II,depth);
      const double viscosity_yield_inverse       = get_viscosity_yield_inverse(pressure,composition,strain_rate_II,depth);
      const double viscosity_peierls_inverse     = get_viscosity_peierls_inverse(temperature,pressure,strain_rate_II,depth);
      viscosity=1./(viscosity_diffusion_inverse + viscosity_dislocation_inverse + viscosity_yield_inverse + viscosity_peierls_inverse);

       // Apply cutoff
       viscosity = std::max(viscosity,viscosity_cutoff_low);
       viscosity = std::min(viscosity,viscosity_cutoff_high);

       /*
       double radius=sqrt(position.square());
       Melt_fraction=Data_Melt.Melting_fraction(temperature,pressure,radius,0.,0.);
       viscosity*=std::exp(-std::log(exponential_melt)*Melt_fraction);
       */

       // Apply compositional difference, it can go outside cutoff
       if(viscosity_difference.size()>0)
         for(unsigned i=0;i<composition.size();i++)
           viscosity*=pow(viscosity_difference[i],composition[i]);
       return viscosity;
    }

    template <int dim>
    double
    Melt<dim>::
    viscosity_ratio (const double temperature,
                     const double pressure,
                     const std::vector<double> &composition,
                     const SymmetricTensor<2,dim> &strain_rate,
                     const Point<dim> &position) const
    {
      const double depth = (this->get_geometry_model()).depth(position);
      const double strain_rate_II=std::max(std::sqrt(std::fabs(second_invariant(strain_rate))),1e-17);

      std::vector<double> viscosity_inverse(4);
      viscosity_inverse[0] = get_viscosity_diffusion_inverse(temperature,pressure,depth);
      viscosity_inverse[1] = get_viscosity_dislocation_inverse(temperature,pressure,strain_rate_II,depth);
      viscosity_inverse[2] = get_viscosity_yield_inverse(pressure,composition,strain_rate_II,depth);
      viscosity_inverse[3] = get_viscosity_peierls_inverse(temperature,pressure,strain_rate_II,depth);
      std::vector<double>::iterator biggest = std::max_element(viscosity_inverse.begin(), viscosity_inverse.end());
      // Return which one is the largest in four inversed viscosity values.
      return std::distance(viscosity_inverse.begin(), biggest);
    }

    template <int dim>
    double
    Melt<dim>::
    density (const double temperature,
             const double pressure,
             const std::vector<double> &compositional_fields,
             const Point<dim> &position) const
    {
      double rho = Steinberger<dim>::density(temperature,pressure,compositional_fields,position);
      if(density_difference.size()!=0 && density_difference.size()==compositional_fields.size())
      {
        for(unsigned i=0;i<compositional_fields.size();i++)
        {
          rho+=density_difference[i]*compositional_fields[i];
        }
      }
      else
      {
        Assert (density_difference.size()==compositional_fields.size(),
            ExcMessage("Size of density differences is different from size of compositional fields."
                       " This can not be intended."));
      }
                
      return rho;
    }

    template <int dim>
    bool
    Melt<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return ((dependence & NonlinearDependence::pressure)
              ||
              (dependence & NonlinearDependence::temperature)
              ||
              (dependence & NonlinearDependence::strain_rate)
              ||
              (dependence & NonlinearDependence::compositional_fields));
    }


    template <int dim>
    bool
    Melt<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return ((dependence & NonlinearDependence::pressure)
              ||
              (dependence & NonlinearDependence::temperature)
              ||
              (dependence & NonlinearDependence::compositional_fields) );
    }

    template <int dim>
    void
    Melt<dim>::
    revise_composition (const double temperature,
                        const double pressure,
                        const std::vector<double> &compositional_fields,
                        const Point<dim> &position,
                        std::vector<double> &new_compositional_fields) const
    {
      Melt_Katz::Melt_Katz melt_calculation(melting_parameters);
      new_compositional_fields.resize(compositional_fields.size());
      double Mcpx,X_H2O,fraction;
      for(unsigned int i=0;i<compositional_fields.size();i++)
        new_compositional_fields[i]=compositional_fields[i];
      if(i_composition_Cpx>=0 && i_composition_Cpx<(int)compositional_fields.size())
        Mcpx=compositional_fields[i_composition_Cpx];
      else
        Mcpx=0.;
      if(i_composition_H2O>=0 && i_composition_H2O<(int)compositional_fields.size())
        X_H2O=compositional_fields[i_composition_H2O];
      else
        X_H2O=0.;
      fraction=melt_calculation.melt_fraction.get_melt_fraction(temperature,pressure,Mcpx,X_H2O);

      if(i_composition_Cpx>=0 && i_composition_Cpx<(int)compositional_fields.size())
        new_compositional_fields[i_composition_Cpx]-=melt_calculation.melt_fraction.get_Mcpx(pressure,Mcpx,fraction)*fraction;
      if(i_composition_H2O>=0 && i_composition_H2O<(int)compositional_fields.size())
        new_compositional_fields[i_composition_H2O]-=std::min(melt_calculation.melt_fraction.get_X(X_H2O,fraction),
                                                              melt_calculation.melt_fraction.get_Xs(pressure))*fraction;
      
/*
			cout<<"Old composition:"<<compositional_fields[0]<<","
				                    <<compositional_fields[1]<<","
									<<compositional_fields[2]<<std::endl
				<<"New composition:"<<new_compositional_fields[0]<<","
				                    <<new_compositional_fields[1]<<","
									<<new_compositional_fields[2]<<std::endl;
                                                                        */
      //Make sure modified compositional field within 0~100
      for(unsigned int i=0;i<3;i++)
      {
        if(new_compositional_fields[i]<0.)
          new_compositional_fields[i]=0.;
        if(new_compositional_fields[i]>1.)
          new_compositional_fields[i]=1.;
      }
      return;
    }

    template <int dim>
    double
    Melt<dim>::
    reaction_term (const double temperature,
                   const double pressure,
                   const std::vector<double> &compositional_fields,
                   const Point<dim> &position,
                   const unsigned int compositional_variable) const
    {
      double delta_C=0.;
      /*
      Melt_Katz::Melt_Katz melt_calculation(melting_parameters);
      double Mcpx,X_H2O,fraction;
      if(i_composition_Cpx>=0 && i_composition_Cpx<(int)compositional_fields.size())
        Mcpx=std::max(0.,compositional_fields[i_composition_Cpx]);
      else
        Mcpx=0.;
      if(i_composition_H2O>=0 && i_composition_H2O<(int)compositional_fields.size())
        X_H2O=std::max(0.,compositional_fields[i_composition_H2O]);
      else
        X_H2O=0.;
      fraction=melt_calculation.melt_fraction.get_melt_fraction(temperature,pressure,Mcpx,X_H2O);
      if((int)compositional_variable==i_composition_Cpx)
        delta_C = -melt_calculation.melt_fraction.get_Mcpx(pressure,Mcpx,fraction);
      if((int)compositional_variable==i_composition_H2O)
        delta_C = -std::min(melt_calculation.melt_fraction.get_X(X_H2O,fraction),
                            melt_calculation.melt_fraction.get_Xs(pressure))*fraction;
      if(delta_C>0. || delta_C+compositional_fields[compositional_variable]<0.)
        cout<<"["<<compositional_variable<<"]P="<<pressure<<",T="<<temperature<<",F="<<fraction<<",Mcpx="<<Mcpx<<",delta_C="<<delta_C<<endl;
      */
      return delta_C;
    }

    template <int dim>
    double
    Melt<dim>::
    entropy_derivative (const double temperature,
                        const double pressure,
                        const std::vector<double> &compositional_fields,
                        const Point<dim> &position,
                        const NonlinearDependence::Dependence dependence) const
    {
      Melt_Katz::Melt_Katz melt_calculation(melting_parameters);
      double Mcpx,X_H2O,fraction;
      if(i_composition_Cpx>=0 && i_composition_Cpx<(int)compositional_fields.size())
          Mcpx=std::max(0.,compositional_fields[i_composition_Cpx]);
      else
          Mcpx=default_Cpx;
      if(i_composition_H2O>=0 && i_composition_H2O<(int)compositional_fields.size())
          X_H2O=std::max(0.,compositional_fields[i_composition_H2O]);
      else
          X_H2O=0.;
     if (dependence == NonlinearDependence::temperature)
       return melt_calculation.melt_fraction.get_melt_entropy_derivative_temperature(
           temperature,pressure,Mcpx,X_H2O);
     else if (dependence == NonlinearDependence::pressure)
       return melt_calculation.melt_fraction.get_melt_entropy_derivative_pressure(temperature,pressure,Mcpx,X_H2O);
     return 0.;
    }

    template <int dim>
    double
    Melt<dim>::
    new_temperature_after_melt(const double temperature,
                               const double pressure,
                               const std::vector<double> &compositional_fields,
                               const Point<dim> &position) const
    {
      Melt_Katz::Melt_Katz melt_calculation(melting_parameters);
      double Mcpx,X_H2O,fraction;
      if(i_composition_Cpx>=0 && i_composition_Cpx<(int)compositional_fields.size())
        Mcpx=compositional_fields[i_composition_Cpx];
      else
        Mcpx=default_Cpx;
      if(i_composition_H2O>=0 && i_composition_H2O<(int)compositional_fields.size())
        X_H2O=compositional_fields[i_composition_H2O];
      else
        X_H2O=0.;
      return melt_calculation.get_modified_temperature(temperature,pressure,Mcpx,X_H2O,fraction);
    }

    template <int dim>
    double
    Melt<dim>::
    melt_fraction(const double temperature,
                  const double pressure,
                  const std::vector<double> &compositional_fields,
                  const Point<dim> &position) const
    {
      Melt_Katz::Melt_Katz melt_calculation(melting_parameters);
      double Mcpx,X_H2O,fraction;
      if(i_composition_Cpx>=0 && i_composition_Cpx<(int)compositional_fields.size())
          Mcpx=std::max(0.,compositional_fields[i_composition_Cpx]);
      else
          Mcpx=default_Cpx;
      if(i_composition_H2O>=0 && i_composition_H2O<(int)compositional_fields.size())
          X_H2O=std::max(0.,compositional_fields[i_composition_H2O]);
      else
          X_H2O=0.;
      //melt_calculation.get_modified_temperature(temperature,pressure,Mcpx,X_H2O,fraction);
      fraction=melt_calculation.melt_fraction.get_melt_fraction(temperature,pressure,Mcpx,X_H2O);
      return fraction;
    }

    template <int dim>
    double
    Melt<dim>::get_viscosity_arrhenius (double temperature, double pressure, double strain_rate_II, double A, double E, double V, double n) const
    {
      return pow(A,-1./n)
            *exp((E+pressure*V)/(n*R_gas*temperature))
            *pow(strain_rate_II,(1.-n)/n);
    }

    template <int dim>
    double
    Melt<dim>::get_viscosity_diffusion_inverse (double temperature, double pressure, double depth) const
    {
      if(is_diffusion_enable)
      {
        if(depth<depth_lower)
          return 1./get_viscosity_arrhenius(temperature,pressure,1.,
                                            prefactor_diffusion_um,  
                                            activation_energy_diffusion_um,
                                            activation_volume_diffusion_um,
                                            1.);
        else
          return 1./get_viscosity_arrhenius(temperature,pressure,1.,        
                                            prefactor_diffusion_lm,         
                                            activation_energy_diffusion_lm, 
                                            activation_volume_diffusion_lm, 
                                            1.);
      }
      else
        return 0.;
    }

    template <int dim>
    double
    Melt<dim>::get_viscosity_dislocation_inverse (double temperature, double pressure, double strain_rate_II, double depth) const
    {
      if(is_dislocation_enable && depth<depth_lower)
        return 1./get_viscosity_arrhenius(temperature,pressure,strain_rate_II,
                                          prefactor_dislocation,
                                          activation_energy_dislocation,
                                          activation_volume_dislocation,
                                          stress_exponent_dislocation);
      else
        return 0.;
    }

    template <int dim> 
    double
    Melt<dim>::get_viscosity_peierls_inverse (double temperature, double pressure, double strain_rate_II, double depth) const
    {
      if(is_peierls_enable && depth<depth_lower)
        return 1./get_viscosity_arrhenius(temperature,pressure,strain_rate_II,
                                         prefactor_peierls,
                                         activation_energy_peierls,
                                         activation_volume_peierls,
                                         stress_exponent_peierls);
      else
        return 0.;
    }

    template <int dim> 
    double
    Melt<dim>::get_viscosity_yield_inverse (double pressure,
                                            const std::vector<double>    &compositional_fields,
                                            double strain_rate_II,
                                            double depth) const
    {
      if(is_yield_enable && depth<depth_lower)
      {
        unsigned int n_compositional_fields=compositional_fields.size();
        double yield_stress=0.;
        double default_composition=1.;
        if(is_yield_dependent_on_composition)
        {
          for(unsigned i=0;i<n_compositional_fields;i++)
          {
            double composition_i=std::min(1.,std::max(compositional_fields[i]*yield_composition_factor[i+1],0.));
            yield_stress += (yield_stress_surface[i+1]+yield_friction[i+1]*pressure)*composition_i;
            default_composition -= composition_i;
          }
          default_composition = std::min(std::max(0.,default_composition),1.);
          yield_stress += (yield_stress_surface[0]+yield_friction[0]*pressure)
                          *default_composition*yield_composition_factor[0];
        }
        else
        {
          yield_stress=yield_stress_surface[0]+yield_friction[0]*pressure;        
        }
        return 1./yield_stress*strain_rate_II*2.0;
      }
      else
        return 0.;
    }

    template <int dim>
    void
    Melt<dim>::declare_parameters (ParameterHandler &prm)
    {
      Steinberger<dim>::declare_parameters(prm);
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt model");
        {
          prm.enter_subsection ("Viscosity");
          {
            prm.declare_entry ("Viscosity cutoff low","1e19",
                               Patterns::Double (0),
                               "The lowest viscosity cut off Unit: Pa s");
            prm.declare_entry ("Viscosity cutoff high","1e24",
                               Patterns::Double (0),
                                "The largest viscosity cut off Unit: Pa s");
            prm.declare_entry ("Exponential Melt", "1e0",
                               Patterns::Double (0),
                               "Multiplication factor or melting fraction exponent");
             prm.declare_entry ("Lower mantle depth","660e3",
                               Patterns::Double (0),
                               "The depth of lower mantle top.");
             prm.declare_entry ("Density difference", "",
                                Patterns::List(Patterns::Double ()),
                                "Density difference of different composition (kg/m^3)");
             prm.declare_entry ("Viscosity factor", "",
                                Patterns::List(Patterns::Double ()),
                                "Viscosity difference of different composition (Pa s)");
             //Yield stress
             prm.declare_entry ("Enable yield", "false",
                                Patterns::Bool (),
                                "Whether the model consider yield. ");
             prm.enter_subsection ("Yield stress");
             {
               
               prm.declare_entry ("Composition factor", "0",
                                  Patterns::List(Patterns::Double (0)),
                                  "Scale factor of different composition for calculate yield stress."
                                  "It contains list of n_compositional_field+1 factors,"
                                  "the fist one is for default composition (0 no yield)");
               prm.declare_entry ("Yield stress surface", "1.17e8",
                                  Patterns::List(Patterns::Double ()),
                                  "Yield stress on surface for different materials, the first one is the "
                                  "default material. There should be n_compostion+1 values in the "
                                  "list. (Pa)");
               prm.declare_entry ("Yield friction", "0.1",
                                  Patterns::List(Patterns::Double (0)),
                                  "Yield friction coefficient for different materials, "
                                  "the first one is for the default material. There should be "
                                  "n_compostion+1 values in the list.");
               prm.declare_entry ("Yield stress max","10e9",
                                  Patterns::List(Patterns::Double ()),
                                  "The max yield stress for different materials, "
                                  "the first one is for the default material. There should be "
                                  "n_compostion+1 values in the list.");
             }
             prm.leave_subsection();

             //Diffusion creep
             prm.declare_entry ("Enable diffusion", "false",
                                Patterns::Bool (),
                                "Whether the model consider diffusion creep. ");
             prm.enter_subsection ("Diffusion creep");
             {
               prm.declare_entry ("E_UM", "335e3",
                                  Patterns::Double (0),
                                  "Activation energy for diffusion creep of upper mantle."
                                  "Unit: J/mol");
               prm.declare_entry ("V_UM", "4.0e-6",
                                  Patterns::Double (0),
                                  "Activation volume for diffusion creep of upper mantle."
                                  "Unit: m^3/mol");
               prm.declare_entry ("A_UM", "1.92e-11",
                                  Patterns::Double (0),
                                  "Prefactor for diffusion creep of upper mantle."
                                  "Yita=1/A*(E+P*V)/(R*T)");

               prm.declare_entry ("E_LM", "335e3",
                                  Patterns::Double (0),
                                  "Activation energy for diffusion creep of lower mantle"
                                  "Unit: J/mol");
               prm.declare_entry ("V_LM", "4.0e-6",
                                  Patterns::Double (0),
                                  "Activation volume for diffusion creep of lower mantle"
                                  "Unit: m^3/mol");
               prm.declare_entry ("A_LM", "1.92e-11",
                                  Patterns::Double (0),
                                  "prefactor for diffusion creep of lower mantle"
                                  "Yita=1/A*(E+P*V)/(R*T)");
             }
             prm.leave_subsection();

             // Dislocation creep
             prm.declare_entry ("Enable dislocation","false",
                                Patterns::Bool (),
                                "Whether the model consider dislocation creep. ");
             prm.enter_subsection ("Dislocation creep");
             {
               prm.declare_entry ("E", "335e3",
                                  Patterns::Double (0),
                                  "Activation energy for dislocation creep"
                                  "Unit: J/mol");
               prm.declare_entry ("V", "4.0e-6",
                                  Patterns::Double (0),
                                  "Activation volume for dislocation creep"
                                  "Unit: m^3/mol");
               prm.declare_entry ("A", "1.92e-11",
                                  Patterns::Double (0),
                                  "Prefactor for dislocation creep "
                                  "Yita=A^(-1/n)*(E+P*V)/(n*R*T)*Strain_rate_II^((1-n)/n)");
               prm.declare_entry ("n", "3.5",
                                  Patterns::Double (0),
                                  "Stress exponent for dislocation creep");
             }
             prm.leave_subsection();

             //Peierls creep
             prm.declare_entry ("Enable peierls","false",
                                Patterns::Bool (),
                                "Whether the model consider peierls creep. ");
             prm.enter_subsection ("Peierls creep");
             {
               prm.declare_entry ("E", "335e3",
                                  Patterns::Double (0),
                                  "Activation energy for peierls creep"
                                  "Unit: J/mol");
               prm.declare_entry ("V", "4.0e-6",
                                 Patterns::Double (0),
                                 "Activation volume for peierls creep"
                                 "Unit: m^3/mol");
               prm.declare_entry ("A", "1.92e-11",
                                  Patterns::Double (0),
                                  "Prefactor for peierls creep "
                                  "Yita=A^(-1/n)*(E+P*V)/(n*R*T)*Strain_rate_II^((1-n)/n)");
               prm.declare_entry ("n", "3.5",
                                  Patterns::Double (0),
                                  "Stress exponent for dislocation creep");
             }
             prm.leave_subsection();
          }
          prm.leave_subsection();
          prm.enter_subsection ("Data");
          {
            prm.declare_entry ("Solidus filename", "",
                               Patterns::Anything(),
                               "The solidus filename.");
            prm.declare_entry ("Liquidus filename", "",
                               Patterns::Anything(),
                               "The liquidus filename.");
            prm.declare_entry ("Latent heat","0",
                               Patterns::Double (),
                               "The latent hear of melt Units: J/kg");
          }
          prm.leave_subsection();

          prm.enter_subsection ("Melting Composition");
          {
            prm.declare_entry ("Cpx","-1",
                Patterns::Integer(),
                "The compositional field index of Cpx");
            prm.declare_entry ("H2O","-1",
                Patterns::Integer(),
                "The compositional field index of H2O");
            prm.declare_entry ("Default Cpx","0.15",
                Patterns::Double(),
                "The default Cpx concentration of Cpx if no compositoinal field is used");
          }
          prm.leave_subsection();
          prm.enter_subsection ("Melting Data");
          {
            prm.declare_entry ("A1","1085.7",
                Patterns::Double (),
                "A1");
            prm.declare_entry ("A2","132.9",
                Patterns::Double (),
                "A2");
            prm.declare_entry ("A3","-5.1",
                Patterns::Double (),
                "A3");
            prm.declare_entry ("B1","1475.0",
                Patterns::Double (),
                "B1");
            prm.declare_entry ("B2","80.0",
                Patterns::Double (),
                "B2");
            prm.declare_entry ("B3","-3.2",
                Patterns::Double (),
                "B3");
            prm.declare_entry ("C1","1780.0",
                Patterns::Double (),
                "C1");
            prm.declare_entry ("C2","45.0",
                Patterns::Double (),
                "C2");
            prm.declare_entry ("C3","-2.0",
                Patterns::Double (),
                "C3");
            prm.declare_entry ("r1","0.5",
                Patterns::Double (),
                "r1");
            prm.declare_entry ("r2","0.08",
                Patterns::Double (),
                "r2");
            prm.declare_entry ("beta1","1.50",
                Patterns::Double (),
                "beta1");
            prm.declare_entry ("beta2","1.50",
                Patterns::Double (),
                "beta2");
            prm.declare_entry ("K","43",
                Patterns::Double (),
                "K");
            prm.declare_entry ("gamma","0.75",
                Patterns::Double (),
                "gamma");
            prm.declare_entry ("D_H2O","0.01",
                Patterns::Double (),
                "D_H2O");
            prm.declare_entry ("X1","12.0",
                Patterns::Double (),
                "X1");
            prm.declare_entry ("X2","1.0",
                Patterns::Double (),
                "X2");        
            prm.declare_entry ("lambda","0.60",
                Patterns::Double (),
                "lambda");
            prm.declare_entry ("Cp","1000",
                Patterns::Double (),
                "Cp");
            prm.declare_entry ("deltaS","300",
                Patterns::Double (),
                "deltaS");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Melt<dim>::parse_parameters (ParameterHandler &prm)
    {
      Steinberger<dim>::parse_parameters(prm);
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt model");
        {
          prm.enter_subsection ("Viscosity");
          {
            viscosity_cutoff_low  = prm.get_double ("Viscosity cutoff low");
            viscosity_cutoff_high = prm.get_double ("Viscosity cutoff high");
            exponential_melt      = prm.get_double ("Exponential Melt");

            Assert (dynamic_cast<const GeometryModel::SphericalShell<dim>*>
                   (&(this->get_geometry_model()))
                    != 0,
                    ExcMessage ("Scaled melting production from 2D to 3D can only be used if the geometry "
                                "is a spherical shell."));
            depth_lower           = prm.get_double ("Lower mantle depth");
            density_difference    = Utilities::string_to_double(Utilities::split_string_list(prm.get("Density difference")));
            viscosity_difference    = Utilities::string_to_double(Utilities::split_string_list(prm.get("Viscosity factor")));
            
            //Yield stress
            is_yield_enable                     = prm.get_bool   ("Enable yield");
            prm.enter_subsection ("Yield stress");
            {
              yield_stress_surface             = Utilities::string_to_double(Utilities::split_string_list(prm.get("Yield stress surface")));
              yield_friction                   = Utilities::string_to_double(Utilities::split_string_list(prm.get("Yield friction")));
              yield_composition_factor         = Utilities::string_to_double(Utilities::split_string_list(prm.get("Composition factor")));
              yield_stress_max                 = Utilities::string_to_double(Utilities::split_string_list(prm.get("Yield stress max")));
            }
             prm.leave_subsection();


            //Diffusion creep
            is_diffusion_enable                = prm.get_bool   ("Enable diffusion");
            prm.enter_subsection ("Diffusion creep");
            {
              activation_energy_diffusion_um   = prm.get_double ("E_UM");
              activation_volume_diffusion_um   = prm.get_double ("V_UM");
              prefactor_diffusion_um           = prm.get_double ("A_UM");
            
              activation_energy_diffusion_lm   = prm.get_double ("E_LM");
              activation_volume_diffusion_lm   = prm.get_double ("V_LM");
              prefactor_diffusion_lm           = prm.get_double ("A_LM");
            }
            prm.leave_subsection();
            
            //Dislocation Creep (Upper Mantle)
            is_dislocation_enable            = prm.get_bool   ("Enable dislocation");
            prm.enter_subsection ("Dislocation creep");
            {
              activation_energy_dislocation    = prm.get_double ("E");
              activation_volume_dislocation    = prm.get_double ("V");
              prefactor_dislocation            = prm.get_double ("A");
              stress_exponent_dislocation      = prm.get_double ("n");
            }
            prm.leave_subsection();
            
            //Peierls Creep (Upper Mantle)
            is_peierls_enable                = prm.get_bool   ("Enable peierls");
            prm.enter_subsection ("Peierls creep");
            {
              activation_energy_peierls        = prm.get_double ("E");
              activation_volume_peierls        = prm.get_double ("V");
              prefactor_peierls                = prm.get_double ("A");
              stress_exponent_peierls          = prm.get_double ("n");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
          prm.enter_subsection ("Data");
          {
            solidus_filename=prm.get ("Solidus filename");
            liquidus_filename=prm.get ("Liquidus filename");
            Lh=prm.get_double ("Latent heat");
          }
          prm.leave_subsection();

          prm.enter_subsection ("Melting Composition");
          {
            i_composition_Cpx = prm.get_integer("Cpx");
            i_composition_H2O = prm.get_integer("H2O");
            default_Cpx       = prm.get_double ("Default Cpx");
          }
          prm.leave_subsection();
          prm.enter_subsection ("Melting Data");
          {
            melting_parameters.A1      = prm.get_double ("A1");
            melting_parameters.A2      = prm.get_double ("A2");
            melting_parameters.A3      = prm.get_double ("A3");
            melting_parameters.B1      = prm.get_double ("B1");
            melting_parameters.B2      = prm.get_double ("B2");
            melting_parameters.B3      = prm.get_double ("B3");
            melting_parameters.C1      = prm.get_double ("C1");
            melting_parameters.C2      = prm.get_double ("C2");
            melting_parameters.C3      = prm.get_double ("C3");
            melting_parameters.r1      = prm.get_double ("r1");
            melting_parameters.r2      = prm.get_double ("r2");
            melting_parameters.beta1   = prm.get_double ("beta1");
            melting_parameters.beta2   = prm.get_double ("beta2");
            melting_parameters.K       = prm.get_double ("K");
            melting_parameters.gamma   = prm.get_double ("gamma");
            melting_parameters.D_H2O   = prm.get_double ("D_H2O");
            melting_parameters.X1      = prm.get_double ("X1");
            melting_parameters.X2      = prm.get_double ("X2");
            melting_parameters.lambda  = prm.get_double ("lambda");
            melting_parameters.Cp      = prm.get_double ("Cp");
            melting_parameters.deltaS  = prm.get_double ("deltaS");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Melt,
                                   "melt",
                                   "A material model that reads tables of pressure and temperature "
                                   "dependent material coefficients from files. The default values for "
                                   "this model's runtime parameters use a material description taken "
                                   "from the paper \\textit{Complex phase distribution and seismic velocity "
                                   "structure of the transition zone: Convection model predictions "
                                   "for a magnesium-endmember olivine-pyroxene mantle} by Michael H.G. "
                                   "Jacobs and Arie P. van den Berg, Physics of the Earth and Planetary "
                                   "Interiors, Volume 186, Issues 1-2, May 2011, Pages 36--48. See "
                                   "\\url{http://www.sciencedirect.com/science/article/pii/S0031920111000422}."
                                   "As well as partial melting calculation.")
  }
}

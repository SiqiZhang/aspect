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

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
   namespace internal
    {
      class MaterialLookup
      {
        public:
          MaterialLookup(const std::string &filename,
                         const bool interpol)
          {

            /* Initializing variables */
            interpolation = interpol;
            delta_press=-1.0;
            min_press=-1.0;
            delta_temp=-1.0;
            min_temp=-1.0;
            numtemp=0;
            numpress=0;

            std::string temp;
            std::ifstream in(filename.c_str(), std::ios::in);
            AssertThrow (in,
                         ExcMessage (std::string("Couldn't open file <") + filename));

            getline(in, temp); // eat first line
            getline(in, temp); // eat next line
            getline(in, temp); // eat next line
            getline(in, temp); // eat next line

            in >> min_temp;
            getline(in, temp);
            in >> delta_temp;
            getline(in, temp);
            in >> numtemp;
            getline(in, temp);
            getline(in, temp);
            in >> min_press;
            min_press *= 1e5;  // conversion from [bar] to [Pa]
            getline(in, temp);
            in >> delta_press;
            delta_press *= 1e5; // conversion from [bar] to [Pa]
            getline(in, temp);
            in >> numpress;
            getline(in, temp);
            getline(in, temp);
            getline(in, temp);

            Assert(min_temp >= 0.0, ExcMessage("Read in of Material header failed (mintemp)."));
            Assert(delta_temp > 0, ExcMessage("Read in of Material header failed (delta_temp)."));
            Assert(numtemp > 0, ExcMessage("Read in of Material header failed (numtemp)."));
            Assert(min_press >= 0, ExcMessage("Read in of Material header failed (min_press)."));
            Assert(delta_press > 0, ExcMessage("Read in of Material header failed (delta_press)."));
            Assert(numpress > 0, ExcMessage("Read in of Material header failed (numpress)."));


            max_temp = min_temp + (numtemp-1) * delta_temp;
            max_press = min_press + (numtemp-1) * delta_press;

            density_values.reinit(numtemp,numpress);
            thermal_expansivity_values.reinit(numtemp,numpress);
            specific_heat_values.reinit(numtemp,numpress);
            vp_values.reinit(numtemp,numpress);
            vs_values.reinit(numtemp,numpress);
            enthalpy_values.reinit(numtemp,numpress);

            unsigned int i = 0;
            while (!in.eof())
              {
                double temp1,temp2;
                double rho,alpha,cp,vp,vs,h;
                in >> temp1 >> temp2;
                in >> rho;
                if (in.fail())
                  {
                    in.clear();
                    rho = density_values[(i-1)%numtemp][(i-1)/numtemp];
                  }
                in >> alpha;
                if (in.fail())
                  {
                    in.clear();
                    alpha = thermal_expansivity_values[(i-1)%numtemp][(i-1)/numtemp];
                  }
                in >> cp;
                if (in.fail())
                  {
                    in.clear();
                    cp = specific_heat_values[(i-1)%numtemp][(i-1)/numtemp];
                  }
                in >> vp;
                if (in.fail())
                  {
                    in.clear();
                    vp = vp_values[(i-1)%numtemp][(i-1)/numtemp];
                  }
                in >> vs;
                if (in.fail())
                  {
                    in.clear();
                    vs = vs_values[(i-1)%numtemp][(i-1)/numtemp];
                  }
                in >> h;
                if (in.fail())
                  {
                    in.clear();
                    h = enthalpy_values[(i-1)%numtemp][(i-1)/numtemp];
                  }

                getline(in, temp);
                if (in.eof())
                  break;

                density_values[i%numtemp][i/numtemp]=rho;
                thermal_expansivity_values[i%numtemp][i/numtemp]=alpha;
                specific_heat_values[i%numtemp][i/numtemp]=cp;
                vp_values[i%numtemp][i/numtemp]=vp;
                vs_values[i%numtemp][i/numtemp]=vs;
                enthalpy_values[i%numtemp][i/numtemp]=h;

                i++;
              }
            Assert(i==numtemp*numpress, ExcMessage("Material table size not consistent with header."));

          }

          double
          specific_heat(double temperature,
                        double pressure) const
          {
            return value(temperature,pressure,specific_heat_values,interpolation);
          }

          double
          density(double temperature,
                  double pressure) const
          {
            return value(temperature,pressure,density_values,interpolation);
          }

          double
          thermal_expansivity(const double temperature,
                              const double pressure) const
          {
            return value(temperature,pressure,thermal_expansivity_values,interpolation);
          }

          double
          seismic_Vp(const double temperature,
                     const double pressure) const
          {
            return value(temperature,pressure,vp_values,false);
          }

          double
          seismic_Vs(const double temperature,
                     const double pressure) const
          {
            return value(temperature,pressure,vs_values,false);
          }

          double
          dHdT (const double temperature,
                const double pressure) const
          {
            const double h = value(temperature,pressure,enthalpy_values,interpolation);
            const double dh = value(temperature+delta_temp,pressure,enthalpy_values,interpolation);
            return (dh - h) / delta_temp;
          }

          double
          dHdp (const double temperature,
                const double pressure) const
          {
            const double h = value(temperature,pressure,enthalpy_values,interpolation);
            const double dh = value(temperature,pressure+delta_press,enthalpy_values,interpolation);
            return (dh - h) / delta_press;
          }
          
		  double
          dRhodp (const double temperature,
                  const double pressure) const
          {
            const double rho = value(temperature,pressure,density_values,interpolation);
            const double drho = value(temperature,pressure+delta_press,density_values,interpolation);
            return (drho - rho) / delta_press;
          }		  

          double
          value (const double temperature,
                 const double pressure,
                 const dealii::Table<2,
                 double> &values,
                 bool interpol) const
          {
            const double nT = get_nT(temperature);
            const unsigned int inT = static_cast<unsigned int>(nT);

            const double np = get_np(pressure);
            const unsigned int inp = static_cast<unsigned int>(np);

            Assert(inT<values.n_rows(), ExcMessage("not in range"));
            Assert(inp<values.n_cols(), ExcMessage("not in range"));

            if (!interpol)
              return values[inT][inp];
            else
              {
                // compute the coordinates of this point in the
                // reference cell between the data points
                const double xi = nT-inT;
                const double eta = np-inp;

                Assert ((0 <= xi) && (xi <= 1), ExcInternalError());
                Assert ((0 <= eta) && (eta <= 1), ExcInternalError());

                // use these coordinates for a bilinear interpolation
                return ((1-xi)*(1-eta)*values[inT][inp] +
                        xi    *(1-eta)*values[inT+1][inp] +
                        (1-xi)*eta    *values[inT][inp+1] +
                        xi    *eta    *values[inT+1][inp+1]);
              }
          }



        private:


          double get_nT(double temperature) const
          {
            temperature=std::max(min_temp, temperature);
            temperature=std::min(temperature, max_temp-delta_temp);
            Assert(temperature>=min_temp, ExcMessage("not in range"));
            Assert(temperature<=max_temp, ExcMessage("not in range"));
            return (temperature-min_temp)/delta_temp;
          }

          double get_np(double pressure) const
          {
            pressure=std::max(min_press, pressure);
            pressure=std::min(pressure, max_press-delta_press);
            Assert(pressure>=min_press, ExcMessage("not in range"));
            Assert(pressure<=max_press, ExcMessage("not in range"));
            return (pressure-min_press)/delta_press;
          }


          dealii::Table<2,double> density_values;
          dealii::Table<2,double> thermal_expansivity_values;
          dealii::Table<2,double> specific_heat_values;
          dealii::Table<2,double> vp_values;
          dealii::Table<2,double> vs_values;
          dealii::Table<2,double> enthalpy_values;


          double delta_press;
          double min_press;
          double max_press;
          double delta_temp;
          double min_temp;
          double max_temp;
          unsigned int numtemp;
          unsigned int numpress;
          bool interpolation;
      };
	}

  
  
	template <int dim>
	void
	Melt<dim>::
	initialize()
	{
		n_material_data = material_file_names.size();
		for (unsigned i = 0; i < n_material_data; i++)
        material_lookup.push_back(std_cxx1x::shared_ptr<internal::MaterialLookup>
                                  (new internal::MaterialLookup(datadirectory+material_file_names[i],interpolation)));
	
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
      const double R=  8.341; //TODO gasconstant (well its constant....)
	  const double depth = (this->get_geometry_model()).depth(position);
	  double depth_nd = depth/(this->get_geometry_model()).maximal_depth();
	  double T_nd = (temperature-reference_T)/reference_dT;
	  //const double StressZ=yield_stress+yield_stress_increase*depth;
	  //const double StressZ=yield_stress+yield_stress_increase*pressure;

    const double strain_rate_II=std::sqrt(std::fabs(second_invariant(strain_rate)));

	  
	  if (!strcmp(viscosity_model.c_str(),"Exponential"))
       {
         viscosity = reference_eta * std::exp(-std::log(exponential_T)*T_nd+std::log(exponential_P)*depth_nd);
       }
      else if (!strcmp(viscosity_model.c_str(),"Diffusion"))
        {
          //viscosity = prefactor_diffusion*exp((activation_energy_diffusion+activation_volume_diffusion*pressure)/(R*temperature));
          if(depth<depth_lower)
          {
              viscosity = prefactor_diffusion_um*exp((activation_energy_diffusion_um+activation_volume_diffusion_um*pressure)/(R*temperature));
          }
          else
          {
              viscosity = prefactor_diffusion_lm*exp((activation_energy_diffusion_lm+activation_volume_diffusion_lm*pressure)/(R*temperature));
          }
        }
      else if (!strcmp(viscosity_model.c_str(),"Dislocation"))
        {
          viscosity = std::min(1e24,std::pow(prefactor_dislocation,-1e0/stress_exponent)*
                               std::pow(std::max(strain_rate.norm(),1e-17),(1e0-stress_exponent)/stress_exponent)*
                               exp((activation_energy_dislocation+activation_volume_dislocation*pressure)/(stress_exponent*R*temperature)));
        }
      else if (!strcmp(viscosity_model.c_str(),"Composite"))
        {
          double smooth_layer=0.1*depth_lower;
          // Upper mantle
          double viscosity_diffusion,viscosity_dislocation;
          viscosity_diffusion = prefactor_diffusion_um*exp((activation_energy_diffusion_um+activation_volume_diffusion_um*pressure)/(R*temperature));
          viscosity_dislocation = std::pow(prefactor_dislocation,1e0/stress_exponent)*
                                std::pow(std::max(strain_rate_II,1e-17),(1.0-stress_exponent)/stress_exponent)*
                                exp((activation_energy_dislocation+activation_volume_dislocation*pressure)/(stress_exponent*R*temperature));
          double viscosity_um = viscosity_dislocation * viscosity_diffusion / (viscosity_diffusion +  viscosity_dislocation);
          // Lower mantle
          double viscosity_lm = prefactor_diffusion_lm*exp((activation_energy_diffusion_lm+activation_volume_diffusion_lm*pressure)/(R*temperature));
          
          double smooth_viscosity=exp(log(viscosity_lm/viscosity_um)*(1.+erf((depth-depth_lower)/smooth_layer))/2.);
          viscosity = viscosity_um*smooth_viscosity;
        }
      else
        {
          viscosity = reference_eta;
        }

      // Apply layer dependent factors
      if(depth<depth_litho)
          viscosity *= viscosity_factor_litho;
      else if(depth<depth_trans)
          viscosity *= 1.e0;
      else if(depth<depth_lower)
          viscosity *= viscosity_factor_trans;
      else
      {
          double smooth_layer=0.5*(depth_lower-depth_trans);
          double smooth_viscosity=exp(log(viscosity_factor_lower/viscosity_factor_trans)*(1.+erf((depth-depth_lower)/smooth_layer))/2.);
          viscosity *= viscosity_factor_trans*smooth_viscosity;
      }

       // Apply cutoff
       viscosity = std::max(viscosity,viscosity_cutoff_low);
       viscosity = std::min(viscosity,viscosity_cutoff_high);

	     // Apply melt
       double Melt_fraction;
       if(composition.size()==3)
           Melt_fraction=composition[2]/100;
       else
       {
           double radius=sqrt(position.square());
           Melt_fraction=Data_Melt.Melting_fraction(temperature,pressure,radius,0.,0.);
       }
       viscosity*=std::exp(-std::log(exponential_melt)*Melt_fraction);

       // Apply compositional difference, it can go outside cutoff
       if(viscosity_difference.size()>0)
         for(unsigned i=0;i<composition.size();i++)
           viscosity*=pow(viscosity_difference[i],composition[i]);

       // Apply yield stress
       if(yield_factor.size()>0 && strain_rate_II>1e-20)
       {
         //double viscosity_yield=std::max(viscosity_cutoff_low,std::min(viscosity,StressZ/std::max(1e-20,strain_rate_II)/2.0));
         //double viscosity_yield=std::max(viscosity_cutoff_low,std::min(viscosity,StressZ/strain_rate_II/2.0));
         double viscosity_new=viscosity;
         double composition_rest=1.0;
         for(unsigned i=0;i<composition.size();i++)
         {
           double yield_stress_p=yield_stress[i+1]+yield_stress_increase[i+1]*pressure;
           double viscosity_yield=std::max(viscosity_cutoff_low,std::min(viscosity,yield_stress_p/strain_rate_II/2.0));
           double composition_i=std::max(0.,std::min(1.,composition[i]));
           viscosity_new*=pow(viscosity_yield/viscosity,yield_factor[i+1]*composition_i);
           composition_rest-=composition_i;
         }
         {
           double yield_stress_p=yield_stress[0]+yield_stress_increase[0]*pressure;
           double viscosity_yield=std::max(viscosity_cutoff_low,std::min(viscosity,yield_stress_p/strain_rate_II/2.0));
           composition_rest=std::max(0.,std::min(1.,composition_rest));
           viscosity_new*=pow(viscosity_yield/viscosity,yield_factor[0]*composition_rest);
         }
         viscosity=viscosity_new;
       }
       
       return viscosity;
    }

    template <int dim>
    double
    Melt<dim>::
    get_deltat (const Point<dim> &position) const
    {
		/*
      if (!(&this->get_adiabatic_conditions()))
        return 0.0;
      static const bool a = this->include_adiabatic_heating();
      return a ? 0.0 : (this->get_adiabatic_conditions().temperature(position)
                        - this->get_adiabatic_surface_temperature());*/
	  return 0.0;
    }
    template <int dim>
    double
    Melt<dim>::
    specific_heat (const double temperature,
                   const double pressure,
                   const std::vector<double> &compositional_fields,
                   const Point<dim> &position) const
    {
      Assert ((n_material_data <= compositional_fields.size()) || (n_material_data == 1),
              ExcMessage("There are more material files provided than compositional"
                         " Fields. This can not be intended."));
      double cp = 0.0;
      if (!latent_heat)
        {
          if (n_material_data == 1)
            cp = material_lookup[0]->specific_heat(temperature, pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                cp += compositional_fields[i] * material_lookup[i]->specific_heat(temperature,pressure);
            }
        }
      else
        {
          if (n_material_data == 1)
            cp = material_lookup[0]->dHdT(temperature+get_deltat(position),pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                cp += compositional_fields[i] * material_lookup[i]->dHdT(temperature,pressure);
              cp = std::max(std::min(cp,6000.0),500.0);
            }
        }
      return cp;
    }

	
	

    template <int dim>
    double
    Melt<dim>::
    viscosity_ratio (const double temperature,
                     const double pressure,
                     const std::vector<double> &, /*composition*/
                     const SymmetricTensor<2,dim> &strain_rate,
                     const Point<dim> &position) const
    {
      const double R=  8.3143; //TODO gasconstant (well its constant....)
      const double depth = (this->get_geometry_model()).depth(position);

      if (viscosity_model != "Composite")
        return 1;
      /*
      const double viscosity_diffusion = std::min(1e22,(1e0/prefactor_diffusion)*
                                                  std::exp((activation_energy_diffusion+
                                                            activation_volume_diffusion*pressure)/(R*temperature)));
      const double viscosity_dislocation = std::min(1e22,std::pow(prefactor_dislocation,-1e0/stress_exponent)*
                                                    std::pow(strain_rate.norm(),(1e0-stress_exponent)/
                                                             stress_exponent)*
                                                    std::exp((activation_energy_dislocation+
                                                              activation_volume_dislocation*pressure)/(stress_exponent*R*temperature)));
      */
      if(depth<depth_lower)
      {
        double viscosity_diffusion,viscosity_dislocation;
        viscosity_diffusion = prefactor_diffusion_um*exp((activation_energy_diffusion_um+activation_volume_diffusion_um*pressure)/(R*temperature));
        viscosity_dislocation = std::pow(prefactor_dislocation,1e0/stress_exponent)*
                                std::pow(strain_rate.norm(),(1.0-stress_exponent)/stress_exponent)*
                                exp((activation_energy_dislocation+activation_volume_dislocation*pressure)/(stress_exponent*R*temperature));
        return std::max(1e17,viscosity_dislocation)/std::max(1e17,viscosity_diffusion);
      }
      else
        return 1.;
    }



    template <int dim>
    double
    Melt<dim>::
    reference_viscosity () const
    {
      return reference_eta;
    }



    template <int dim>
    double
    Melt<dim>::
    reference_density () const
    {
      return reference_rho;
    }



    template <int dim>
    double
    Melt<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return reference_alpha;
    }

    template <int dim>
    double
    Melt<dim>::
    reference_cp () const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    Melt<dim>::
    reference_thermal_diffusivity () const
    {
      return k_value/(reference_rho*reference_specific_heat);
    }

    template <int dim>
    double
    Melt<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      // this model assumes that the thermal conductivity is in fact constant
      return k_value;
    }

    template <int dim>
    double
    Melt<dim>::
    density (const double temperature,
             const double pressure,
             const std::vector<double> &compositional_fields,
             const Point<dim> &position) const
    {
      Assert ((n_material_data <= compositional_fields.size()) || (n_material_data == 1),
              ExcMessage("There are more material files provided than compositional"
                         " Fields. This can not be intended."));
      double rho = 0.0;
      if (n_material_data == 1)
      {
        rho = material_lookup[0]->density(temperature,pressure);
        /*
        if(compositional_fields.size()==1)
          rho+=density_difference*compositional_fields[0];
        if(compositional_fields.size()==3)
            rho+=density_difference*compositional_fields[1];
            */
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
                
      }
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            rho += compositional_fields[i] * material_lookup[i]->density(temperature,pressure);
        }
      return rho;
    }

    template <int dim>
    double
    Melt<dim>::
    thermal_expansion_coefficient (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &position) const
    {
      Assert ((n_material_data <= compositional_fields.size()) || (n_material_data == 1),
              ExcMessage("There are more material files provided than compositional"
                         " Fields. This can not be intended."));
      double alpha = 0.0;
      if (!latent_heat)
        {
          if (n_material_data == 1)
            alpha = material_lookup[0]->thermal_expansivity(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                alpha += compositional_fields[i] * material_lookup[i]->thermal_expansivity(temperature,pressure);
            }
        }
      else
        {
          double dHdp = 0.0;
          if (n_material_data == 1)
            dHdp += material_lookup[0]->dHdp(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                dHdp += compositional_fields[i] * material_lookup[i]->dHdp(temperature,pressure);
            }
          alpha = (1 - density(temperature,pressure,compositional_fields,position) * dHdp) / temperature;
          alpha = std::max(std::min(alpha,1e-3),1e-5);
        }
      return alpha;
    }

    template <int dim>
    double
    Melt<dim>::
    seismic_Vp (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &position) const
    {
      Assert ((n_material_data <= compositional_fields.size()) || (n_material_data == 1),
              ExcMessage("There are more material files provided than compositional"
                         " Fields. This can not be intended."));
      double vp = 0.0;
      if (n_material_data == 1)
        vp += material_lookup[0]->seismic_Vp(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            vp += compositional_fields[i] * material_lookup[i]->seismic_Vp(temperature,pressure);
        }
      return vp;
    }



    template <int dim>
    double
    Melt<dim>::
    seismic_Vs (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &position) const
    {
      Assert ((n_material_data <= compositional_fields.size()) || (n_material_data == 1),
              ExcMessage("There are more material files provided than compositional"
                         " Fields. This can not be intended."));
      double vs = 0.0;
      if (n_material_data == 1)
        vs += material_lookup[0]->seismic_Vs(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            vs += compositional_fields[i] * material_lookup[i]->seismic_Vs(temperature,pressure);
        }
      return vs;
    }



    template <int dim>
    double
    Melt<dim>::
    compressibility (const double temperature,
                     const double pressure,
                     const std::vector<double> &compositional_fields,
                     const Point<dim> &position) const
    {
      //return 0.0;
      double dRhodp = 0.0;
      if (n_material_data == 1)
        dRhodp += material_lookup[0]->dRhodp(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            dRhodp += compositional_fields[i] * material_lookup[i]->dRhodp(temperature,pressure);
        }
      const double rho = density(temperature,pressure,compositional_fields,position);
      return (1/rho)*dRhodp;	  
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
    bool
    Melt<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return ((dependence & NonlinearDependence::pressure)
              ||
              (dependence & NonlinearDependence::temperature)
			  ||
			  (dependence & NonlinearDependence::compositional_fields));
    }



    template <int dim>
    bool
    Melt<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return ((dependence & NonlinearDependence::pressure)
              ||
              (dependence & NonlinearDependence::temperature));
    }


    template <int dim>
    bool
    Melt<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence) const
    {
      // this coefficient is in fact constant in this model
      return false;
    }



    template <int dim>
    bool
    Melt<dim>::
    is_compressible () const
    {
      //return false;
	  return model_is_compressible;
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
		new_compositional_fields.resize(compositional_fields.size());
		for(unsigned int i=0;i<compositional_fields.size();i++)
			new_compositional_fields[i]=compositional_fields[i];
		if(compositional_fields.size()==4)
		{
	    	//Depletion,  compsitional field[0]
			// If melting fraction>0.01, the change of depletion = melting fracton extracted
			if(compositional_fields[2]>1. )
				if(pressure<24e9)
          new_compositional_fields[0]=compositional_fields[0]+compositional_fields[2];
			//Water, compsitional field[1]
			// If melting fraction>0.01, Water will be fully extracted
			if(compositional_fields[2]>1.)
				  new_compositional_fields[1]=0.;
			//Melt compsitional field[2]
			//static aspect::melting::Melting_data Data_Melt(solidus_filename,liquidus_filename);
			new_compositional_fields[2]=Data_Melt.Melting_fraction(temperature,pressure,sqrt(position.square()),compositional_fields[1],compositional_fields[0]);
      if(pressure<1e9)new_compositional_fields[2]=0.;
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
				if(new_compositional_fields[i]>100.)
					new_compositional_fields[i]=100.;
			}
		}
	    return;
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
      //static aspect::melting::Melting_data Data_Melt(solidus_filename,liquidus_filename);
      double radius=sqrt(position.square());
      double T_solidus,T_liquidus,water,depletion;
      if(compositional_fields.size()==3)
      {
        water=compositional_fields[1];
        depletion=compositional_fields[0];
      }
      else
      {
        water=0.;
        depletion=0.;
      }
      T_solidus=Data_Melt.get_solidus(pressure,radius,water,depletion);
      T_liquidus=Data_Melt.get_liquidus(pressure,radius,water,depletion);
      if(T_solidus<temperature && temperature<T_liquidus)
        return (Lh/temperature/(T_liquidus-T_solidus));
      else
        return 0.0;
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
      if(i_composition_Cpx>0 && i_composition_Cpx<(int)compositional_fields.size())
        Mcpx=compositional_fields[i_composition_Cpx];
      else
        Mcpx=default_Cpx;
      if(i_composition_H2O>0 && i_composition_H2O<(int)compositional_fields.size())
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
    void
    Melt<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Gravity", "30",
                             Patterns::Double (0),
                             "The value of the gravity constant."
                             "Units: $m/s^2$.");
          prm.declare_entry ("Composition", "standard",
                             Patterns::Anything (),
                             "The Composition of the model. ");
          prm.declare_entry ("Path to model data", "data/material-model/table/",
                             Patterns::DirectoryName (),
                             "The path to the model data. ");
          prm.declare_entry ("ComputePhases", "false",
                             Patterns::Bool (),
                             "whether to compute phases. ");
          prm.declare_entry ("Compressible", "true",
                             Patterns::Bool (),
                             "whether the model is compressible. ");
          prm.enter_subsection ("Viscosity");
          {
            prm.declare_entry ("Viscosity Model", "Exponential",
                               Patterns::Anything (),
                               "Viscosity Model");
            prm.declare_entry ("Reference Viscosity", "5e24",
                               Patterns::Double (0),
                               "The value of the constant viscosity. Units: $kg/m/s$.");
			      prm.declare_entry ("Viscosity increase lower mantle", "1e0",
                               Patterns::Double (0),
                               "The Viscosity increase (jump) in the lower mantle.");
            prm.declare_entry ("Yield stress", "1.17e8",
                               Patterns::List(Patterns::Double (0)),
                               "Yield stress for different materials, the first one is the "
                               "default material. There should be n_compostion+1 values in the "
                               "list. (Pa)");			
             prm.declare_entry ("Yield stress increase", "0.1",
                               Patterns::List(Patterns::Double (0)),
                               "Yield stress increase with pressure for different materials, "
                               "the first one is for the default material. There should be "
                               "n_compostion+1 values in the list.");
             prm.declare_entry ("Viscosity cutoff low","1e19",
                               Patterns::Double (0),
                                "The lowest viscosity cut off Unit: Pa s");
             prm.declare_entry ("Viscosity cutoff high","1e24",
                                Patterns::Double (0),
                                "The largest viscosity cut off Unit: Pa s");
             prm.declare_entry ("Exponential Melt", "1e0",
					                      Patterns::Double (0),
                                "Multiplication factor or melting fraction exponent");
             prm.declare_entry ("Viscosity factor lithosphere","1",
                               Patterns::Double (0),
                               "The Viscosity factor of lithosphere");
             prm.declare_entry ("Viscosity factor transition zone","1",
                               Patterns::Double (0),
                               "The Viscosity factor of transition zone.");
             prm.declare_entry ("Viscosity factor lower mantle","1",
                               Patterns::Double (0),
                               "The Viscosity factor of lower mantle.");
             prm.declare_entry ("Lithosphere depth","100e3",
                               Patterns::Double (0),
                               "The depth of lithosphere bottom.");
             prm.declare_entry ("Transition zone depth","410e3",
                               Patterns::Double (0),
                               "The depth of transition zone top.");
             prm.declare_entry ("Lower mantle depth","660e3",
                               Patterns::Double (0),
                               "The depth of lower mantle top.");
             prm.declare_entry ("Density difference", "",
                                Patterns::List(Patterns::Double ()),
                                "Density difference of different composition (kg/m^3)");
             prm.declare_entry ("Viscosity factor", "",
                                Patterns::List(Patterns::Double ()),
                                "Viscosity difference of different composition (Pa s)");
             prm.declare_entry ("Yield stress factor", "0",
                                Patterns::List(Patterns::Double (0)),
                                "Yield stress factor of different composition."
                                "It contains list of n_compositional_field+1 factors,"
                                "the fist one is for default composition (0 no yield, 1 have yield)");             

             prm.enter_subsection ("Exponential");
             {
              prm.declare_entry ("Exponential T", "1e0",
                                 Patterns::Double (0),
                                 "Multiplication factor or Temperature exponent");
              prm.declare_entry ("Exponential P", "1e0",
                                 Patterns::Double (0),
                                 "Multiplication factor or Pressure exponent");
              prm.declare_entry ("Mantle Temperature Jump", "4000.0",
                                 Patterns::Double (0),
                                 "Mantle temperature jump (K)");
             /*
			 prm.declare_entry ("Viscosity factor lithosphere","1",
                               Patterns::Double (0),
                               "The Viscosity factor of lithosphere");			
             prm.declare_entry ("Viscosity factor transition zone","1",
                               Patterns::Double (0),
                               "The Viscosity factor of transition zone.");
             prm.declare_entry ("Viscosity factor lower mantle","1",
                               Patterns::Double (0),
                               "The Viscosity factor of lower mantle.");
             prm.declare_entry ("Lithosphere depth","100e3",
                               Patterns::Double (0),
                               "The depth of lithosphere bottom.");
             prm.declare_entry ("Transition zone depth","410e3",
                               Patterns::Double (0),
                               "The depth of transition zone top.");
             prm.declare_entry ("Lower mantle depth","660e3",
                               Patterns::Double (0),
                               "The depth of lower mantle top.");			 
			 prm.declare_entry ("Viscosity cutoff low","1e-2",
					           Patterns::Double (0),
							    "The lowest viscosity cut off in the times of reference");
			 prm.declare_entry ("Viscosity cutoff high","1e2",
					           Patterns::Double (0),
							   "The largest viscosity cut off in the times of reference");*/
			}
            prm.leave_subsection();
            prm.enter_subsection ("Diffusion");
            {
              prm.declare_entry ("Activation energy diffusion upper mantle", "335e3",
                                 Patterns::Double (0),
                                 "activation energy for diffusion creep");
              prm.declare_entry ("Activation volume diffusion upper mantle", "4.0e-6",
                                 Patterns::Double (0),
                                 "activation volume for diffusion creep");
              prm.declare_entry ("Prefactor diffusion upper mantle", "1.92e-11",
                                 Patterns::Double (0),
                                 "prefactor for diffusion creep "
                                 "(1e0/prefactor)*exp((activation\\_energy+activation\\_volume*pressure)/(R*temperature))");
              prm.declare_entry ("Activation energy diffusion lower mantle", "335e3",
                                 Patterns::Double (0),
                                 "activation energy for diffusion creep");
              prm.declare_entry ("Activation volume diffusion lower mantle", "4.0e-6",
                                 Patterns::Double (0),
                                 "activation volume for diffusion creep");
              prm.declare_entry ("Prefactor diffusion lower mantle", "1.92e-11",
                                 Patterns::Double (0),
                                 "prefactor for diffusion creep "
                                 "(1e0/prefactor)*exp((activation\\_energy+activation\\_volume*pressure)/(R*temperature))");			  
            }
            prm.leave_subsection();
            prm.enter_subsection ("Dislocation");
            {
              prm.declare_entry ("Activation energy dislocation", "335e3",
                                 Patterns::Double (0),
                                 "activation energy for dislocation creep");
              prm.declare_entry ("Activation volume dislocation", "4.0e-6",
                                 Patterns::Double (0),
                                 "activation volume for dislocation creep");
              prm.declare_entry ("Prefactor dislocation", "1.92e-11",
                                 Patterns::Double (0),
                                 "prefactor for dislocation creep "
                                 "prefactor*exp((activation\\_energy+activation\\_volume*pressure)/(R*temperature))");
              prm.declare_entry ("Stress exponent", "3.5",
                                 Patterns::Double (0),
                                 "stress exponent for dislocation creep");
            }
            prm.leave_subsection();
            prm.enter_subsection ("Composite");
            {
              prm.declare_entry ("Activation energy diffusion upper mantle", "335e3",
                                 Patterns::Double (0),
                                 "activation energy for diffusion creep");
              prm.declare_entry ("Activation volume diffusion upper mantle", "4.0e-6",
                                 Patterns::Double (0),
                                 "activation volume for diffusion creep");
              prm.declare_entry ("Prefactor diffusion upper mantle", "1.92e-11",
                                 Patterns::Double (0),
                                 "prefactor for diffusion creep "
                                 "(1e0/prefactor)*exp((activation\\_energy+activation\\_volume*pressure)/(R*temperature))");
              prm.declare_entry ("Activation energy diffusion lower mantle", "335e3",
                                 Patterns::Double (0),
                                 "activation energy for diffusion creep");
              prm.declare_entry ("Activation volume diffusion lower mantle", "4.0e-6",
                                 Patterns::Double (0),
                                 "activation volume for diffusion creep");
              prm.declare_entry ("Prefactor diffusion lower mantle", "1.92e-11",
                                 Patterns::Double (0),
                                 "prefactor for diffusion creep "
                                 "(1e0/prefactor)*exp((activation\\_energy+activation\\_volume*pressure)/(R*temperature))");
              prm.declare_entry ("Activation energy dislocation", "540e3",
                                 Patterns::Double (0),
                                 "activation energy for dislocation creep");
              prm.declare_entry ("Activation volume dislocation", "14.0e-6",
                                 Patterns::Double (0),
                                 "activation volume for dislocation creep");
              prm.declare_entry ("Prefactor dislocation", "2.42e-10",
                                 Patterns::Double (0),
                                 "prefactor for dislocation creep "
                                 "(1e0/prefactor)*exp((activation\\_energy+activation\\_volume*pressure)/(R*temperature))");
              prm.declare_entry ("Stress exponent", "3.5",
                                 Patterns::Double (0),
                                 "stress exponent for dislocation creep");
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
		prm.enter_subsection("Steinberger model");
		{
		  prm.declare_entry ("Data directory", "data/material-model/steinberger/",
							 Patterns::DirectoryName (),
							 "The path to the model data. ");
		  prm.declare_entry ("Material file names", "pyr-ringwood88.txt",
							 Patterns::List (Patterns::Anything()),
							 "The file names of the material data. "
							 "List with as many components as active"
							 "compositional fields (material data is assumed to"
							 "be in order with the ordering of the fields). ");
		  prm.declare_entry ("Bilinear interpolation", "true",
							 Patterns::Bool (),
							 "whether to use bilinear interpolation to compute "
							 "material properties (slower but more accurate).");
		  prm.declare_entry ("Latent heat", "false",
							 Patterns::Bool (),
							 "whether to include latent heat effects in the"
							 "calculation of thermal expansivity and specific heat."
							 "Following the approach of Nakagawa et al. 2009.");
		  prm.leave_subsection();
		}

     	}
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Melt<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt model");
        {
          reference_rho         = prm.get_double ("Reference density");
          reference_T           = prm.get_double ("Reference temperature");
          k_value                 = prm.get_double ("Thermal conductivity");
          reference_specific_heat = prm.get_double ("Reference specific heat");
          reference_alpha       = prm.get_double ("Thermal expansion coefficient");
          composition           = prm.get ("Composition");
          data_directory        = prm.get ("Path to model data") + "/" + composition +"/";
          compute_phases        = prm.get_bool ("ComputePhases");
          model_is_compressible = prm.get_bool ("Compressible");
          prm.enter_subsection ("Viscosity");
          {
            viscosity_model      = prm.get ("Viscosity Model");
            reference_eta       = prm.get_double ("Reference Viscosity");
            increase_lower_mantle   = prm.get_double ("Viscosity increase lower mantle");
            viscosity_cutoff_low  = prm.get_double ("Viscosity cutoff low");
            viscosity_cutoff_high = prm.get_double ("Viscosity cutoff high");
            //yield_stress          = prm.get_double ("Yield stress");
            //yield_stress_increase = prm.get_double ("Yield stress increase");
            exponential_melt      = prm.get_double ("Exponential Melt");

			Assert (dynamic_cast<const GeometryModel::SphericalShell<dim>*>
                   (&(this->get_geometry_model()))
                    != 0,
                    ExcMessage ("Scaled melting production from 2D to 3D can only be used if the geometry "
                                "is a spherical shell."));
            //earth_radius = (dynamic_cast<const GeometryModel::SphericalShell<dim>&> (this->get_geometry_model())).outer_radius();
			//mantle_thickness = earth_radius - dynamic_cast<const GeometryModel::SphericalShell<dim>&> (this->get_geometry_model()).inner_radius();
			//reference_dT = (this->get_boundary_temperature()).maximal_temperature()
			//	             - (this->get_boundary_temperature()).minimal_temperature();
			//std::cout<<earth_radius<<","<<mantle_thickness<<","<< reference_dT<<std::endl;

			viscosity_factor_litho= prm.get_double ("Viscosity factor lithosphere");
			viscosity_factor_trans= prm.get_double ("Viscosity factor transition zone");
			viscosity_factor_lower= prm.get_double ("Viscosity factor lower mantle");
			depth_litho           = prm.get_double ("Lithosphere depth");
			depth_trans           = prm.get_double ("Transition zone depth");
			depth_lower           = prm.get_double ("Lower mantle depth");
      density_difference    = Utilities::string_to_double(Utilities::split_string_list(prm.get("Density difference")));
      yield_stress          = Utilities::string_to_double(Utilities::split_string_list(prm.get("Yield stress")));
      yield_stress_increase = Utilities::string_to_double(Utilities::split_string_list(prm.get("Yield stress increase")));
      viscosity_difference    = Utilities::string_to_double(Utilities::split_string_list(prm.get("Viscosity factor")));
      yield_factor            = Utilities::string_to_double(Utilities::split_string_list(prm.get("Yield stress factor")));

			if (viscosity_model == "Exponential")
              {
                prm.enter_subsection ("Exponential");
                {
                  exponential_T         = prm.get_double ("Exponential T");
                  exponential_P         = prm.get_double ("Exponential P");
                  reference_dT          = prm.get_double ("Mantle Temperature Jump");
                }
                prm.leave_subsection();
              }
            if (viscosity_model == "Diffusion")
              {
                prm.enter_subsection ("Diffusion");
                {
                  activation_energy_diffusion_um   = prm.get_double ("Activation energy diffusion upper mantle");
                  activation_volume_diffusion_um   = prm.get_double ("Activation volume diffusion upper mantle");
                  prefactor_diffusion_um           = prm.get_double ("Prefactor diffusion upper mantle");
                  activation_energy_diffusion_lm   = prm.get_double ("Activation energy diffusion lower mantle");
                  activation_volume_diffusion_lm   = prm.get_double ("Activation volume diffusion lower mantle");
                  prefactor_diffusion_lm           = prm.get_double ("Prefactor diffusion lower mantle");				  
                }
                prm.leave_subsection();
              }
            if (viscosity_model == "Dislocation")
              {
                prm.enter_subsection ("Dislocation");
                {
                  activation_energy_dislocation = prm.get_double ("Activation energy dislocation");
                  activation_volume_dislocation = prm.get_double ("Activation volume dislocation");
                  prefactor_dislocation         = prm.get_double ("Prefactor dislocation");
                  stress_exponent                = prm.get_double ("Stress exponent");
                }
                prm.leave_subsection();
              }
            if (viscosity_model == "Composite")
              {
                prm.enter_subsection ("Composite");
                {
                  activation_energy_diffusion_um   = prm.get_double ("Activation energy diffusion upper mantle");
                  activation_volume_diffusion_um   = prm.get_double ("Activation volume diffusion upper mantle");
                  prefactor_diffusion_um           = prm.get_double ("Prefactor diffusion upper mantle");
                  activation_energy_diffusion_lm   = prm.get_double ("Activation energy diffusion lower mantle");
                  activation_volume_diffusion_lm   = prm.get_double ("Activation volume diffusion lower mantle");
                  prefactor_diffusion_lm           = prm.get_double ("Prefactor diffusion lower mantle");
                  activation_energy_dislocation    = prm.get_double ("Activation energy dislocation");
                  activation_volume_dislocation    = prm.get_double ("Activation volume dislocation");
                  prefactor_dislocation            = prm.get_double ("Prefactor dislocation");
                  stress_exponent                  = prm.get_double ("Stress exponent");
                }
                prm.leave_subsection();
              }
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
        prm.enter_subsection("Steinberger model");
        {
          datadirectory        = prm.get ("Data directory");
          material_file_names  = Utilities::split_string_list
                                 (prm.get ("Material file names"));
          interpolation        = prm.get_bool ("Bilinear interpolation");
          latent_heat          = prm.get_bool ("Latent heat");

          prm.leave_subsection();
		  
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
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

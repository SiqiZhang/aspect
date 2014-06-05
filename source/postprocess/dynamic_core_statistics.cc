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
/*  $Id: heat_flux_statistics.cc 1790 2013-06-16 22:14:04Z bangerth $  */


#include <aspect/postprocess/dynamic_core_statistics.h>
#include <aspect/simulator_access.h>
#include <aspect/boundary_temperature/dynamic_core.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {

    template <int dim>
    DynamicCoreStatistics<dim>::DynamicCoreStatistics()
    {
      set_CMB_heat_flux(0.0);
    }

    template <int dim>
    double
    DynamicCoreStatistics<dim>::get_CMB_heat_flux() const
    {
      return(CMB_heat_flux);
    }

    template <int dim>
    void
    DynamicCoreStatistics<dim>::set_CMB_heat_flux(double heat_flux)
    {
      CMB_heat_flux=heat_flux;
    }

    template <int dim>
    std::pair<std::string,std::string>
    DynamicCoreStatistics<dim>::execute (TableHandler &statistics)
    {
      // now add all of the computed heat fluxes to the statistics object
      // and create a single string that can be output to the screen
      std::ostringstream screen_text;
		  /*
	      double T_CMB=SimulatorAccess<dim>::Get_T_CMB(),
				 R_i=SimulatorAccess<dim>::Get_R_i();*/
      AssertThrow (dynamic_cast<const BoundaryTemperature::Dynamic_core<dim>*> (&(SimulatorAccess<dim>::get_boundary_temperature()))
			  !=0,
			  ExcMessage ("Dynamic core statistics has to be working with dynamic core boundary conditions."));
      const struct BoundaryTemperature::_Core_Data* core_data
        =(dynamic_cast<const BoundaryTemperature::Dynamic_core<dim>&> 
		  (SimulatorAccess<dim>::get_boundary_temperature())).get_core_data();
      // Calculate the heat flux at the top and bottom boundaries
	  // copied from heat_flux_statistics.cc
      // create a quadrature formula based on the temperature element alone.
      const QGauss<dim-1> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.temperature).degree+1);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula,
                                        update_gradients      | update_values |
                                        update_normal_vectors |
                                        update_q_points       | update_JxW_values);

      std::vector<Tensor<1,dim> > temperature_gradients (quadrature_formula.size());
      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      std::map<types::boundary_id, double> local_boundary_fluxes;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      typename MaterialModel::Interface<dim>::MaterialModelInputs in(fe_face_values.n_quadrature_points, this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out(fe_face_values.n_quadrature_points, this->n_compositional_fields());

      // for every surface face on which it makes sense to compute a
      // heat flux and that is owned by this processor,
      // integrate the normal heat flux given by the formula
      //   j =  - k * n . grad T
      //
      // for the spherical shell geometry, note that for the inner boundary,
      // the normal vector points *into* the core, i.e. we compute the flux
      // *out* of the mantle, not into it. we fix this when we add the local
      // contribution to the global flux
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->at_boundary(f))
              {
                fe_face_values.reinit (cell, f);
                fe_face_values[this->introspection().extractors.temperature].get_function_gradients (this->get_solution(),
                    temperature_gradients);
                fe_face_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                    in.temperature);
                fe_face_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                               in.pressure);
                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  fe_face_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                      composition_values[c]);

                in.position = fe_face_values.get_quadrature_points();

                // since we are not reading the viscosity and the viscosity
                // is the only coefficient that depends on the strain rate,
                // we need not compute the strain rate. set the corresponding
                // array to empty, to prevent accidental use and skip the
                // evaluation of the strain rate in evaluate().
                in.strain_rate.resize(0);

                for (unsigned int i=0; i<fe_face_values.n_quadrature_points; ++i)
                  {
                    for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                      in.composition[i][c] = composition_values[c][i];
                  }

                this->get_material_model().evaluate(in, out);


                double local_normal_flux = 0;
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    const double thermal_conductivity
                      = out.thermal_conductivities[q];

                    local_normal_flux
                    +=
                      -thermal_conductivity *
                      (temperature_gradients[q] *
                       fe_face_values.normal_vector(q)) *
                      fe_face_values.JxW(q);
                  }

                local_boundary_fluxes[cell->face(f)->boundary_indicator()]
                += local_normal_flux;
              }

      // now communicate to get the global values
      std::map<types::boundary_id, double> global_boundary_fluxes;
      {
        // first collect local values in the same order in which they are listed
        // in the set of boundary indicators
        const std::set<types::boundary_id>
        boundary_indicators
          = this->get_geometry_model().get_used_boundary_indicators ();
        std::vector<double> local_values;
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p)
          local_values.push_back (local_boundary_fluxes[*p]);

        // then collect contributions from all processors
        std::vector<double> global_values;
        Utilities::MPI::sum (local_values, this->get_mpi_communicator(), global_values);

        // and now take them apart into the global map again
        unsigned int index = 0;
        // Scale 2D shell to 3D
        double scale_factor_2D_shell=1.;
        {
          const GeometryModel::SphericalShell<2>* shperical_shell_geometry =
              dynamic_cast<const GeometryModel::SphericalShell<2>*> (&(this->get_geometry_model()));
          if(shperical_shell_geometry!=NULL)
            scale_factor_2D_shell=2. * (shperical_shell_geometry->R0);
        }
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p, ++index)
          global_boundary_fluxes[*p] = global_values[index]*scale_factor_2D_shell;
        set_CMB_heat_flux(-global_boundary_fluxes[0]);
      }

      // now add all of the computed heat fluxes to the statistics object
      // and create a single string that can be output to the screen
      unsigned int index = 0;
      for (std::map<types::boundary_id, double>::const_iterator
           p = global_boundary_fluxes.begin();
           p != global_boundary_fluxes.end(); ++p, ++index)
        {
          const std::string name = "Outward heat flux through boundary with indicator "
                                   + Utilities::int_to_string(p->first)
                                   + " (TW)";
          statistics.add_value (name, p->second/1e12);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.set_precision (name, 4);
          statistics.set_scientific (name, false);

          // finally have something for the screen
          screen_text.precision(3);
          screen_text << p->second/1e12 << " TW,";

	    }
		

          const std::string name1 = "CMB Temperature (K)";
          //statistics.add_value (name1, SimulatorAccess<dim>::postprocess_dynamic_core.Tcmb);
		  statistics.add_value (name1, core_data->Ti);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.set_precision (name1, 2);
          statistics.set_scientific (name1, false);
		  
		  const std::string name2 = "Inner core radius (km)";
		  //statistics.add_value (name2, SimulatorAccess<dim>::postprocess_dynamic_core.Ri*1e-3);
		  statistics.add_value (name2, core_data->Ri*1e-3);
		  // also make sure that the other columns filled by the this object
		  // all show up with sufficient accuracy and in scientific notation
		  statistics.set_precision (name2, 2);
		  statistics.set_scientific (name2, false);

		  const std::string name8 = "Light element concentration (%)";
		  //statistics.add_value (name8, SimulatorAccess<dim>::postprocess_dynamic_core.Xi*100);
		  statistics.add_value (name8, core_data->Xi*100);
		  statistics.set_precision (name8, 4);
		  statistics.set_scientific (name8, false);
		  
           const std::string name3 = "Es (MW/K)";
           //statistics.add_value (name3, SimulatorAccess<dim>::postprocess_dynamic_core.Es*1e-6);
		   statistics.add_value (name3, core_data->Es*core_data->dT_dt*1e-6);
           // also make sure that the other columns filled by the this object
           // all show up with sufficient accuracy and in scientific notation
           statistics.set_precision (name3, 2);
           statistics.set_scientific (name3, false);

           const std::string name4 = "Er (MW/K)";
           //statistics.add_value (name4, SimulatorAccess<dim>::postprocess_dynamic_core.Er*1e-6);
           statistics.add_value (name4, core_data->Er*1e-6);
		   // also make sure that the other columns filled by the this object
           // all show up with sufficient accuracy and in scientific notation
           statistics.set_precision (name4, 2);
           statistics.set_scientific (name4, false);

           const std::string name5 = "El (MW/K)";
           //statistics.add_value (name5, SimulatorAccess<dim>::postprocess_dynamic_core.El*1e-6);
		   statistics.add_value (name5, core_data->El*core_data->dR_dt*1e-6);
           // also make sure that the other columns filled by the this object
           // all show up with sufficient accuracy and in scientific notation
           statistics.set_precision (name5, 2);
           statistics.set_scientific (name5, false);

           const std::string name6 = "Eg (MW/K)";
           //statistics.add_value (name6, SimulatorAccess<dim>::postprocess_dynamic_core.Eg*1e-6);
           statistics.add_value (name6, core_data->Eg*core_data->dR_dt*1e-6);
		   // also make sure that the other columns filled by the this object
           // all show up with sufficient accuracy and in scientific notation
           statistics.set_precision (name6, 2);
           statistics.set_scientific (name6, false);

           const std::string name7 = "Ek (MW/K)";
           //statistics.add_value (name7, SimulatorAccess<dim>::postprocess_dynamic_core.Ek*1e-6);
           statistics.add_value (name7, core_data->Ek*1e-6);
		   // also make sure that the other columns filled by the this object
           // all show up with sufficient accuracy and in scientific notation
           statistics.set_precision (name7, 2);
           statistics.set_scientific (name7, false);



          // finally have something for the screen
          screen_text.precision(5);
		  screen_text <<core_data->Ti <<" K,"
			          <<core_data->Ri*1e-3<<" km";

      return std::pair<std::string, std::string> ("Core data (Q_CMB/Q_surface/Tc/Ri)",
                                                  screen_text.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(DynamicCoreStatistics,
                                  "dynamic core statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the dynamic core.")
  }
}

/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/postprocess/heat_flux_map.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/lac/vector.h>
#include <deal.II/dofs/dof_tools.h>
//#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/fe/fe_system.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
      template <int dim>
      std::vector<std::vector<std::pair<double, double> > >
      compute_heat_flux_through_boundary_faces (const SimulatorAccess<dim> &simulator_access)
      {
        std::vector<std::vector<std::pair<double, double> > > heat_flux_and_area(simulator_access.get_triangulation().n_active_cells(),
                                                                                 std::vector<std::pair<double, double> >(GeometryInfo<dim>::faces_per_cell,
                                                                                     std::pair<double,double>()));

        // create a quadrature for the boundary cell volume
        const QGauss<dim> quadrature_formula (simulator_access.get_fe().base_element(simulator_access.introspection().base_elements.temperature).degree+1);
        FEValues<dim> fe_values (simulator_access.get_mapping(),
                                 simulator_access.get_fe(),
                                 quadrature_formula,
                                 update_gradients      | update_values |
                                 update_q_points       | update_JxW_values);

        // create a quadrature for the boundary faces
        const QGauss<dim-1> quadrature_formula_face (simulator_access.get_fe().base_element(simulator_access.introspection().base_elements.temperature).degree+1);

        FE_DGQ<dim> heat_flow_fe(simulator_access.get_fe().base_element(simulator_access.introspection().base_elements.temperature).degree);
        FEFaceValues<dim> heat_flow_fe_face_values (simulator_access.get_mapping(),
                                                    heat_flow_fe,
                                                    quadrature_formula_face,
                                                    update_normal_vectors |
                                                    update_values |
                                                    update_JxW_values);
        DoFHandler<dim> heat_flow_dof_handler(simulator_access.get_triangulation());
        heat_flow_dof_handler.distribute_dofs (heat_flow_fe);

        std::vector< Vector<double> >
          local_heat_flow_values_at_qpoints(dim),
          local_heat_flow_fe_values (dim);

        std::vector< std::vector<double> >
          values_at_face_qpoints(dim);

        for(unsigned int i=0;i<dim;i++)
        {
          local_heat_flow_values_at_qpoints[i].reinit(quadrature_formula.size());
          local_heat_flow_fe_values[i].reinit(heat_flow_fe.dofs_per_cell);
          values_at_face_qpoints[i].resize(quadrature_formula_face.size());
        }

        typename MaterialModel::Interface<dim>::MaterialModelInputs in(fe_values.n_quadrature_points, simulator_access.n_compositional_fields());
        typename MaterialModel::Interface<dim>::MaterialModelOutputs out(fe_values.n_quadrature_points, simulator_access.n_compositional_fields());

        std::vector<Tensor<1,dim> > temperature_gradients (quadrature_formula.size());

        FullMatrix<double> qpoint_to_dof_matrix (heat_flow_fe.dofs_per_cell,
                                                 quadrature_formula.size());
        FETools::compute_projection_from_quadrature_points_matrix
                    (heat_flow_fe,
                     quadrature_formula, quadrature_formula,
                     qpoint_to_dof_matrix);
        // loop over all of the surface cells and evaluate the heat flux
        typename DoFHandler<dim>::active_cell_iterator
        cell = simulator_access.get_dof_handler().begin_active(),
        endc = simulator_access.get_dof_handler().end(),
        dg_cell = heat_flow_dof_handler.begin_active();;

        for (; cell!=endc; ++cell, ++dg_cell)
          if (cell->is_locally_owned() && cell->at_boundary())
            {
              fe_values.reinit (cell);
              // Set use_strain_rates to false since we don't need viscosity
              in.reinit(fe_values, cell, simulator_access.introspection(), simulator_access.get_solution(), true);
              simulator_access.get_material_model().evaluate(in, out);

              // If we use a reference density for the temperature equation, use it to compute the heat flux
              // otherwise use full density.
              if (simulator_access.get_parameters().formulation_temperature_equation ==
                  Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
                {
                  for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
                    out.densities[q] = simulator_access.get_adiabatic_conditions().density(in.position[q]);
                }

              MaterialModel::MaterialAveraging::average (simulator_access.get_parameters().material_averaging,
                                                         cell,
                                                         fe_values.get_quadrature(),
                                                         fe_values.get_mapping(),
                                                         out);

              // Get the temperature gradients from the solution.
              fe_values[simulator_access.introspection().extractors.temperature].get_function_gradients (simulator_access.get_solution(),
                  temperature_gradients);

              Tensor<1,dim> heat_flux;
              for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
              {
                heat_flux = (-out.thermal_conductivities[q] * temperature_gradients[q] +
                            in.velocity[q] * out.densities[q] * out.specific_heat[q] * in.temperature[q]);
                for(unsigned int i=0;i<dim;i++)
                  local_heat_flow_values_at_qpoints[i](q)=heat_flux[i];
              }
              for(unsigned int i=0;i<dim;i++)
              {
                qpoint_to_dof_matrix.vmult (local_heat_flow_fe_values[i],
                                            local_heat_flow_values_at_qpoints[i]);
              }

              for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                if (cell->at_boundary(f))
                  {
                    heat_flow_fe_face_values.reinit (dg_cell, f);

                    Tensor<1,dim> heat_flow_vector;
                    double face_area = 0;
                    double face_heat_flow = 0.;

                    for (unsigned int q=0; q<heat_flow_fe_face_values.n_quadrature_points; ++q)
                      {
                        for(unsigned int i=0;i<dim;i++)
                        {
                          heat_flow_vector[i] = 0;
                          for(unsigned int i_dof=0;i_dof<heat_flow_fe.dofs_per_cell;i_dof++)
                            heat_flow_vector[i] += heat_flow_fe_face_values.shape_value(i_dof,q) * 
                                                   local_heat_flow_fe_values[i](i_dof);
                        }
                        face_heat_flow += heat_flow_fe_face_values.normal_vector(q) *
                                          heat_flow_vector * 
                                          heat_flow_fe_face_values.JxW(q);

                        face_area += heat_flow_fe_face_values.JxW(q);
                      }

                    heat_flux_and_area[cell->active_cell_index()][f].first = face_heat_flow;
                    heat_flux_and_area[cell->active_cell_index()][f].second = face_area;
                  }
            }

        return heat_flux_and_area;
      }
    }

    template <int dim>
    std::pair<std::string,std::string>
    HeatFluxMap<dim>::execute (TableHandler &)
    {
      std::vector<std::vector<std::pair<double, double> > > heat_flux_and_area =
        internal::compute_heat_flux_through_boundary_faces (*this);

      // have a stream into which we write the data. the text stream is then
      // later sent to processor 0
      std::ostringstream output;
      std::vector<std::pair<Point<dim>,double> > stored_values;

      // loop over all of the surface cells and evaluate the heat flux
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->at_boundary(f) &&
                (this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top" ||
                 this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "bottom"))
              {
                // evaluate position of heat flow to write into output file
                const Point<dim> midpoint_at_surface = cell->face(f)->center();

                const double flux_density = heat_flux_and_area[cell->active_cell_index()][f].first /
                                            heat_flux_and_area[cell->active_cell_index()][f].second;

                // store final position and heat flow
                stored_values.push_back (std::make_pair(midpoint_at_surface, flux_density));
              }


      // Write the solution to an output stream
      for (unsigned int i=0; i<stored_values.size(); ++i)
        {
          output << stored_values[i].first
                 << ' '
                 << stored_values[i].second
                 << std::endl;
        }



      std::string filename = this->get_output_directory() +
                             "heat_flux." +
                             Utilities::int_to_string(this->get_timestep_number(), 5);
      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        filename.append("." + Utilities::int_to_string (this->get_nonlinear_iteration(), 4));

      const unsigned int max_data_length = Utilities::MPI::max (output.str().size()+1,
                                                                this->get_mpi_communicator());
      const unsigned int mpi_tag = 567;

      // on processor 0, collect all of the data the individual processors send
      // and concatenate them into one file
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::ofstream file (filename.c_str());

          file << "# "
               << ((dim==2)? "x y" : "x y z")
               << " heat flux" << std::endl;

          // first write out the data we have created locally
          file << output.str();

          std::string tmp;
          tmp.resize (max_data_length, '\0');

          // then loop through all of the other processors and collect
          // data, then write it to the file
          for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
            {
              MPI_Status status;
              // get the data. note that MPI says that an MPI_Recv may receive
              // less data than the length specified here. since we have already
              // determined the maximal message length, we use this feature here
              // rather than trying to find out the exact message length with
              // a call to MPI_Probe.
              MPI_Recv (&tmp[0], max_data_length, MPI_CHAR, p, mpi_tag,
                        this->get_mpi_communicator(), &status);

              // output the string. note that 'tmp' has length max_data_length,
              // but we only wrote a certain piece of it in the MPI_Recv, ended
              // by a \0 character. write only this part by outputting it as a
              // C string object, rather than as a std::string
              file << tmp.c_str();
            }
        }
      else
        // on other processors, send the data to processor zero. include the \0
        // character at the end of the string
        {
          MPI_Send (&output.str()[0], output.str().size()+1, MPI_CHAR, 0, mpi_tag,
                    this->get_mpi_communicator());
        }

      return std::pair<std::string,std::string>("Writing heat flux map:",
                                                filename);
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(HeatFluxMap,
                                  "heat flux map",
                                  "A postprocessor that computes the (conductive) heat flux "
                                  "density across each boundary. The heat density flux is computed in "
                                  "outward direction, i.e., from the domain to the outside, "
                                  "using the formula $-k \\nabla T \\cdot \\mathbf n$, where "
                                  "$k$ is the thermal conductivity as reported by the material "
                                  "model, $T$ is the temperature, and $\\mathbf n$ is the outward "
                                  "normal. Note that the quantity so computed does not include "
                                  "any energy transported across the boundary by material "
                                  "transport in cases where $\\mathbf u \\cdot \\mathbf n \\neq 0$. "
                                  "The integrated heat flux for each boundary can be obtained "
                                  "from the heat flux statistics postprocessor.")
  }
}

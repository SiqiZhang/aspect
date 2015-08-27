#include <aspect/postprocess/melt_statistics.h>
#include <aspect/simulator_access.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/material_model/melt.h>
#include <aspect/melting.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    MeltStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the temperature element alone.
      // be defensive about determining that what we think is the temperature
      // element is it in fact
      
      //Assert (this->n_compositional_fields()==3,ExecMessage("Number of compositional fields should be 3."));
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(2).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);
      //Melt data
      const MaterialModel::Melt<dim> *material_model=dynamic_cast<const MaterialModel::Melt<dim> *>(&this->get_material_model());
      AssertThrow(material_model != 0,ExcMessage("This postprocess can only be worked with Melt material model."));
      //std::cout<<"Pass melt data"<<std::endl;

      std::vector<double> melting_fractions(n_q_points),
                          temperature(n_q_points),
                          pressure(n_q_points),
                          water(n_q_points),
                          depletion(n_q_points);
      std::vector<std::vector<double> > composition_fields;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      double local_melting_integral = 0;

      // compute the integral quantities by quadrature
      // Compositional_field[2] is melting fraction
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[this->introspection().extractors.temperature].get_function_values (
              this->get_solution(),temperature);
          fe_values[this->introspection().extractors.pressure].get_function_values (
              this->get_solution(),pressure);
          for(unsigned i=0;i<this->n_compositional_fields();i++)
          {
            std::vector<double> composition_i(n_q_points);
            fe_values[this->introspection().extractors.compositional_fields[i]].get_function_values (
                this->get_solution(), composition_i);
            composition_fields.push_back(composition_i);
          }
          //std::cout<<"Composition set"<<std::endl;

          for(unsigned int q=0;q<n_q_points;q++)
          {
            std::vector<double> composition_q(this->n_compositional_fields());
            for(unsigned j=0;j<this->n_compositional_fields();j++)
              composition_q[j] = composition_fields[j][q];

            // Correct the temperature and pressure for incompressible case
            double corrected_temperature,
                   corrected_pressure;
            if(material_model->is_compressible())
            {
              corrected_temperature = temperature[q];
              corrected_pressure    = pressure[q];
            }
            else
            {
              corrected_temperature = temperature[q] + this->get_adiabatic_conditions().temperature(fe_values.quadrature_point(q));
              corrected_pressure    = this->get_adiabatic_conditions().pressure(fe_values.quadrature_point(q));
            }

            melting_fractions[q]=material_model->melt_fraction(corrected_temperature,corrected_pressure,composition_q,fe_values.quadrature_point(q));
          }
          //std::cout<<"Melt fraction set"<<std::endl;

          //fe_values[this->introspection().extractors.compositional_fields[2]].get_function_values (this->get_solution(),
          //                                                                             melting_fractions);
          for (unsigned int q=0; q<n_q_points; ++q)
          {
            //Melting only be extracted to surface when melting fraction > 1%
            if(melting_fractions[q]>0.01)
              local_melting_integral += melting_fractions[q]*fe_values.JxW(q);
          }
        }
      double global_melting_integral
        = Utilities::MPI::sum (local_melting_integral, this->get_mpi_communicator());
      global_melting_integral*=1.0e-9/this->get_timestep()*year_in_seconds;// Change melt production to km^3/year
      // Scale 2D melt production to 3D
      if(dim==2)
      {
        Assert (dynamic_cast<const GeometryModel::SphericalShell<dim>*>
            (this->get_geometry_model())
            != 0,
            ExcMessage ("Scaled melting production from 2D to 3D can only be used if the geometry "
                        "is a spherical shell."));
        const double R0 = dynamic_cast<const GeometryModel::SphericalShell<dim>&> (this->get_geometry_model()).inner_radius();
        const double R1 = dynamic_cast<const GeometryModel::SphericalShell<dim>&> (this->get_geometry_model()).outer_radius();
        global_melting_integral*=4./3.*(pow(R1,3)-pow(R0,3))/(pow(R1,2)-pow(R0,2));
      }

      statistics.add_value ("Melting production (km^3/year)",
          global_melting_integral);
      statistics.set_precision ("Melting production (km^3/year)", 6);
      statistics.set_scientific ("Melting production (km^3/year)", true);

      std::ostringstream output;
      output.precision(6);
      output << global_melting_integral<<" km^3/year";
      return std::pair<std::string, std::string> ("Melting production rate:",
          output.str());

    }
  }
}
// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MeltStatistics,
                                  "melt statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the melting production rate.")
  }
}

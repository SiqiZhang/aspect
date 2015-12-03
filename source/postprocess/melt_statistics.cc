#include <aspect/postprocess/melt_statistics.h>
#include <aspect/simulator_access.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/material_model/melt.h>
#include <aspect/melting.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <string>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    MeltStatistics<dim>::execute (TableHandler &statistics)
    {
      //Melt data
      const MaterialModel::Melt<dim> *material_model=dynamic_cast<const MaterialModel::Melt<dim> *>(&this->get_material_model());
      AssertThrow(material_model != 0,ExcMessage("This postprocess can only be worked with Melt material model."));
      
      //Only works with 2D spherical shell geometry.
      const GeometryModel::SphericalShell<2>* spherical_shell_geometry =
            dynamic_cast<const GeometryModel::SphericalShell<2>*> (&(this->get_geometry_model()));
      if(spherical_shell_geometry != NULL)
      {
        double R1 = spherical_shell_geometry->R1;
        double melt_grid_depth = material_model->get_extraction_depth();
        double R0 = std::max(spherical_shell_geometry->R0, R1 - melt_grid_depth);
        melt_grid.set_grid(R0,R1,melt_grid_nr,melt_grid_nh);
      }

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
              corrected_temperature = temperature[q] + this->get_adiabatic_conditions().temperature(fe_values.quadrature_point(q))
                                                     - this->get_adiabatic_surface_temperature();
              corrected_pressure    = this->get_adiabatic_conditions().pressure(fe_values.quadrature_point(q));
            }

            melting_fractions[q]=material_model->melt_extraction(corrected_temperature,corrected_pressure,composition_q,fe_values.quadrature_point(q));
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

          // Melting treatment
          // TODO: this only works in 2D
          if(melt_grid.is_grid_set())
            if(dim==2)
          {
            double cell_volume=0.;
            double cell_melt=0.;
            for (unsigned int q=0; q<n_q_points; ++q)
            {
              cell_volume += fe_values.JxW(q);
              //cell_melt   += 0.01*fe_values.JxW(q);
              cell_melt   += melting_fractions[q]*fe_values.JxW(q);
            }
            Tensor<1,dim> p0 = transfer_coord(cell->vertex(0));
            Tensor<1,dim> p1 = transfer_coord(cell->vertex(1));
            Tensor<1,dim> p2 = transfer_coord(cell->vertex(2));
            struct MeltGrid::MeltCell cell0;
            cell0.r0=std::min(p0[1],std::min(p1[1],p2[1]));
            cell0.r1=std::max(p0[1],std::max(p1[1],p2[1]));
            cell0.h0=std::min(p0[0],std::min(p1[0],p2[0]));
            cell0.h1=std::max(p0[0],std::max(p1[0],p2[0]));
            if((cell0.h1-cell0.h0)>M_PI)
            {
              std::swap(cell0.h0,cell0.h1);
              cell0.h1 += 2.*M_PI;
            }
            //std::cout<<"0 r="<<p0[1]<<", theta="<<p0[0]<<std::endl;
            //std::cout<<"1 r="<<p1[1]<<", theta="<<p1[0]<<std::endl;
            //std::cout<<"2 r="<<p2[1]<<", theta="<<p2[0]<<std::endl;
            melt_grid.add_melt(cell0.h0,cell0.h1,cell0.r0,cell0.r1,cell_melt/cell_volume);
          }

        }

      if(melt_grid.is_grid_set())
      {
        melt_grid.calculate_compression();
        if(melt_grid_output)
        {
          unsigned int time_step_num = this->get_timestep_number ();
          melt_grid.output_vtk(time_step_num);
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

    template <int dim>
    Tensor<1,dim>
    MeltStatistics<dim>::transfer_coord(const Point<dim> &p)
    {
      Tensor<1,dim> scoord;
      if(dim==2)
      {
        scoord[0] = std::atan2(p[0],p[1]);      //Theta
        if(scoord[0]<0.)
          scoord[0] += 2*M_PI;
        scoord[1] = std::sqrt(p.norm_square()); //R
      }
      else
      {
        //Not implimented
      }
      return scoord;
    }

    template<int dim>
    double
    MeltStatistics<dim>::get_compression(const Point<dim> &p) const
    {
      return melt_grid.get_compression(p);
    }


    template <int dim>
    MeltStatistics<dim>::MeltGrid::MeltGrid()
    {
      grid_set = false;
    }

    template <int dim>
    bool
    MeltStatistics<dim>::MeltGrid::is_grid_set() const
    {
      return grid_set;
    }


    template <int dim>
    void
    MeltStatistics<dim>::MeltGrid::set_grid(double R0, double R1, unsigned nr, unsigned nh)
    {
      this->R0=R0;
      this->R1=R1;
      this->nr=nr;
      this->nh=nh;

      dr=(R1-R0)/nr;
      dh=2.*M_PI/nh;

      melt_fraction.assign(nr*nh,0.);
      matrix_compression.assign((nr+1)*(nh+1),0.);

      grid_set = true;

    }

    template <int dim>
    void
    MeltStatistics<dim>::MeltGrid::add_melt(double h0, double h1, double r0, double r1, double f)
    {
      int ir_start = std::max<int>(0,std::floor((r0-R0)/dr));
      int ir_end   = std::min<int>(nr,std::ceil((r1-R0)/dr));
      if(h1<h0) h1 += 2.*M_PI;
      int ih_start = std::floor(h0/dh);
      int ih_end   = std::ceil(h1/dh);
      
      struct MeltCell cell_a;
      cell_a.h0=h0;
      cell_a.h1=h1;
      cell_a.r0=r0;
      cell_a.r1=r1;
      for(int i=ih_start;i<ih_end;i++)
      {
        struct MeltCell cell_b;
        cell_b.h0=dh*(i%nh);
        cell_b.h1=dh*(i%nh+1);
        for(int j=ir_start;j<ir_end;j++)
        {
          cell_b.r0=R0+dr*j;
          cell_b.r1=R0+dr*(j+1);
          /*
          std::cout<<"Adding cell["<<cell_a.h0<<"-"<<cell_a.h1<<"]["<<cell_a.r0<<"-"<<cell_a.r1<<
                     "] ["<<ih_start<<"-"<<ih_end<<"]["<<ir_start<<"-"<<ir_end<<"] "<<
                     "to cell [" <<cell_b.h0<<","<<cell_b.h1<<"]-["<<cell_b.r0<<","<<cell_b.r1<<
                     "] ["<<i<<"]["<<j<<"]" <<std::endl;
                     */
          struct MeltCell cell_c=get_overlapping_cell(cell_a,cell_b);
          melt_fraction[cell_index(i,j)]+=f*get_volume(cell_c)/get_volume(cell_b);
          //std::cout<<"Volume overlapping "<<get_volume(cell_c)<<", cell volume "<<get_volume(cell_b)<<std::endl;
        }
      }
    }

    template <int dim>
    struct MeltStatistics<dim>::MeltGrid::MeltCell
    MeltStatistics<dim>::MeltGrid::get_overlapping_cell(struct MeltCell a, struct MeltCell b) const
    {
      //a.h1 = modify_radian(a.h0,a.h1);
      //b.h0 = modify_radian(a.h0,b.h0);
      //b.h1 = modify_radian(a.h0,b.h1);
      struct MeltCell c;
      c.r0=std::max(a.r0,b.r0);
      c.r1=std::max(std::min(a.r1,b.r1),c.r0);
      c.h0=std::max(a.h0,b.h0);
      c.h1=std::max(std::min(a.h1,b.h1),c.h0);
      return c;
    }

    template <int dim>
    double
    MeltStatistics<dim>::MeltGrid::get_volume(struct MeltCell a) const
    {
      return .5*(a.h1-a.h0)*(a.r1*a.r1-a.r0*a.r0);
    }

    template <int dim>
    double
    MeltStatistics<dim>::MeltGrid::modify_radian(double ref_angle, double angle) const
    {
      double diff_angle=angle-ref_angle;
      if(diff_angle<0.)
        return angle+2.*M_PI;
      else if(diff_angle>=2.*M_PI)
        return angle-2.*M_PI;
      return angle;
    }

    template<int dim>
    void
    MeltStatistics<dim>::MeltGrid::calculate_compression() 
    {
      //matrix_compression.assign((nr+1)*(nh+1),0.);
      for(int i=0;i<=(int)nh;i++)
        for(int j=1;j<=(int)nr;j++)
        {
          double fraction=.5*(melt_fraction[cell_index(i%nh,(j-1))]+melt_fraction[cell_index(((i-1)%(int)nh+nh)%nh,(j-1))]);
          double dr0 = matrix_compression[point_index(i,j-1)];
          double r0  = (R0+dr*(j-1));
          double r1  = (R0+dr*j);
          double dr1 = r1 - sqrt((r0-dr0)*(r0-dr0)+(1.-fraction)*(r1*r1-r0*r0));
          matrix_compression[point_index(i,j)] = dr1;
          //if(i==0 || i==(int)nh)std::cout<<j<<" cell["<<i%nh<<","<<(j-1)<<"] ("<<cell_index(i%nh,(j-1)) <<")"
          //  <<" & ["<<((i-1)%(int)nh+nh)%nh<<","<<(j-1)<<"] ("<<cell_index(((i-1)%(int)nh+nh)%nh,(j-1)) <<")"<<std::endl;
        }
    }

    template<int dim>
    double
    MeltStatistics<dim>::MeltGrid::get_compression(const Point<dim> &p) const
    {
      Tensor<1,dim> p0 = MeltStatistics<dim>::transfer_coord(p);
      if(p0[1]>R0 && p0[1]<R1)
      {
        // Linear interpolation
        unsigned i0=p0[0]/dh;
        unsigned j0=(p0[1]-R0)/dr;
        double w1=fmod(p0[0],dh)/dh;
        double w2=fmod(p0[1],dr)/dr;
        double compression =
            matrix_compression[point_index(i0,j0)] * (1.-w1)*(1.-w2)
          + matrix_compression[point_index((i0+1)%nh,j0)] * w1*(1.-w2)
          + matrix_compression[point_index(i0,j0+1)] * (1.-w1)*w2
          + matrix_compression[point_index((i0+1)%nh,j0+1)] * w1*w2;
        return compression;
      }
      else
        return 0.;
    }

    template <int dim>
    void
    MeltStatistics<dim>::MeltGrid::output_vtk(unsigned int time_step_num = 0) const
    {
      //unsigned int time_step_num = this->get_timestep_number ();
      std::string str_time_step = static_cast<std::ostringstream*>( &(std::ostringstream() << time_step_num) )->str();
      std::string file_name = std::string("melt_gird.") + str_time_step + ".vtk";
      FILE *fp=fopen(file_name.c_str(),"w");
      fprintf(fp,"# vtk DataFile Version 2.0\n");
      fprintf(fp,"Melt grid\n");
      fprintf(fp,"ASCII\n");
      fprintf(fp,"DATASET STRUCTURED_GRID\n");
      fprintf(fp,"DIMENSIONS %u %u %u\n",nr+1,nh+1,1);
      //fprintf(fp,"ORIGIN %f %f %f\n",0.,0.,0.);
      //fprintf(fp,"SPACING %e %e %e\n",dr,dh,0.);
      fprintf(fp,"POINTS %u float\n",(nr+1)*(nh+1));
      for(unsigned i=0;i<=nh;i++)
        for(unsigned j=0;j<=nr;j++)
          fprintf(fp,"%e %e %e\n",(R0+dr*j)*sin(dh*i),(R0+dr*j)*cos(dh*i),0.);
      fprintf(fp,"CELL_DATA %u\n",nh*nr);
      fprintf(fp,"SCALARS melt_fraction float 1\n");
      fprintf(fp,"LOOKUP_TABLE default\n");
      for(unsigned i=0;i<melt_fraction.size();i++)
        fprintf(fp,"%e\n",melt_fraction[i]);
      fprintf(fp,"SCALARS matrix_compression float 1\n");
      fprintf(fp,"LOOKUP_TABLE default\n");
      for(unsigned i=0;i<nh;i++)
          for(unsigned j=0;j<nr;j++)
          {
            Point<dim> p;
            p[0]=(R0+dr*(j+.5))*sin(dh*(i+.5));
            p[1]=(R0+dr*(j+.5))*cos(dh*(i+.5));
            fprintf(fp,"%e\n",get_compression(p));
          }

      fprintf(fp,"POINT_DATA %u\n",(nh+1)*(nr+1));
      fprintf(fp,"SCALARS matrix_compression float 1\n");
      fprintf(fp,"LOOKUP_TABLE default\n");
      for(unsigned i=0;i<matrix_compression.size();i++)
        fprintf(fp,"%e\n",matrix_compression[i]);

      fclose(fp);
    }

    template <int dim>
    void
    MeltStatistics<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Melt statistics");
        {
          prm.declare_entry("Grid points horizontal","360",
                            Patterns::Integer(0),
                            "Seperate uniforme grid is used for melt treatment. "
                            "This defines the number of points used in horizontal.");
          prm.declare_entry("Grid points vertical","50",
                            Patterns::Integer(0),
                            "Seperate uniforme grid is used for melt treatment. "
                            "This defines the number of points used in vertical.");
          prm.declare_entry("Output grid", "false",
                            Patterns::Bool(),
                            "Seperate uniforme grid is used for melt treatment. "
                            "Set this to true will output the grid each timestep to check."); 
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    MeltStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Melt statistics");
        {
          melt_grid_nh      = prm.get_integer("Grid points horizontal");
          melt_grid_nr      = prm.get_integer("Grid points vertical");
          melt_grid_output  = prm.get_bool   ("Output grid");
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
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MeltStatistics,
                                  "melt statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the melting production rate.")
  }
}

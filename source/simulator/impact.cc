#include <aspect/impact.h>
#include <aspect/global.h>

#include <aspect/simulator.h>
#include <deal.II/base/point.h>
#include <deal.II/fe/fe_values.h>

#include <iostream>
#include <stdio.h>

//using namespace dealii;
namespace aspect
{
  namespace impact
  {
    template<int dim>
    ImpactFunction<dim>::ImpactFunction()
    :
      rotate()
    {
    }
    
    template<int dim>
    void
    ImpactFunction<dim>::Initialize(const std::string &filename,
                                    double time0,
                                    double R0,
                                    const Tensor<1,2> &surface_point_one,
                                    const Tensor<1,2> &surface_point_two,
                                    ConditionalOStream &pcout)
    {
      old_time=time0;
      //Read data & fill in vector Impacts_all;
      double lat,lon;
      const double D2P=acos(-1.)/180.;
      //const double D2P = number::PI/180.;
      struct Impact_Data Data1;
      char temp[256];
      this->R0=R0;
      std::ifstream in(filename.c_str(), std::ios::in);
      AssertThrow(in, ExcMessage(std::string("Couldn't open file <")+filename));
      in.getline(temp,256);
      //cout<<temp<<endl;
      while(!in.eof())
      {
        in>>Data1.time>>lon>>lat>>Data1.radius>>Data1.velocity;
        if(!in.fail())
        {
          Data1.time*=year_in_seconds*1e6; // Ma to seconds
          Data1.radius*=1e3;               // km to meters
          Data1.velocity*=1e3;             // km/s to m/s
          double Dc=0.305*Data1.radius*pow(Data1.velocity/1.e3,0.361),
                 r=R0-Dc;
          /*
             switch (dim){
             case 2:
             Data1.position(0)=r*sin(D2P*lon);
             Data1.position(1)=r*cos(D2P*lon);
             break;
             case 3:
             Data1.position(0)=r*cos(D2P*lat)*sin(D2P*lon);
             Data1.position(1)=r*cos(D2P*lat)*cos(D2P*lon);
             Data1.position(2)=r*sin(D2P*lat);
             break;
             default:
             AssertThrow(dim==2||dim==3, ExcMessage(std::string("Impact only works on 2D & 3D.")));
             break;
             }
             */
          Data1.position(0)=r*cos(D2P*lat)*cos(D2P*lon);
          Data1.position(1)=r*cos(D2P*lat)*sin(D2P*lon);
          Data1.position(2)=r*sin(D2P*lat);
          Impacts_all.push_back(Data1);
        }
        in.getline(temp,256);
        //cout<<temp<<endl;

       /* 
        if(Impacts_all.size()<10)
          std::cout<<Data1.time/year_in_seconds<<
                 ", "<<Data1.position(0)<<
                 ", "<<Data1.position(1)<<
                 ", "<<Data1.position(2)<<
                 ", "<<Data1.radius<<
                 ", "<<Data1.velocity<<std::endl;
                 */
      }
      {
        rotate.build_rotation_matrix(surface_point_one,surface_point_two);
        rotate.screen_output<dim>(surface_point_one,surface_point_two,pcout);
      }
      /*
      {
        FILE *fp = fopen("impact.vtk","w");
        fprintf(fp,"# vtk DataFile Version 2.0\n");
        fprintf(fp,"Impacts\n");
        fprintf(fp,"ASCII\n");
        fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
        fprintf(fp,"POINTS %d float\n",Impacts_all.size());
        for(unsigned i=0;i<Impacts_all.size();i++)
          fprintf(fp,"%e  %e  %e\n",Impacts_all[i].position(0),Impacts_all[i].position(1),Impacts_all[i].position(2));
        fclose(fp);
      }
      */
    }

    template<int dim>
    void
    ImpactFunction<dim>::update_time(double time)
    {
      //Update vector Impacts_active;
      Impacts_active.clear();
      for(unsigned i=0;i<Impacts_all.size();i++)
        if(old_time<=Impacts_all[i].time && Impacts_all[i].time<time)
          Impacts_active.push_back(Impacts_all[i]);
      //std::cout<<Impacts_active.size()<<" impact sources active at "<<time/year_in_seconds<<std::endl;
      //for(unsigned i=0;i<Impacts_active.size();i++)
      //  std::cout<<Impacts_active[i].position[0]<<","<<Impacts_active[i].position[1]<<","<<Impacts_active[i].position[2]<<std::endl;

      old_time=time;
    }    

    template<int dim>
    double
    ImpactFunction<dim>::deltaT(const Point<dim> &p, double pressure, double temperature)
    {
      double dT=0.;
      const double C=7.24e3,
                   S=1.25,
                   Cp=1.25e3,
                   Rho0=3.5e3;

      double Beta,f,Ps,Pd,Rc,r,Uc;
      for(unsigned i=0;i<Impacts_active.size();i++)
      {
        /*
           if(sqrt((Impacts_active[i].position-p).square())<Impacts_active[i].radius)
           dT=300.;*/
        AssertThrow(dim==2||dim==3, ExcMessage(std::string("Impact only works on 2D & 3D.")));
        Point<3> p1 = rotate.get_new_coord(p);

        Uc=Impacts_active[i].velocity/2.;
        r=sqrt((Impacts_active[i].position-p1).norm_square());
        Rc=0.451*Impacts_active[i].radius*pow(Impacts_active[i].velocity/1.e3,0.211);
        if(r<=Rc)
        {
          Ps=Rho0*(C+S*Uc)*Uc;
          //std::cout<<"["<<p<<"],["<<p1<<"] r="<<r<<std::endl;
        }
        else
	{
          //Ps=Rho0*(C+S*Uc)*Uc*pow(Rc/r,-1.84+2.61*log10(Impacts_active[i].velocity/1.e3));
	  Ps=Rho0*(C+S*Uc)*Uc*pow(Rc/r,-1.68+2.74*log10(Uc/1.e3));
          //Ps=Rho0*(C+S*Uc)*Uc*pow(Rc/r,1.25+0.625*log10(Impacts_active[i].velocity/1.e3));
	}
        Pd=Ps;
        if(Pd>0.)
        {
          //Beta=pow(C,2)*Rho0/(2.*S);
          //f=-Pd/Beta/(1-sqrt(2.*Pd/Beta+1.));
	  f= -2*S*Pd/(C*C*Rho0) * pow((1 - sqrt(4*S*Pd/(C*C*Rho0) + 1.0)),-1);
          //dT+=(Pd/(2.5*Rho0*S)*(1.-1./f)-pow(C/S,2)*(f-log(f)-1))/Cp;
	   /// Check ln vs log
          dT+=(Pd/(2.5*Rho0*S)*(1.-1./f)-pow(C/S,2)*(f-log(f)-1))/Cp;
        }
      }
      //if(dT==0)return(temperature);
      //std::cout<<"["<<p<<"] r="<<r<<", Pd="<<Pd<<",T0="<<temperature<<",dT="<<dT<<std::endl;
      double T_solid;
      if(pressure <12.e9)
        T_solid=1374.+130.e-9*pressure-5.6e-18*pow(pressure,2);
      else
        T_solid=1374.+130.e-9*12.e9-5.6e-18*pow(12.e9,2)+(pressure-12.e9)*15.e-9;
      if(temperature+dT>T_solid)
        return(std::max(T_solid,temperature));
      else
        return(temperature+dT);
    }

    namespace internal
    {
      void
      Rotate::build_rotation_matrix(const Tensor<1,2> &surface_point_one,
                                    const Tensor<1,2> &surface_point_two)
      {
        // get the Cartesian coordinates of the points the 2D model will lie in
        // this computation is done also for 3D since it is not expensive and the
        // template dim is currently not used here. Could be changed.
        const Tensor<1,3> point_one = cartesian_surface_coordinates(convert_tensor<2,3>(surface_point_one));
        const Tensor<1,3> point_two = cartesian_surface_coordinates(convert_tensor<2,3>(surface_point_two));
        const double normal[3] = {0.0,0.0,1.0};
        const Tensor<1,3> unrotated_normal_vector (normal);
        Tensor<1,3> rotated_normal_vector;
        cross_product(rotated_normal_vector,point_one,point_two);

        rotated_normal_vector /= rotated_normal_vector.norm();
        if ((rotated_normal_vector - unrotated_normal_vector).norm() > 1e-3)
        {
          Tensor<1,3> rotation_axis;
          cross_product(rotation_axis,unrotated_normal_vector,rotated_normal_vector);
          rotation_axis /= rotation_axis.norm();
          const double rotation_angle = std::acos(rotated_normal_vector*unrotated_normal_vector);
          rotation_matrix = rotation_matrix_from_axis(rotation_axis,rotation_angle);
          // Now apply the rotation that will project point_one onto the known point
          // (0,1,0).
          const Tensor<1,3> rotated_point_one = transpose(rotation_matrix) * point_one;

          const double point_one_coords[3] = {0.0,1.0,0.0};
          const Tensor<1,3> final_point_one (point_one_coords);

          const double second_rotation_angle = std::acos(rotated_point_one*final_point_one);
          Tensor<1,3> second_rotation_axis;
          cross_product(second_rotation_axis,final_point_one,rotated_point_one);
          second_rotation_axis /= second_rotation_axis.norm();

          const Tensor<2,3> second_rotation_matrix = rotation_matrix_from_axis(second_rotation_axis,second_rotation_angle);

          // The final rotation used for the model will be the combined
          // rotation of the two operation above. This is achieved by a
          // matrix multiplication of the rotation matrices.
          // This concatenation of rotations is the reason for using a
          // rotation matrix instead of a combined rotation_axis + angle
          rotation_matrix = rotation_matrix * second_rotation_matrix;
        }
        else
        {
          rotation_matrix[0][0] = 1.0;
          rotation_matrix[1][1] = 1.0;
          rotation_matrix[2][2] = 1.0;
        }
      }

      template<int dim>
      Point<3> Rotate::get_new_coord(const Point<dim> &p) const
      {
        Point<3> p1;
        Tensor<1,3> tensor_p_rotated;
        Tensor<1,dim> tensor_position;
        for(unsigned j=0;j<dim;j++)tensor_position[j]=p[j];
        if(dim==2)
        {
          tensor_p_rotated = rotation_matrix * convert_tensor<dim,3>(tensor_position);
        }
        else
        {
          tensor_p_rotated = convert_tensor<dim,3>(tensor_position);
        }
        for(unsigned j=0;j<3;j++)p1[j]=tensor_p_rotated[j];
        return p1;
      }

      template<int dim>
      void Rotate::screen_output(const Tensor<1,2> &surface_point_one,
                                 const Tensor<1,2> &surface_point_two,
                                 const ConditionalOStream &pcout) const
      {
        const Tensor<1,3> point_one = cartesian_surface_coordinates(convert_tensor<2,3>(surface_point_one));
        const Tensor<1,3> point_two = cartesian_surface_coordinates(convert_tensor<2,3>(surface_point_two));

        std::ostringstream output;

        output << std::setprecision (3) << std::setw(3) << std::fixed << std::endl
               << "   Setting up Impact."  << std::endl
               << std::endl;
        if (dim == 2)
        {
          Tensor<1,3> rotation_axis;
          const double rotation_angle = rotation_axis_from_matrix(rotation_axis,rotation_matrix);

          std_cxx1x::array<double,3> angles = angles_from_matrix(rotation_matrix);
          std_cxx1x::array<double,3> back_angles = angles_from_matrix(transpose(rotation_matrix));

          output << "   Input point 1 spherical coordinates: " << surface_point_one  << std::endl
                 << "   Input point 1 normalized cartesian coordinates: " << point_one  << std::endl
                 << "   Input point 1 rotated model coordinates: " << transpose(rotation_matrix) * point_one  << std::endl
                 << "   Input point 2 spherical coordinates: " << surface_point_two  << std::endl
                 << "   Input point 2 normalized cartesian coordinates: " << point_two  << std::endl
                 << "   Input point 2 rotated model coordinates: " << transpose(rotation_matrix) * point_two << std::endl
                 << std::endl <<  std::setprecision(2)
                 << "   Model will be rotated by " << -rotation_angle*180/numbers::PI
                 << " degrees around axis " << rotation_axis << std::endl
                 << "   The ParaView rotation angles are: " << angles[0] << " " << angles [1] << " " << angles[2] << std::endl
                 << "   The inverse ParaView rotation angles are: " << back_angles[0] << " " << back_angles [1] << " " << back_angles[2]

                 << std::endl;
        }

        pcout << output.str();
      }

      Tensor<1,3>
      Rotate::rotate_around_axis (const Tensor<1,3> &position,
                          const Tensor<1,3> &rotation_axis,
                          const double angle) const
      {
        Tensor<1,3> cross;
        cross_product(cross,rotation_axis,position);
        const Tensor<1,3> newpos = (1-std::cos(angle)) * rotation_axis*(rotation_axis*position) +
                                   std::cos(angle) * position + std::sin(angle) * cross;
        return newpos;

      }

      std_cxx1x::array<double,3>
      Rotate::angles_from_matrix (const Tensor<2,3> &rotation_matrix) const
      {
        std_cxx1x::array<double,3> orientation;
        // first rotate about y axis
        const double x2 = rotation_matrix[2][0];
        const double y2 = rotation_matrix[2][1];
        const double z2 = rotation_matrix[2][2];

        const double x3 = rotation_matrix[1][0];
        const double y3 = rotation_matrix[1][1];
        const double z3 = rotation_matrix[1][2];

        double d1 = sqrt(x2*x2 + z2*z2);

        double cosTheta, sinTheta;
        if (d1 < std::numeric_limits<double>::min())
        {
          cosTheta = 1.0;
          sinTheta = 0.0;
        }
        else
        {
          cosTheta = z2/d1;
          sinTheta = x2/d1;
        }

        double theta = atan2(sinTheta, cosTheta);
        orientation[1] = - theta * 180 / numbers::PI;

        // now rotate about x axis
        double d = sqrt(x2*x2 + y2*y2 + z2*z2);

        double sinPhi, cosPhi;
        if (d < std::numeric_limits<double>::min())
        {
          sinPhi = 0.0;
          cosPhi = 1.0;
        }
        else if (d1 < std::numeric_limits<double>::min())
        {
          sinPhi = y2/d;
          cosPhi = z2/d;
        }
        else
        {
          sinPhi = y2/d;
          cosPhi = (x2*x2 + z2*z2)/(d1*d);
        }

        double phi = atan2(sinPhi, cosPhi);
        orientation[0] = phi * 180 / numbers::PI;

        // finally, rotate about z
        double x3p = x3*cosTheta - z3*sinTheta;
        double y3p = - sinPhi*sinTheta*x3 + cosPhi*y3 - sinPhi*cosTheta*z3;
        double d2 = sqrt(x3p*x3p + y3p*y3p);

        double cosAlpha, sinAlpha;
        if (d2 < std::numeric_limits<double>::min())
        {
          cosAlpha = 1.0;
          sinAlpha = 0.0;
        }
        else
        {
          cosAlpha = y3p/d2;
          sinAlpha = x3p/d2;
        }

        double alpha = atan2(sinAlpha, cosAlpha);
        orientation[2] = alpha * 180 / numbers::PI;
        return orientation;
      }

      double
      Rotate::rotation_axis_from_matrix (Tensor<1,3> &rotation_axis,
                                         const Tensor<2,3> &rotation_matrix) const
      {
        double rotation_angle = std::acos(0.5 * (rotation_matrix[0][0] + rotation_matrix[1][1] + rotation_matrix[2][2] - 1));

        if (rotation_angle > std::numeric_limits<double>::min())
        {
          rotation_axis[0] = (rotation_matrix[2][1] - rotation_matrix[1][2]) / (2*std::sin(rotation_angle));
          rotation_axis[1] = (rotation_matrix[0][2] - rotation_matrix[2][0]) / (2*std::sin(rotation_angle));
          rotation_axis[2] = (rotation_matrix[1][0] - rotation_matrix[0][1]) / (2*std::sin(rotation_angle));
        }
        else
        {
          rotation_axis[0] = 0.0;
          rotation_axis[1] = 0.0;
          rotation_axis[2] = 1.0;
        }

        return rotation_angle;
      }

      Tensor<2,3>
      Rotate::rotation_matrix_from_axis (const Tensor<1,3> &rotation_axis,
                                         const double rotation_angle) const
      {
        Tensor<2,3> rotation_matrix;
        rotation_matrix[0][0] = (1-std::cos(rotation_angle)) * rotation_axis[0]*rotation_axis[0] + std::cos(rotation_angle);
        rotation_matrix[0][1] = (1-std::cos(rotation_angle)) * rotation_axis[0]*rotation_axis[1] - rotation_axis[2] * std::sin(rotation_angle);
        rotation_matrix[0][2] = (1-std::cos(rotation_angle)) * rotation_axis[0]*rotation_axis[2] + rotation_axis[1] * std::sin(rotation_angle);
        rotation_matrix[1][0] = (1-std::cos(rotation_angle)) * rotation_axis[1]*rotation_axis[0] + rotation_axis[2] * std::sin(rotation_angle);
        rotation_matrix[1][1] = (1-std::cos(rotation_angle)) * rotation_axis[1]*rotation_axis[1] + std::cos(rotation_angle);
        rotation_matrix[1][2] = (1-std::cos(rotation_angle)) * rotation_axis[1]*rotation_axis[2] - rotation_axis[0] * std::sin(rotation_angle);
        rotation_matrix[2][0] = (1-std::cos(rotation_angle)) * rotation_axis[2]*rotation_axis[0] - rotation_axis[1] * std::sin(rotation_angle);
        rotation_matrix[2][1] = (1-std::cos(rotation_angle)) * rotation_axis[2]*rotation_axis[1] + rotation_axis[0] * std::sin(rotation_angle);
        rotation_matrix[2][2] = (1-std::cos(rotation_angle)) * rotation_axis[2]*rotation_axis[2] + std::cos(rotation_angle);
        return rotation_matrix;
      }

      template <int in, int out>
      Tensor<1,out>
      Rotate::convert_tensor (const Tensor<1,in> &old_tensor) const
      {
        Tensor<1,out> new_tensor;
        for (unsigned int i = 0; i < out; i++)
          if (i < in) new_tensor[i] = old_tensor[i];
          else new_tensor[i] = 0.0;

        return new_tensor;
      }

      Tensor<1,3>
      Rotate::cartesian_surface_coordinates(const Tensor<1,3> &sposition) const
      {
        Tensor<1,3> ccoord;

        ccoord[0] = std::sin(sposition[0]) * std::cos(sposition[1]); // X
        ccoord[1] = std::sin(sposition[0]) * std::sin(sposition[1]); // Y
        ccoord[2] = std::cos(sposition[0]); // Z
        return ccoord;
      }
    }
  }

  template<int dim>
  void 
  Simulator<dim>::set_impacts()
  {
    //std::cout<<"Setting up impacts."<<std::endl;
    LinearAlgebra::BlockVector impacts_solution;
    // base element in the finite element is 2 for temperature 
    const unsigned int base_element = 2; 
    impacts_solution.reinit(system_rhs,false);
    // get the temperature support points
    const std::vector<Point<dim> > support_points
      = finite_element.base_element(base_element).get_unit_support_points();
    Assert (support_points.size() != 0,
        ExcInternalError());
    // create an FEValues object with just the temperature element
    FEValues<dim> fe_values (mapping, finite_element,
        support_points,
        update_values |
        update_quadrature_points | 
        update_JxW_values);
    std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
    for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
        cell != dof_handler.end(); ++cell)
      if (cell->is_locally_owned())
      {
        fe_values.reinit (cell);
        //fe_values_pressure.reinit (cell);
        // go through the temperature dofs and set modify their global values
        // according to the impacts
        cell->get_dof_indices (local_dof_indices);
        std::vector<double> pressure_values(finite_element.base_element(base_element).dofs_per_cell);
        fe_values[introspection.extractors.pressure].get_function_values(solution,pressure_values);
        /*
           for(unsigned int i=0;i<pressure_value.size();i++)
           cout<<" "<<pressure_value[i]<<" ";
           cout<<endl;*/
        for (unsigned int i=0; i<finite_element.base_element(base_element).dofs_per_cell; ++i)
        {
          const unsigned int system_local_dof
            = finite_element.component_to_system_index(/*temperature component=*/ dim+1,
                /*dof index within component=*/ i);
          double value = impacts.deltaT(fe_values.quadrature_point(i),pressure_values[i],solution(local_dof_indices[system_local_dof]));
          /*
             if(value>300.)
             cout<<"deltaT="<<value<<" at point ("<<fe_values.quadrature_point(i)(0)<<","
							 <<fe_values.quadrature_point(i)(1)<<")"<<endl;*/
          impacts_solution(local_dof_indices[system_local_dof])=value;
          //impacts_solution(local_dof_indices[system_local_dof])=pressure_values[i];
        }
      }
    impacts_solution.compress(VectorOperation::insert);
    // we should not have written at all into any of the blocks with
    // the exception of the current temperature block
    for (unsigned int b=0; b<impacts_solution.n_blocks(); ++b)
      if (b != 2)
        Assert (impacts_solution.block(b).l2_norm() == 0,
            ExcInternalError());
    current_constraints.distribute(impacts_solution);
    solution.block(2) = impacts_solution.block(2);
  }
  template class impact::ImpactFunction<2>;
  template class impact::ImpactFunction<3>;
  template class Simulator<2>;
  template class Simulator<3>;
};

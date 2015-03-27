#include <aspect/impact.h>
#include <aspect/global.h>

#include <aspect/simulator.h>
#include <deal.II/base/point.h>
#include <deal.II/fe/fe_values.h>

//using namespace dealii;
namespace aspect
{
  template<int dim>
  void
  ImpactFunction<dim>::Initialize(const std::string &filename, double time0, double R0)
  {
    old_time=time0;
    //Read data & fill in vector Impacts_all;
    double lat,lon;
    const double D2P=acos(-1.)/180.;
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
        Impacts_all.push_back(Data1);
      }
      in.getline(temp,256);
      //cout<<temp<<endl;

      /*
      std::cout<<Data1.time/year_in_seconds<<
				", "<<Data1.position(0)<<
				", "<<Data1.position(1)<<
				", "<<Data1.radius<<
				", "<<Data1.velocity<<std::endl;
      */
    }


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

    old_time=time;
  }

  template<int dim>
  double
  ImpactFunction<dim>::deltaT(const Point<dim> &p, double pressure, double temperature)
  {
    double dT=0.;
    const double C=7.24e3,
                 S=1.25,
                 Cp=1e3,
                 Rho0=3.5e3;

    double Beta,f,Ps,Pd,Rc,r,Uc;
    for(unsigned i=0;i<Impacts_active.size();i++)
    {
      /*
			if(sqrt((Impacts_active[i].position-p).square())<Impacts_active[i].radius)
				dT=300.;*/
      Uc=Impacts_active[i].velocity/2.;
      r=sqrt((Impacts_active[i].position-p).square());
      Rc=0.451*Impacts_active[i].radius*pow(Impacts_active[i].velocity/1.e3,0.211);
      if(r<=Rc)
        Ps=Rho0*(C+S*Uc)*Uc;
      else
        Ps=Rho0*(C+S*Uc)*Uc*pow(Rc/r,-1.84+2.61*log10(Impacts_active[i].velocity/1.e3));
      Pd=Ps-pressure;
      if(Pd>0.)
      {
        Beta=pow(C,2)*Rho0/(2.*S);
        f=-Pd/Beta/(1-sqrt(2.*Pd/Beta+1.));
        dT+=(Pd/(2.*Rho0*S)*(1.-1./f)-pow(C/S,2)*(f-log(f)-1))/Cp;
      }
    }
    if(dT==0)return(temperature);
    double T_solid=1374.+130.e-9*pressure-5.6e-18*pow(pressure,2);
    if(temperature+dT>T_solid && pressure<12.e9)
      return(std::max(T_solid,temperature));
    else
      return(temperature+dT);
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
    constraints.distribute(impacts_solution);
    solution.block(2) = impacts_solution.block(2);
  }
  template class ImpactFunction<2>;
  template class ImpactFunction<3>;
  template class Simulator<2>;
  template class Simulator<3>;
};

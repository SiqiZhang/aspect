#include <aspect/melting.h>
#include <aspect/simulator.h>
#include <deal.II/base/point.h>
#include <deal.II/fe/fe_values.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>

namespace aspect
{
	namespace melting
	{
		Melting_curve::Melting_curve(){}
		void Melting_curve::read(const std::string &filename)
		{
			std::ifstream in(filename.c_str(), std::ios::in);
			char temp[256];
			std::string T_Unit,P_Unit;
          	Num_points=0;
			if(in.fail())return;
			//std::cout<<"Start reading "<<filename<<std::endl;
			in.getline(temp,256);
			in>>T_Unit>>P_Unit;
			//std::cout<<T_Unit<<","<<P_Unit<<std::endl;
			in.getline(temp,256);
          	while(!in.eof())
          	{
            	double T,p;
                in>>T>>p;
                if(!in.fail())
                {
					//Unit switching
					if(T_Unit=="C")    T+=273.15; // Degree C to K
					if(P_Unit=="kbar") p*=1.e8;   // kbar to Pa
					if(P_Unit=="GPa")  p*=1.e9;   // GPa to Pa
					is_radius=false;              // Second column is pressure
					if(P_Unit=="km")
					{
						is_radius=true;           // Second column in radius instead of pressure
						p*=1.e3;                  // km to meters
					}
                	T_array.push_back(T);
                	P_array.push_back(p);
					//std::cout<<"T="<<T<<", P="<<p<<std::endl;
					Num_points++;
              	}
				in.getline(temp,256);
				//std::cout<<"TEMP="<<temp<<std::endl;
          	}
      	}
		
		double Melting_curve::T(const double p, const double radius) const
		{
			double T_value,P_value=is_radius?radius:p;
			if(T_array.size()==0)return(0.);
			//std::cout<<"Num_points="<<Num_points<<", size of array="<<T_array.size()<<std::endl;
			for(unsigned i=1;i<Num_points;i++)
			{
				//std::cout<<"Array="<<P_array[i]<<", Value="<<P_value<<std::endl;
				if(     (i==Num_points-1) || 
						(is_radius && P_value>P_array[i]) ||
						(!is_radius && P_value<P_array[i]) )
				{
					/*
					std::cout<<"P_value="
						     <<"T[i-1]="<<T_array[i-1]<<", "
						     <<"T[i]="  <<T_array[i]<<", "
							 <<"P[i-1]="<<P_array[i-1]<<", "
							 <<"P[i]="  <<P_array[i]<<std::endl;*/
					T_value=T_array[i-1]+(T_array[i]-T_array[i-1])/(P_array[i]-P_array[i-1])*(P_value-P_array[i-1]);
					break;
				}
			}
			//std::cout<<"P="<<p<<", Ts/Tl="<<T_value<<std::endl;
			return(T_value);
		}
		
		Melting_data::Melting_data()
		:
			Solidus(),
			Liquidus()
		{}

		void Melting_data::read(const std::string solidus_file, const std::string liquidus_file)
		{
			Solidus.read(solidus_file);
			Liquidus.read(liquidus_file);
		}

		double
		Melting_data::Melting_fraction(const double T, const double p, const double radius,
				                  const double water, const double depletion) const
		{
			double T_solidus=Solidus.T(p,radius),
				   T_liquidus=Liquidus.T(p,radius),
				   deltaT,
				   fraction=0.;

			//Depletion will increase solidus
			T_solidus+=depletion*(T_liquidus-T_solidus);

			// Olivine liquidus depression. (water in the units of wt%)
			// From MÃ©dard, E. and T. L. Grove (2007). "The effect of H2O on the olivine liquidus of basaltic melts: 
			// experiments and thermodynamic models." Contributions to Mineralogy and Petrology 155(4): 417-432.
			deltaT=-(40.4*water-2.97*pow(water,2)+0.0761*pow(water,3));
			T_solidus+=deltaT;
			T_liquidus+=deltaT;
			// Avoid solidus>liquidus
			if(T_solidus>T_liquidus)T_solidus=T_liquidus;
			//Depletion will reduce melt production capability
			if(T<=T_solidus) 
				fraction=0.;
			else if(T>T_liquidus) 
				fraction=100.;
			else 
				fraction=(T-T_solidus)/(T_liquidus-T_solidus)*(100.-depletion);
            /* 
			if(T_solidus<T && T<T_liquidus)
                 std::cout<<"T="<<T<<","
                     <<"P="<<p<<","
                     <<"R="<<radius<<","
                     <<"Water="<<water<<","
                     <<"depletion="<<depletion<<","
                     <<"Ts="<<T_solidus<<","
                     <<"Tl="<<T_liquidus<<","
					 <<"Melt_fraction="<<fraction<<std::endl;*/
			 return(fraction);
		}
		double
		Melting_data::get_solidus(const double p, const double radius,const double water, const double depletion) const
		{
			double T_solidus=Solidus.T(p,radius),
				   T_liquidus=Liquidus.T(p,radius),
				   deltaT;
			T_solidus+=depletion*(T_liquidus-T_solidus);
			deltaT=-(40.4*water-2.97*pow(water,2)+0.0761*pow(water,3));
    		T_solidus+=deltaT;
			T_liquidus+=deltaT;
			return(T_solidus);
		}

        double
        Melting_data::get_liquidus(const double p, const double radius,const double water, const double depletion) const
        {
            double T_solidus=Solidus.T(p,radius),
                   T_liquidus=Liquidus.T(p,radius),
                   deltaT;
            T_solidus+=depletion*(T_liquidus-T_solidus);
            deltaT=-(40.4*water-2.97*pow(water,2)+0.0761*pow(water,3));
            T_solidus+=deltaT;
            T_liquidus+=deltaT;
            return(T_liquidus);
        }		

  }

  template<int dim>
  void
  Simulator<dim>::revise_composition_melt()
  {
     LinearAlgebra::BlockVector melt_solution;
     const unsigned int base_element = 3;
     melt_solution.reinit(system_rhs,false);
     const std::vector<Point<dim> > support_points
             = finite_element.base_element(base_element).get_unit_support_points();
     Assert (support_points.size() != 0,
          ExcInternalError());
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
         cell->get_dof_indices (local_dof_indices);
         //Get presure value on composition nodes
         std::vector<double> pressure_values(finite_element.base_element(base_element).dofs_per_cell);
         fe_values[introspection.extractors.pressure].get_function_values(solution,pressure_values);
         //Get temperature value on composition nodes
         std::vector<double> temperature_values(finite_element.base_element(base_element).dofs_per_cell);
         fe_values[introspection.extractors.temperature].get_function_values(solution,temperature_values);
         for (unsigned int i=0; i<finite_element.base_element(base_element).dofs_per_cell; ++i)
         {
           unsigned int system_local_dof[parameters.n_compositional_fields];
           std::vector<double> composition_values(parameters.n_compositional_fields);
           std::vector<double> new_composition_value(parameters.n_compositional_fields);
           for(unsigned j=0;j<parameters.n_compositional_fields;j++)
           {
             system_local_dof[j]
               = finite_element.component_to_system_index(/*compositional component=*/ dim+2+j,
                   /*dof index within component=*/ i);
             composition_values[j]=solution[local_dof_indices[system_local_dof[j]]];
           }
           this->material_model.get()->revise_composition(temperature_values[i],pressure_values[i],
                                               composition_values,fe_values.quadrature_point(i),
                                               new_composition_value);
           for(unsigned j=0;j<parameters.n_compositional_fields;j++)
             melt_solution(local_dof_indices[system_local_dof[j]])=new_composition_value[j];
         }
       }
     melt_solution.compress(VectorOperation::insert);
     // we should not have written at all into any of the blocks with
     // the exception of the current compostion block
     for (unsigned int b=0; b<melt_solution.n_blocks(); ++b)
       if (b < 3)
         Assert (melt_solution.block(b).l2_norm() == 0,
             ExcInternalError());
     //constraints.distribute(melt_solution);
     for(unsigned j=0;j<parameters.n_compositional_fields;j++)
       solution.block(3+j) = melt_solution.block(3+j);
  }     
}
 // explicit instantiation of the functions we implement in this file
 namespace aspect
 {
 #define INSTANTIATE(dim) \
   template void Simulator<dim>::revise_composition_melt();

   ASPECT_INSTANTIATE(INSTANTIATE)
 }


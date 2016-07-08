/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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



#include <aspect/heating_model/impact.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/material_model/melt.h>

//TODO: Remove the dependency on the old impact code
#include <aspect/impact.h>

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    Impact<dim>::Impact ()
    {}

    template <int dim>
    void
    Impact<dim>::initialize ()
    {
      //Read the impactor list
      double lat,lon;
      const double D2P=acos(-1.)/180.;
      struct Impact_Data Data1;
      char temp[256];
      //TODO: setup R0
      const double R0=6371e3;
      std::ifstream in(filename.c_str(), std::ios::in);
      AssertThrow(in, ExcMessage(std::string("Couldn't open file <")+filename));
      in.getline(temp,256);
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
          Data1.position(0)=r*cos(D2P*lat)*cos(D2P*lon);
          Data1.position(1)=r*cos(D2P*lat)*sin(D2P*lon);
          Data1.position(2)=r*sin(D2P*lat);
          Impacts_all.push_back(Data1);
        }
        // Remove the reset of the line, so in the data file additional comments 
        // on each impact event is allowed.
        in.getline(temp,256);
        {
          rotate.build_rotation_matrix(surface_point_one,surface_point_two);
          // TODO pcout
          //rotate.screen_output<dim>(surface_point_one,surface_point_two,pcout);
        }
      }

    }

    template <int dim>
    void
    Impact<dim>::update()
    {
      const double time1 = this->get_time();
      const double time0 = time1 - this->get_timestep();
      // Update the active impactor list
      Impacts_active.clear();
      for(unsigned i=0;i<Impacts_all.size();i++)
        if(time0<=Impacts_all[i].time && Impacts_all[i].time<time1)
          Impacts_active.push_back(Impacts_all[i]);
      //--------------------------------
      //if(Impacts_active.size()>0)
      {
        const ConditionalOStream &pcout=this->get_pcout();
        pcout<<"Time: ["<<time0/year_in_seconds/1.e6<<" - "<<time1/year_in_seconds/1.e6<<"] Myr, "
             <<Impacts_active.size()<<" impactors is active."<<std::endl;
      }
      //--------------------------------
    }


    template <int dim>
    void
    Impact<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                           const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                           HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));
      const double timestep = this->get_timestep();
      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
      {
        if(timestep>0)
          heating_model_outputs.heating_source_terms[q] = deltaT(material_model_inputs.position[q],
                                                                 material_model_inputs.pressure[q],
                                                                 material_model_inputs.temperature[q],
                                                                 material_model_inputs.composition[q])
                                                          * material_model_outputs.densities[q]
                                                          * material_model_outputs.specific_heat[q]
                                                          / timestep;
        else
          heating_model_outputs.heating_source_terms[q] = 0;
        heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
      }
    }

    template <int dim>
    double
    Impact<dim>::
    deltaT (const Point<dim> &p, 
            double           pressure, 
            double           temperature,
            const std::vector<double> &compositional_fields) const
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
      if(dT < 0.0) {
        dT = 0.0;
      }
      //if(dT==0)return(temperature);
      //std::cout<<"["<<p<<"] r="<<r<<", Pd="<<Pd<<",T0="<<temperature<<",dT="<<dT<<std::endl;
      double T_solidus, T_liquidus;
      const MaterialModel::Melt<dim> *material_model=dynamic_cast<const MaterialModel::Melt<dim> *>(&(this->get_material_model()));
      AssertThrow(material_model != 0,ExcMessage("This impact module can only be worked with Melt material model."));

      if(material_model->Data_Melt.is_initialized())
      {
        const double water     = material_model->get_water(compositional_fields);
        const double depletion =  material_model->get_depletion(compositional_fields);
        T_solidus  = material_model->Data_Melt.get_solidus(pressure,sqrt(p.square()),water,depletion);
        T_liquidus = material_model->Data_Melt.get_liquidus(pressure,sqrt(p.square()),water,depletion);
      }
      else
      {
        if(pressure <12.e9)
          T_solidus=1374.+130.e-9*pressure-5.6e-18*pow(pressure,2);
        else
          T_solidus=1374.+130.e-9*12.e9-5.6e-18*pow(12.e9,2)+(pressure-12.e9)*15.e-9;
        T_liquidus = T_solidus;
      }
      double T_melt = T_solidus + (T_liquidus - T_solidus) * super_solidus_ratio;
      if(temperature+dT>T_melt)
        return(std::max(T_melt - temperature,0.));
      else
        return(dT);   
    }

    template <int dim>
    void
    Impact<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Impacts");
        {
          prm.declare_entry ("Data file","",
                             Patterns::Anything (),
                             "The data for impact events data, in the following format: "
                             "Time(Ma) Lon(degree) Lat(degree) Radius(km) Velocity(km/s)");
          prm.declare_entry ("Point one", "0.0,0.0",
                             Patterns::Anything (),
                             "Point that determines the plane in which a 2D model lies in. "
                             "Has to be in the format 'a,b' where a and b are latitude  "
                             " and longitude in degree.");
          prm.declare_entry ("Point two", "0.0,90.0",
                             Patterns::Anything (),
                             "Point that determines the plane in which a 2D model lies in. "
                             "Has to be in the format 'a,b' where a and b are latitude "
                             " and longitude in degree.");
          prm.declare_entry ("Super solidus ratio","0",
                             Patterns::Double(0,1),
                             "Large impact events will result strong melt. However the magma "
                             "pool cools very fast compare to mantle convection process. So "
                             "in the mantle convection simulation here, the temperature perturbation "
                             "introduced by impact events are limited by a temperature between mantle "
                             "solidus and liquidus. This parameter set this temperature as a temperature "
                             "ratio between mantle solidus and liquidus (eg. 0 means solidus, 1 means "
                             "liquidus). Any value greater than 0 will allow some partial melt left in the "
                             "impact zones.");
                             

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Impact<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Impacts");
        {
          filename=prm.get("Data file");
          aspect::Utilities::replace_path(filename);
          std::vector<double> point_one, point_two;
          point_one = dealii::Utilities::string_to_double
                      (dealii::Utilities::split_string_list(prm.get("Point one")));
          point_two = dealii::Utilities::string_to_double
                      (dealii::Utilities::split_string_list(prm.get("Point two")));
          AssertThrow(point_one.size()==2 && point_two.size()==2, ExcMessage("Wrong formate of the two points"));
          const double D2P=acos(-1.)/180.;
          surface_point_one[0] = (90. - point_one[0]) * D2P;
          surface_point_two[0] = (90. - point_two[0]) * D2P;
          surface_point_one[1] = point_one[1] * D2P;
          surface_point_two[1] = point_two[1] * D2P;
          if (dim == 2)
            Assert (surface_oint_one != surface_point_two,
                    ExcMessage ("To define a plane for the 2D model the two assigned points "
                                "may not be equal."));
          super_solidus_ratio = prm.get_double("Super solidus ratio");
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
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(Impact,
                                  "impact",
                                  "Implementation of a model in which account for the meteorite "
                                  "impact as a heating source.")
  }
}



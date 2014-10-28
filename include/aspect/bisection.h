#ifndef BISECTION_H
#define BISECTION_H

#include<string>
#include<cmath>

namespace Bisection
{
    template <class F>
    double bisecion_solve(F a,double min_value,double max_value,double max_error=1e-6,bool screen_out=false)
    {
        double x0=min_value,x2=max_value;
        unsigned int max_steps=log(max_error/(max_value-min_value))/log(0.5);
        double e0=a.error_fun(x0),e2=a.error_fun(x2);

        if(e0==0)return(x0);
        if(e2==0)return(x2);
        for(unsigned int i=0; i<max_steps; i++)
        {
            if(e0*e2>0) //No solution in the section
            {
                if(screen_out==true)
                    printf("Step: %3d  No solution found in section [% 8E,% 8E]\n",i,x0,x2);
                throw -1;
            }
            double x1=0.5*(x0+x2);
            double e1=a.error_fun(x1);
            if(screen_out==true)
            {
                printf("Step: %3d [% 6E,% 6E] Mid-point error= % 6E\n",i,x0,x2,e1);
            }
            if(e1==0)return(x1);
            if(e0*e1<0)
            {
                x2=x1;
                e2=e1;
            }
            else //  e1*e2<0
            {
                x0=x1;
                e0=e1;
            }
        }
        return 0.5*(x0+x2);
    };
}


#endif // BISECTION_H

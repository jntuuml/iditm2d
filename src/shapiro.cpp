#include "param.h"
#include "funcdef.h"

void shapiro(int n, double *z, int type)
{
    int i;
    double *y = new double[n];
    
    for (i = 0; i < n; i++) y[i]=z[i];

    if (type == 0) {
        for (i = 1; i < n-1; i++) z[i]=(1.0-cshap)*y[i]+0.5*cshap*(y[i-1]+y[i+1]);
    }
    else {
        for (i = 0; i < n; i++) z[i]=(1.0-cshap)*y[i]+0.5*cshap*(y[i-1]+y[i+1]);
    }

    delete[] y;

    return;
}

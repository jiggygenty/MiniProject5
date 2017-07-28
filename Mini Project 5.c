#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#define M 5
#define Q 10
double f(double x)
{
    return(exp(-Q*x));
}
void main()
{
    static double arraycons[M][M];
    static double alpha[M];
    static double g[M];
    static double u[M];
    int i,j,count;
    static double error[M];
    static double exfunction[M];
    static double epsilon[M];
    double lamda=1,x;
    double h=0.1;
    static double function[M];
    static double arrayres[M];
    double U0=0;
    double U1=0;
    double L=1;

    for(i=0;i<M;i++)
    {
        for(j=0;j<M;j++)
        {
            arraycons[i][j]=0;
        }
    }

    for(i=0;i<M;i++)
    {
        for(j=0;j<M;j++)
        {
            if(i==j)
            arraycons[i][j]=2.0+lamda*lamda*h*h;
            else if( i==j-1 || i==j+1)
            {
            arraycons[i][j]=-1;
            }
        }
    }

    for(i=0;i<M;i++)
    {
    function[i]=exp(-Q*f(x*i));
    }

    for(i=0;i<M;i++)
    {
    exfunction[i]=U0*sinh(lamda*(L-i*h))/sinh(lamda*L);
    exfunction[i]-=((f(i*h)/lamda)*((lamda*cosh(lamda*i*h)+Q*sinh(lamda*i*h)-lamda)/(lamda*lamda-Q*Q)));
    exfunction[i]+=(sinh(lamda*i*h)/sinh(lamda*L))*(U1+((f(L)/(lamda*lamda-Q*Q))))*(lamda*cosh(lamda*L)+Q*sinh(lamda*L)-lamda);
    }

    arrayres[0]=function[0]+U0;
    arrayres[M-1]=function[M-1]+U1;
    for(i=1;i<M-1;i++)
    {
    arrayres[i]=function[i];
    }

    alpha[0]=arraycons[0][0];
    for(i=1;i<M;i++)
    {
        alpha[i]=arraycons[i][i]-(arraycons[i][i-1]*arraycons[i-1][i])/alpha[i-1];
    }

    g[0]=arrayres[0];
    for(i=1;i<M;i++)
    {
        g[i]=arrayres[i]-((arraycons[i][i-1]*g[i-1])/alpha[i-1]);
    }


    u[M-1]=g[M-1]/alpha[M-1];
    for(i=M-2;i>=0;i--)
    {
        j=-1*(i-M)-1;
        u[i]=(g[i]-arraycons[i][i+1]*u[i+1])/alpha[i];
    }

    for(i=1;i<M-1;i++)
    {
        epsilon[i]=(arraycons[i][i-1]*u[i-1]+arraycons[i][i]*u[i]+arraycons[i-1][i]*u[i+1]-arrayres[i])/(fabs(arraycons[i][i-1]*u[i-1])+fabs(arraycons[i][i]*u[i])+fabs(arraycons[i-1][i]*u[i+1])+fabs(arrayres[i]));
        if (epsilon[i]>0.0001)
        printf("epsilon too large\n");
    }


    for(i=0;i<M;i++)
    {
        error[i]=fabs((exfunction[i]-u[i])/exfunction[i]);
    }

    printf("\tU Values:\tExact Values:\tEpsilon Values\n");
    count=0;
    for(i=0;i<M;i++)
    {
    count++;
        printf("%d\t%e\t%e\t%e\t%e\n",count,u[i],exfunction[i],epsilon[i],error[i]);
    }
}

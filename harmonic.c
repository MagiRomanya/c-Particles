#include <stdio.h>

void copy_array(double* array1, double array2[], int size)
{
	int i;
	for (i=0; i<size; i++)
	{
		array2[i] = *(array1 + i);
	}
}

double mitjana(double array[], int size)
{
	double avg = 0;
	int i;
	for (i=0; i<size; i++){
		avg += array[i];
	}
	return avg/size;
} 
//	EXEMPLE OSCIL·LADOR HARMONIC
//	
//	eq diferencial  dv_z/dt = -kz     (on v_z = dz/dt)
//		        dz/dt = v_z
//	
//	=> y = (v_z, z)	=> dy/dx = f(x,y) = (dv_z/dz, dz/dx) = ( -kz, v_z))
//	   x = t
//
// 		   Variables   Variable    Nombre de      Funció dy/dx = f(x,y)		            Pas		OUTPUT
// 		   dependents  independent variables						  Integració	|
// 		   (posicions,   (temps)   dependents 							|	|
// 		   velocitats)                                     					|	|
// 		         					   x       y[]  	  f(x,y)	|	|
void rungekutta4(const double y[],double x ,int size, int (*func)(double, const double*, double*), double h, double* y1)
{
	int i;
	double K1[size], K2[size], K3[size], K4[size], ycache[size];
	
	// Sets K1 = f(x,y)
	func(x,y, K1);

	// Sets K2 = f(x + h/2, y + h/2*K1)
	for (i=0; i<size; i++)
	{
		ycache[i] = y[i] + h/2*K1[i];
	}
	func(x + h/2, ycache, K2);

	// Sets K3 = f(x + h/2, y + h/2*K2)
	for (i=0; i<size; i++)
	{
		ycache[i] = y[i] + h/2*K2[i];
	}
	func(x + h/2, ycache, K3);

	// Sets K4 = f(x + h, y + h*K3)
	for (i=0; i<size; i++)
	{
		ycache[i] = y[i] + h*K3[i];
	}
	func(x + h, ycache, K4);

	// Computes the final answer
	for (i=0; i<size; i++)
	{
		y1[i] = y[i] + h/6*(K1[i] + 2*(K2[i]+K3[i]) + K4[i]);
	}
}

int harmonic(double x,const double y[], double f[])
{
	double k = 9.8; // Constant recuperadora de la molla
	f[0] = y[1];
	f[1] = -k*y[0];
	return 0;
}

int main()
{
	double y[2], y1[2], t, h;
	int i, iterations;
	t = 0;
	h = 0.01;
	iterations = 500;
	y[0] = 1;
	y[1] = 0;
	for (i=0; i < iterations; i++)
	{
		rungekutta4(y, t, 2, harmonic,h,y1);
		printf("%f, %f, %f\n",t, y[0], y[1]);
		copy_array(y1, y, 2);
		t+=h;
	}
	return 0;
}

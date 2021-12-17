#include <stdio.h>
#include <math.h>

void copy_array(double* array1, double array2[], int size)
{
	int i;
	for (i=0; i<size; i++)
	{
		array2[i] = *(array1 + i);
	}
}


float InvSqrt(float x)
{
        float xhalf = 0.5f * x;
        int i = *(int*)&x;            // store floating-point bits in integer
        i = 0x5f3759df - (i >> 1);    // initial guess for Newton's method
        x = *(float*)&i;              // convert new bits into float
        x = x*(1.5f - xhalf*x*x);     // One round of Newton's method
        return x;
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

// Calcula un pas de runge kutta 4
//
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


// Calcula la força d'un sol oscil·lador harmonic (serveix per testejar rungekutta)
int harmonic(double x,const double y[], double f[])
{
	double k = 9.8; // Constant recuperadora de la molla
	f[0] = y[1];
	f[1] = -k*y[0];
	return 0;
}


// Calcula la força de la gravetat per N cossos de la mateixa massa
// en DIM dimensions
int gravetat(double x, const double y[], double f[])
{
	const int N = 2; // Nombre de particules
	const int DIM = 2; // # Dimensions
	const int SIZE = N*DIM*2; // Tenim N particules amb DIM coordenades i DIM velocitats-> Serà el nombre de dimensions que ha de tenir y[] i f[]
	const int SIZE2 = N*DIM;
	int i, j, k;
	const double G = 1; // Constant que regula la força
	double sum2, dist;
	
	// Initialize f
	for (i=0; i<SIZE2; i++)
	{
		f[i] = y[SIZE2+i]; // dz/dx = v_z
		f[SIZE2+i] = 0;
	}
	for (i=0; i<SIZE2; i+=DIM)
	{
		for (j=0; j<i; j+=DIM)
		{
			sum2 = 0;
			for (k=0; k<DIM; k++)
			{
				sum2 +=(y[j+k]-y[i+k])*(y[j+k]-y[i+k]); // x**2 + y**2 + z**2
			}
			dist = InvSqrt(sum2);
			for (k=0; k<DIM; k++)
			{
				f[SIZE2+i+k] += G*(dist*dist*dist)*(y[j+k]-y[i+k]);
				f[SIZE2+j+k] += -G*(dist*dist*dist)*(y[j+k]-y[i+k]);
			}
		}
	}
	return 0;
}

int main()
{
	int size = 8;
	double y[size], y1[size], t, h;
	int i, iterations;
	t = 0;
	h = 0.1;
	iterations = 1000000;
	for (i=0; i < size; i++)
	{
		y[i] = 0;
	}
	y[3] = 300;
	y[6] = -0.01;
	y[4] = 0.01;
	for (i=0; i < iterations; i++)
	{
		rungekutta4(y, t, size, gravetat,h,y1);
		copy_array(y1, y, size);
		printf("%f, %f, %f, %f, %f, %f, %f, %f, %f\n",t, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);
		t+=h;
	}
	return 0;
}

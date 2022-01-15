#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <SDL2/SDL.h>

void copy_array(double* array1, double array2[], int size)
{
	int i;
	for (i=0; i<size; i++)
	{
		array2[i] = *(array1 + i);
	}
}

void print_array(const double v[], int len)
{
	int j;
	for (j=0; j < len-1; j++)
	{
		printf("%f, ", v[j]);
	}
	printf("%f\n", v[j]);

}

// Stolen Quake3 fast inverse square root function
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

// Calculates a step of RK4
//
//	EXAMPLE: Harmonic oscilator
//	
//	diferential eq  dv_z/dt = -kz     (where v_z = dz/dt)
//		        dz/dt = v_z
//	
//	=> y = (v_z, z)	=> dy/dx = f(x,y) = (dv_z/dz, dz/dx) = ( -kz, v_z))
//	   x = t
//
// 		   Dependent  Indepenedent    Number of      Function dy/dx = f(x,y)	          Integration	OUTPUT
// 		   variables  	variable      dependent						  	step	|
// 		   (positions,   (time)       variables							|	|
// 		   velocities)                                     					|	|
// 		         					   x       y[]  	  f(x,y)	|	|
void rungekutta4(const double y[], double x, int size, int (*func)(double, const double*, double*), double h, double* y1)
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


// Calculates a simple harmonic oscilator for testing RK4 porpuses
int harmonic(double x,const double y[], double f[])
{
	double k = 9.8; // Elasticity regulator constant
	f[0] = y[1];
	f[1] = -k*y[0];
	return 0;
}

#define DIM 2
#define N_PARTICLES 2

// Calculates the gravitational force for N bodies of the same mass in DIM dimensions
// Takes y[] as an argument y[] = {x1,...,xn,v1,...,vn}
// f[] is the output and has the form f[] = {v1,...,vn,a1,...,an}
// where v are the velocities and a the acceleration
int gravity(double x, const double y[], double f[])
{
	const int N = N_PARTICLES; // Particle number
	const int SIZE = N*DIM*2; // We have N particles, DIM coordinates and DIM velocities => number of dimentions of y[] f[]
	const int SIZE2 = N*DIM; // SIZE/2
	int i, j, k;
	const double G = 1.0f; // Gravitational constant multiplied by m
	double sum2, dist;
	
	// Initialize f
	for (i=0; i<SIZE2; i++)
	{
		f[i] = y[SIZE2+i]; // dz/dx = v_z
		f[SIZE2+i] = 0;
	}
	for (i=0; i<SIZE2 -DIM; i+=DIM)
	{
		for (j=i+DIM; j<SIZE2; j+=DIM)
		{
			sum2 = 0;
			// CALCULATE THE DISTANCE
			for (k=0; k<DIM; k++)
			{
				sum2 +=(y[j+k]-y[i+k])*(y[j+k]-y[i+k]); // x**2 + y**2 + z**2
			}
			dist = InvSqrt(sum2);
			// COMPUTES THE ACCELERATION
			for (k=0; k<DIM; k++)
			{
				double f_module = G*(dist*dist*dist)*(y[j+k]-y[i+k]);
				f[SIZE2+i+k] += f_module;
				f[SIZE2+j+k] += -f_module;
			}
		}
	}
	return 0;
}

// Draws the particles in the wondow in 2D
void DrawParticles(SDL_Renderer *renderer, double coord[])
{
	SDL_SetRenderDrawColor(renderer, 200, 200, 200, 0); // Particle color
	double x, y, radius;
	radius = 15;
	for (int i=0; i<N_PARTICLES; i++){
		for (int j=0; j<DIM*N_PARTICLES; j+=2){
			x = coord[j]*5 +500;	// Coordinates here we should do a transformation from coordiantes to pixels
			y = coord[j+1]*5 +500;
			SDL_Rect rect = {
			.x = x,
			.y = y,
			.w = radius,
			.h = radius,
			};
			SDL_RenderFillRect(renderer, &rect);
		}
	}
}

int main()
{
	int size = DIM * N_PARTICLES * 2;
	double y[size], y1[size], t, h;
	int i, j, iterations;
	double coord[] = {342.000000, 426.000000, 809.000000, 702.000000, 0.000000, 0.000000, 0.000000, 0.000000};
	// FINITE ELEMENTS VARIABLES
	h = 0.01;
	iterations = 0;
	
	// INITIAL CONDITIONS
	t = 0; 
	for (i=0; i < size; i++)
	{
		y[i] = 0;
	}
	y[2] = 10;
	y[5] = -10;
	y[6] = 0.2;
	y[8] = -0.2;
	y[11] = 0.2;
	for (i=0; i < size; i++)
	{
		y[i] = coord[i];
	}
	// RUNGE KUTTA LOOP
	for (i=0; i < iterations; i++)
	{
		rungekutta4(y, t, size, gravity,h,y1); // Runge Kutta step
		copy_array(y1, y, size); // Updates y to y1

		// Writes the result of each step in the terminal
		printf("%f, ", t);
		print_array(y, size);
		t+=h;
	}
	
	// GRAPHICAL INTERFACE WITH SDL
	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		fprintf(stderr, "ERROR: Could not initialize SDL: %s\n",SDL_GetError());
		exit(1);
	}
	SDL_Window *window =  SDL_CreateWindow("window test", 0, 0, 800, 600, SDL_WINDOW_RESIZABLE);
	if ( window == NULL ) {
		fprintf(stderr, "ERROR: Could not create a window %s\n", SDL_GetError());
		exit(1);
	}
	
	SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
	if ( renderer == NULL ){
		fprintf(stderr, "ERROR: Could not create a renderer %s\n", SDL_GetError());
		exit(1);
	}
	
	// Window Main loop
	bool quit=false;
	while (!quit) {
		// Event handeler
		SDL_Event event;
		while (SDL_PollEvent(&event)){
			switch(event.type){
			case SDL_QUIT:
				quit=true;
			}
		}
		SDL_SetRenderDrawColor(renderer, 48, 40, 60 ,0);
		SDL_RenderClear(renderer);
		
		SDL_SetRenderDrawColor(renderer, 255, 0, 0, 0);
		SDL_Rect rect = {
		.x = 50,
		.y = 50,
		.w = 50,
		.h = 50,
		};
		SDL_RenderFillRect(renderer, &rect);

		rungekutta4(y, t, size, gravity, h, y1); // Runge Kutta step
		copy_array(y1, y, size); // Updates y to y1
	
		DrawParticles(renderer, y);

		SDL_RenderPresent(renderer);
	}
	SDL_Quit();
	return 0;
}


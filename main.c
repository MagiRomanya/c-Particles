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
	if (len) printf("%f\n", v[j]);
	else printf("Print_array: Array with no size\n");

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
// 		         					   x       y[]  	  f(x,y) size	|	|
void rungekutta4(const double y[], double x, int size, int (*func)(double, const double*, double*, int), double h, double* y1)
{
	int i;
	double K1[size], K2[size], K3[size], K4[size], ycache[size];
	
	// Sets K1 = f(x,y)
	func(x,y,K1,size);

	// Sets K2 = f(x + h/2, y + h/2*K1)
	for (i=0; i<size; i++)
	{
		ycache[i] = y[i] + h/2*K1[i];
	}
	func(x + h/2, ycache, K2, size);

	// Sets K3 = f(x + h/2, y + h/2*K2)
	for (i=0; i<size; i++)
	{
		ycache[i] = y[i] + h/2*K2[i];
	}
	func(x + h/2, ycache, K3, size);

	// Sets K4 = f(x + h, y + h*K3)
	for (i=0; i<size; i++)
	{
		ycache[i] = y[i] + h*K3[i];
	}
	func(x + h, ycache, K4, size);

	// Computes the final answer
	for (i=0; i<size; i++)
	{
		y1[i] = y[i] + h/6*(K1[i] + 2*(K2[i]+K3[i]) + K4[i]);
	}
}


// Calculates a simple harmonic oscilator for testing RK4 porpuses
int harmonic(double x,const double y[], double f[], int size)
{
	double k = 9.8; // Elasticity regulator constant
	f[0] = y[1];
	f[1] = -k*y[0];
	return 0;
}

#define DIM 2
#define MAX_SIZE 8000
#define G_STR 800.0

// Calculates the gravitational force for N bodies of the same mass in DIM dimensions
// Takes y[] as an argument y[] = {x1,...,xn,v1,...,vn} (NULL terminated)
// f[] is the output and has the form f[] = {v1,...,vn,a1,...,an} 
// where v are the velocities and a the acceleration
int gravity(double x, const double y[], double f[], int size)
{
	int i, j, k;
	const int SIZE = size; 
	const int SIZE2 = size/2; // SIZE/2
	const int N=SIZE/(2*DIM); // We have N particles, DIM coordinates and DIM velocities => number of dimentions of y[] f[]
	const double G = G_STR; // Gravitational constant multiplied by m
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
void DrawParticles(SDL_Renderer *renderer, double coord[], int N)
{
	SDL_SetRenderDrawColor(renderer, 200, 200, 200, 0); // Particle color
	double x, y, radius, radius2;
	radius = 14;
	radius2 = 7; // half the radius 
	for (int j=0; j<N*DIM; j=j+2){
		x = coord[j];
		y = coord[j+1];
		SDL_Rect rect = {
		.x = x -radius2,
		.y = y -radius2,
		.w = radius,
		.h = radius,
		};
		SDL_RenderFillRect(renderer, &rect);
	}
}

int Add_Particle(double y[], const double coord[], int *size){
	double y_copy[*size];
	int size2 = (*size)/2; // half of the size
	if (*size + 4< MAX_SIZE){
		for (int i=0; i<*size; i++) y_copy[i]=y[i]; // Copy the array
		for (int i=0; i<(*size)/2; i+=2){
			y[i]=y_copy[i];
			y[i+1]=y_copy[i+1];
			y[size2+2+i]=y_copy[size2+i];
			y[size2+3+i]=y_copy[size2+i+1];
		}
		y[size2] = coord[0]; 
		y[size2+1] = coord[1]; 
		y[*size+2] = coord[2]; 
		y[(*size)+3] = coord[3]; 
		*size += 4;
		return 0;
		for (int i=0; i<4; i++){
			y[*size+i] = coord[i];
		}
		*size += 4;
		return 0;
	}
	else{
		printf("Exceeded maximum number of particles");
		return -1;
	}
}
#define get_size(x) DIM*x*2

int main()
{
	int N_particles = 0;
	int size = get_size(N_particles);
	double y[MAX_SIZE], y1[MAX_SIZE], t, h;
	double coord[4] = { 0., 0., 0., 0. };
	int i, j, iterations;
	
	// FINITE ELEMENTS VARIABLES
	h = 0.1;
	
	// INITIAL CONDITIONS
	t = 0; 
	for (i=0; i < size; i++)
	{
		y[i] = 0;
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
	int mouseX, mouseY;
	while (!quit) {
		// EVENT LOOP
		SDL_Event event;
		while (SDL_PollEvent(&event)){
			switch(event.type){
			case SDL_QUIT:
				quit=true;
			case SDL_MOUSEBUTTONDOWN:
				SDL_GetMouseState(&mouseX, &mouseY); // Where is the mouse?
				coord[0] = mouseX;
				coord[1] = mouseY;
				Add_Particle(y, coord, &size);
				N_particles++;
				break;
			}
		}
		SDL_Delay(17); // Waits to display frames
		// BACKGROUND COLOR
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
		
		// TRAJECTORY CALCULATIONS
		rungekutta4(y, t, size, gravity, h, y1); // Runge Kutta step
		copy_array(y1, y, size); // Updates y to y1

		// DISPLAY THE PARTICLES
		DrawParticles(renderer, y, N_particles);

		SDL_RenderPresent(renderer);
	}
	SDL_Quit();
	return 0;
}


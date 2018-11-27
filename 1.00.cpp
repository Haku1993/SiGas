
#include <iostream>
#include <cstdlib>
#include <fstream>

// global constants
const int N = 1;
const double G = 9.81; // m/s^2
const double DT = 0.01; // m/s^2
const double K = 300.234; // N/m

// body class
struct Body {
  double Rxold = 0, Ryold = 0, Rx = 0, Ry = 0, Vx = 0, Vy = 0, Fx = 0, Fy = 0;
  double mass = 0;
  double rad = 0;
  void arranque(double dt);
  void timestep(double dt);
};

void Body::arranque(double dt)
{
  Rxold = Rx - dt*Vx + Fx*dt*dt/(2*mass);   
  Ryold = Ry - dt*Vy + Fy*dt*dt/(2*mass);   
}

void Body::timestep(double dt)
{
  double tmp;

  tmp = Rx;
  Rx = 2*Rx - Rxold + Fx*dt*dt/(mass);   
  Vx = (Rx - Rxold)/(2*dt);
  Rxold = tmp;

  tmp = Ry;
  Ry = 2*Ry - Ryold + Fy*dt*dt/(mass);   
  Vy = (Ry - Ryold)/(2*dt);
  Ryold = tmp;
}

// function declarations

void set_masses(Body bodies[]);
 
void compute_forces(Body bodies[]);
void start(Body bodies[], double dt);
void evolve(Body bodies[], double dt);

void init_gnuplot(void);
void print_to_gnuplot(Body bodies[]);

// main
int main(void)
{
  srand(0); //random seed

  std::ofstream fout("datos.txt");

  Body bodies[N];

  bodies[0].rad = 0.45;
  bodies[0].Rx = 5;
  bodies[0].Ry = 4.545;
  bodies[0].Vx = 0.4; // you can modify the initial speed
  
  set_masses(bodies);
  compute_forces(bodies);

  init_gnuplot();

  start(bodies, DT);

  for (int it = 0; it < 5000; ++it) {
    fout << DT*it << " , " << bodies[0].Ry << std::endl; 
    compute_forces(bodies);
    evolve(bodies, DT);
    print_to_gnuplot(bodies);
  }

  fout.close();

  return 0;
}

void set_masses(Body bodies[])
{
  int ii;
  for (ii = 0; ii < N; ++ii) {
    bodies[ii].mass = 1 + double(rand())/RAND_MAX;
  }  
}


// function definitions
void compute_forces(Body bodies[])
{
  int ii;
  double delta;

  // reset forces
  for (ii = 0; ii < N; ++ii) {
    bodies[ii].Fx = 0.0;
    bodies[ii].Fy = 0.0;
  }

  // add gravitational force
  for (ii = 0; ii < N; ++ii) {
    bodies[ii].Fy += -bodies[ii].mass*G;
  }

  // add force with bottom wall y
  for (ii = 0; ii < N; ++ii) {
    delta = bodies[ii].rad - bodies[ii].Ry;
    if (delta > 0) {
      bodies[ii].Fy += K*delta;
     }
   }

   // add force with bottom wall x left
  for (ii = 0; ii < N; ++ii) {
    delta = bodies[ii].rad - bodies[ii].Rx;
    if (delta > 0) {
      bodies[ii].Fx += K*delta;
     }
   }
   

 
    // add force with bottom wall x rigth
   for (ii = 0; ii < N; ++ii) {
    delta = bodies[ii].rad - (10-bodies[ii].Rx);
     if (delta > 0) {
       bodies[ii].Fx += -K*delta;
      }
    }
}


void start(Body bodies[], double dt)
{
  int ii;
  for (ii = 0; ii < N; ++ii) {
    bodies[ii].arranque(dt);
  }
}

void evolve(Body bodies[], double dt)
{
  int ii;
  for (ii = 0; ii < N; ++ii) {
    bodies[ii].timestep(dt);
  }
}


void init_gnuplot(void)
{
  std::cout << "set size ratio -1" << std::endl;
  std::cout << "set parametric" << std::endl;
  std::cout << "set trange [0:1]" << std::endl;
  std::cout << "set xrange [-1:10]" << std::endl;
  std::cout << "set yrange [-1:11]" << std::endl;
}


void print_to_gnuplot(Body bodies[])
{
  std::cout << "plot "; 
  for (int ii = 0; ii < N; ++ii) {
    std::cout << bodies[ii].Rx << " + " << bodies[ii].rad << "*cos(t*2*pi) , "
	      << bodies[ii].Ry << " + " << bodies[ii].rad << "*sin(t*2*pi) , ";
  }
  std::cout << " 0, 0"; 
  std::cout << std::endl;
}

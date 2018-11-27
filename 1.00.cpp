
  for (int it = 0; it < 2000; ++it) {
    fout << DT*it << "  " << bodies[0].Ry << std::endl; 
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

  // add force with bottom wall
  for (ii = 0; ii < N; ++ii) {
    delta = bodies[ii].rad - bodies[ii].Ry;
    if (delta > 0) {
      bodies[ii].Fy += K*delta;
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
  std::cout << "set xrange [-1:4]" << std::endl;
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

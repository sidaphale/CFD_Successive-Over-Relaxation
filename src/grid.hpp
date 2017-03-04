//Header file to generate a structured grid

#ifndef GRIDHEADERDEF
#define GRIDHEADERDEF

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

class grid
{
private:
  int N; //Number of nodes in both x and y direction
  double xmax, xmin; //Upper and Lower bound in x direction
  double ymax, ymin; //Upper and lower bound in y direction
  double dx, dy;
  std::vector<double> grid_x;
  std::vector<double> grid_y;

public:
  
  grid()
  {
    N = 0.0;
    xmax = 0.0; xmin = 0.0;
    ymax = 0.0; ymin = 0.0;
    dx = 0.0; dy = 0.0;
    grid_x.resize(N, 0.0);
    grid_y.resize(N, 0.0);
    
  }

  grid(int nodes, double x_max, double x_min, double y_max, double y_min) : N(nodes), xmax(x_max), xmin(x_min), ymax(y_max), ymin(y_min)
  {
    grid_x.resize(N, 0.0);
    grid_y.resize(N, 0.0);

    dx = (xmax - xmin) / (N - 1);
    dy = (ymax - ymin) / (N - 1);

    grid_x[0] = xmin;
    grid_y[0] = ymin;

    for (int i = 1; i < N; i++)
      {
	grid_x[i] = grid_x[i-1] + dx;
	grid_y[i] = grid_y[i-1] + dy;
      }
  }

  std::vector<double> GetXGrid()
  {
    return grid_x;
  }

  std::vector<double> GetYGrid()
  {
    return grid_y;
  }

  double GridSpacingX()
  {
    return (grid_x[1] - grid_x[0]);
  }

  double GridSpacingY()
  {
    return (grid_y[1] - grid_y[0]);
  }
  
};

#endif
    

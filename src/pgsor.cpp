// pgsor.cpp : Defines the entry point for the console application.
//

#define _USE_MATH_DEFINES
#include <cmath>
//#include "stdafx.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include "grid.hpp"

using namespace std;


int main(int argc, char* argv[])
{
  std::vector<double> N{ 40, 50 };
	for (int L = 0; L < N.size(); L++)
	{
	  //double p = 0.0, q = 0.0;
		std::vector<double> x;
		std::vector<double> y;
		grid A(N[L], 2*M_PI, 0, 2*M_PI, 0);

		x = A.GetXGrid();
		y = A.GetYGrid();

	        
		/*p = double(2 * M_PI) / double(N[L] - 1);
		q = double(2) / double(N[L] - 1);

		for (int i = 1; i < x.size(); i++)
		{
			x[i] = x[i - 1] + p;
			y[i] = y[i - 1] + q;
			}*/
		/*for (int i = 0; i < y.size(); i++)
		{
			std::cout << y[i] << "\n";
			}*/

		//Calculting the Grid Spacing
	        double del_x, del_y;
		del_x = A.GridSpacingX();
		del_y = A.GridSpacingY();

		//Calculating the value of Beta
		std::vector<double> B(L, 0.0);
		B.push_back((double(pow(del_x, 2))) / ((double(1) - double(0.5*0.5))*double(pow(del_y, 2))));
		//std::cout << "Value of B = " << B[L] << "\n";

		//Calculating the Analytical Solution
		std::vector<std::vector<double>> phi(x.size(), std::vector<double>(y.size(), 0.0));
		for (int i = 0; i < x.size(); i++)
		{
			for (int j = 0; j < y.size(); j++)
			{
				phi[i][j] = (double(-1 * 0.1*cos(x[i])) / double(sqrt(1 - 0.5)))*exp(-y[j] * sqrt(1 - 0.5));
			}
		}

		//Calculating the Aold values
		std::vector<std::vector<double>> Aold(x.size(), std::vector<double>(y.size(), 0.0));
		std::vector<std::vector<double>> Anew(x.size(), std::vector<double>(y.size(), 0.0));
		std::vector<double> d(L, 0.0);
		d.push_back(1 / (2 * (1 + B[L])));
		double w = 1.68; //Over relaxation factor
		std::vector<std::vector<double>> Asor(x.size(), std::vector<double>(y.size(), 0.0));
		std::vector<std::vector<double>> Aerr(x.size(), std::vector<double>(y.size(), 0.0));
		double error = 0.0;
		std::vector<std::vector<double>> u(x.size(), std::vector<double>(y.size(), 0.0));
		std::vector<std::vector<double>> v(x.size(), std::vector<double>(y.size(), 0.0));

		for (int k = 0; k < 100000; k++)
		{
			//1st point of the grid at x=0 and y=0
			for (int i = 0; i < 1; i++)
			{
				for (int j = 0; j < 1; j++)
				{
					Anew[i][j] = d[L] * (2 * Aold[i][j + 1] + 2 * B[L] * Aold[i + 1][j] - 2 * B[L] * del_y * 0.1*cos(x[0]));
				}
			}

			//Boundary points from Pt. 2 to N-1
			for (int i = 0; i < 1; i++)
			{
				for (int j = 1; j < y.size() - 1; j++)
				{
					Anew[i][j] = d[L] * (Aold[i][j + 1] + Anew[i][j - 1] + 2 * B[L] * Aold[i + 1][j] - 2 * B[L] * del_y * 0.1*cos(x[j]));
				}
			}

			//Point (0,N)
			for (int i = 0; i < 1; i++)
			{
				for (int j = y.size() - 1; j < y.size(); j++)
				{
					Anew[i][j] = d[L] * (2 * Anew[i][j - 1] + 2 * B[L] * Aold[i + 1][j] - 2 * B[L] * del_y * 0.1*cos(x[j]));
				}
			}

			//x=0 boundary to y = 2 to N-1
			for (int i = 1; i < x.size() - 1; i++)
			{
				for (int j = 0; j < 1; j++)
				{
					Anew[i][j] = d[L] * (2 * Aold[i][j + 1] + B[L] * Aold[i + 1][j] + B[L] * Anew[i - 1][j]);
				}
			}

			//x=0 boundary to y = N
			for (int i = x.size() - 1; i < x.size(); i++)
			{
				for (int j = 0; j < 1; j++)
				{
					Anew[i][j] = d[L] * (2 * Aold[i][j + 1] + 2 * B[L] * Anew[i - 1][j]);
				}
			}

			//For all the inner grid points
			for (int i = 1; i < x.size() - 1; i++)
			{
				for (int j = 1; j < y.size() - 1; j++)
				{
					Anew[i][j] = d[L] * (Aold[i][j + 1] + Anew[i][j - 1] + B[L] * Aold[i + 1][j] + B[L] * Anew[i - 1][j]);
				}
			}

			//x=N boundary, y=2 to N-1 points
			for (int i = 1; i < x.size() - 1; i++)
			{
				for (int j = y.size() - 1; j < y.size(); j++)
				{
					Anew[i][j] = d[L] * (2 * Anew[i][j - 1] + B[L] * Aold[i + 1][j] + B[L] * Anew[i - 1][j]);
				}
			}

			//y=N boundary, x = 2 to N-1 grid points
			for (int i = x.size() - 1; i < x.size(); i++)
			{
				for (int j = 1; j < y.size() - 1; j++)
				{
					Anew[i][j] = d[L] * (Aold[i][j + 1] + Anew[i][j - 1] + 2 * B[L] * Anew[i - 1][j]);
				}
			}

			//y=N, x=N point
			Anew[x.size() - 1][y.size() - 1] = d[L] * (2 * Anew[x.size() - 2][y.size() - 2] + 2 * B[L] * Anew[x.size() - 1][y.size() - 1]);

			//Introducing the Over Relaxation term

			for (int i = 0; i < x.size(); i++)
			{
				for (int j = 0; j < y.size(); j++)
				{
					Asor[i][j] = (1 - w)*Aold[i][j] + w*Anew[i][j];
				}
			}

			//Calculating the error
			for (int i = 0; i < x.size(); i++)
			{
				for (int j = 0; j < y.size(); j++)
				{
					Aerr[i][j] = std::abs(Asor[i][j] - Aold[i][j]);
				}
			}

			//Checking the convergence criteria
			double sum1 = 0.0, sum2 = 0.0;
			for (int i = 0; i < x.size(); i++)
			{
				for (int j = 0; j < y.size(); j++)
				{
					sum1 += Aerr[i][j];
					sum2 += std::abs(Aold[i][j]);
				}
			}
			
			//Calculating the error 
			error = sum1 / sum2;
			

			if (error < 1.e-8)
			{
				std::cout << "The solution is converged" << "\n";
				break;
			}

			for (int i = 0; i < x.size(); i++)
			{
				for (int j = 0; j < y.size(); j++)
				{
					Aold[i][j] = Asor[i][j];
				}
			}
		}

		//Writing the values of phi and Asor to a dat file

		std::ofstream write_phi("phi.dat");
		write_phi.setf(std::ios::scientific);
		write_phi.precision(16);

		assert(write_phi.is_open());
		for (int i = 0; i < x.size(); i++)
		{
			for (int j = 0; j < y.size(); j++)
			{
				write_phi << phi[i][j] << " ";
			}
			write_phi << "\n";
		}
		write_phi.close();

		std::ofstream write_asor("Asor.dat");
		write_asor.setf(std::ios::scientific);
		write_asor.precision(16);

		assert(write_asor.is_open());
		for (int i = 0; i < x.size(); i++)
		{
			for (int j = 0; j < y.size(); j++)
			{
				write_asor << Asor[i][j] << " ";
			}
			write_asor << "\n";
		}
		write_asor.close();

		
		//Calculating the inner velocities u and v
		for (int i = 1; i < x.size() - 1; i++)
		{
			for (int j = 0; j < y.size(); j++)
			{
				u[i][j] = (Asor[i + 1][j] - Asor[i - 1][j]) / (2 * del_x);
			}
		}
		for (int i = 0; i < x.size(); i++)
		{
			for (int j = 1; j < y.size() - 1; j++)
			{
				v[i][j] = (Asor[i][j + 1] - Asor[i][j - 1]) / (2 * del_y);
			}
		}

		//Calculating the velocities at the boudaries

		//u velocities
		for (int i = 0; i < 1; i++)
		{
			for (int j = 0; j < y.size(); j++)
			{
				u[i][j] = (Asor[i + 1][j] - Asor[i][j]) / (del_x);
			}
		}
		for (int i = x.size() - 1; i < x.size(); i++)
		{
			for (int j = 0; j < y.size(); j++)
			{
				u[i][j] = (Asor[i][j] - Asor[i - 1][j]) / (del_x);
			}
		}

		//v velocities	
		for (int i = 0; i < x.size(); i++)
		{
			for (int j = 0; j < 1; j++)
			{
				v[i][j] = (Asor[i][j + 1] - Asor[i][j]) / (del_y);
			}
		}
		for (int i = 0; i < x.size(); i++)
		{
			for (int j = y.size() - 1; j < y.size(); j++)
			{
				v[i][j] = (Asor[i][j] - Asor[i][j - 1]) / (del_y);
			}
		}

		//Writing the velocities u and v to dat file
		std::ofstream write_u("u.dat");
		write_u.setf(std::ios::scientific);
		write_u.precision(16);

		assert(write_u.is_open());
		for (int i = 0; i < x.size(); i++)
		{
			for (int j = 0; j < y.size(); j++)
			{
				write_u << u[i][j] << " ";
			}
			write_u << "\n";
		}
		write_u.close();

		std::ofstream write_v("v.dat");
		write_v.setf(std::ios::scientific);
		write_v.precision(16);

		assert(write_v.is_open());
		for (int i = 0; i < x.size(); i++)
		{
			for (int j = 0; j < y.size(); j++)
			{
				write_v << v[i][j] << " ";
			}
			write_v << "\n";
		}
		write_v.close();

		std::ofstream write_x("x.dat");
		write_x.setf(std::ios::scientific);
		write_x.precision(16);

		assert(write_x.is_open());
		for (int i = 0; i < x.size(); i++)
		{
			write_x << x[i] << " ";
		}
		write_x.close();

		std::ofstream write_y("y.dat");
		write_y.setf(std::ios::scientific);
		write_y.precision(16);

		assert(write_y.is_open());
		for (int i = 0; i < y.size(); i++)
		{
			write_y << y[i] << " ";
		}
		write_y.close();



	}

	return 0;
}


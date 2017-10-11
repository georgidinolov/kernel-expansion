## This script generates figures which compare the profile of the
## likelihood with respect to resolution in the ``small'' direction
rm(list=ls());
library("mvtnorm");
source("../2-d-solution.R");
source("../../classical-solution/2-d-solution.R");

PLOT.SOLUTION = TRUE;
dx = 1.0/400
dy = 1.0/400
K.prime = 12;

sigma_x=0.001320 * 1.0; ## ``small'' direction
sigma_y=0.032554 * 1.0; ## ``big'' direction

ax=-0.018268;
x_T=-0.017759;
bx=0.007022; 
## ##
ay=-0.014171;
y_T=-0.013896;
by=0.008332;

Lx = bx-ax
Ly = by-ay

problem.parameters.generate.data = NULL;
problem.parameters.generate.data$t <- 2.0 * ( (sigma_y/Ly)^2 );
problem.parameters.generate.data$sigma.2.x <- ( (sigma_x/Lx)^2 ) / ( (sigma_y/Ly)^2 ); 
problem.parameters.generate.data$sigma.2.y <- ( (sigma_y/Ly)^2 ) / ( (sigma_y/Ly)^2 );
problem.parameters.generate.data$rho <- 0.6;
problem.parameters.generate.data$x.ic <- 0;
problem.parameters.generate.data$y.ic <- 0;
dt <- problem.parameters.generate.data$t/1000;

datum = NULL;
datum$ax = ax/Lx;
datum$x.fc = x_T/Lx;
datum$bx = bx/Lx;

datum$ay = ay/Ly;
datum$y.fc = y_T/Ly;
datum$by = by/Ly;

problem.parameters <- datum;

x <- seq(0,1,by=dx);
y <- seq(0,1,by=dy);

sigma.x=0.40;
sigma.y=0.40;
std.dev.factor=1;

sigma2.x=sigma.x^2;
sigma2.y=sigma.y^2;
l=1;
function.list <-
    basis.functions.normal.kernel(rho=problem.parameters.generate.data$rho,
                                  l=l,
                                  sigma2.x=sigma2.x,
                                  sigma2.y=sigma2.y,
                                  dx,dy,
                                  std.dev.factor=std.dev.factor);
orthonormal.function.list <- function.list;
orthonormal.function.list <- orthonormal.functions(function.list,
                                                   dx,dy,x,y,
                                                   TRUE);
system.mats <- system.matrices(orthonormal.function.list,
                               dx,dy);


    a.indeces = c( 0, -1);
    b.indeces = c(-1,  1);
    c.indeces = c( 0, -1);
    d.indeces = c(-1,  1);
    
    a.power=0;
    b.power=0;
    c.power=0;
d.power=0;

h.ay = exp(-2.7);
    h.by = exp(-2.7);
    h.bx = exp(-2.7);
h.ax = exp(-2.7);

                derivative = 0;
                a.power=0;
                b.power=0;
                c.power=0;
                d.power=0;
                
                for ( i in c(1,2)) {
                    if (i==1) { a.power=1; } else { a.power=0; };
                    
                    for ( j in c(1,2)) {
                        if (j==1) { b.power=1; } else { b.power=0; };
                        
                        for ( k in c(1,2)) {
                            if (k==1) { c.power=1; } else { c.power=0; };
                            
                            for ( l in c(1,2)) {
                                if (l==1) { d.power=1; } else { d.power=0; };
                                
                                problem.parameters.original <- datum
                                
                                problem.parameters.original$ax <-
                                    problem.parameters.original$ax + a.indeces[i]*h.ax;
                                problem.parameters.original$bx <-
                                    problem.parameters.original$bx + b.indeces[j]*h.bx;
                                
                                problem.parameters.original$ay <-
                                    problem.parameters.original$ay + c.indeces[k]*h.ay;
                                problem.parameters.original$by <-
                                    problem.parameters.original$by + d.indeces[l]*h.by;
                                
                                current.sol  <- blackbox(function.list,
                                                         orthonormal.function.list,
                                                         system.mats,
                                                         problem.parameters.original,
                                                         dx,dy,
                                                         FALSE, FALSE);
                                
                                L.x <- (problem.parameters.original$bx-
                                        problem.parameters.original$ax);
                                L.y <- (problem.parameters.original$by-
                                        problem.parameters.original$ay);
                                
                                current.sol = current.sol *
                                    1.0/(L.x * L.y);
                                
                                derivative = derivative +
                                    current.sol * ( (-1)^a.power *
                                                         (-1)^b.power *
                                                              (-1)^c.power *
                                                                   (-1)^d.power );
                            }
                        }
                    }
                }

                derivative <- derivative /
                    (h.ax * 2*h.bx * h.ay * 2*h.by);
                print(c(derivative, h));

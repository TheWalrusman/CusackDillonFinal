#ifndef BOUNDARYFUNCTIONS_H
#define BOUNDARYFUNCTIONS_H

//lower Dirichlet boundary function for x = 0
double xLower (const double x, const double y);

//upper Dirichlet boundary function for x = 1
double xUpper (const double x, const double y);

//lower Dirichlet boundary function for y = 0
double yLower (const double x, const double y);

//upper Dirichlet boundary function for y = 1
double yUpper (const double x, const double y);

#endif
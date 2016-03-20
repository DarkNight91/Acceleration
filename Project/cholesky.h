#ifndef CHOLESKY_H
#define CHOLESKY_H

int cholesky(double **orig, int n, float **aug, int mcol, double **chol, float **cholaug, int ofs);
int cholesky(float **orig, int n, float **aug, int mcol, float **chol, float **cholaug, int ofs);

#endif

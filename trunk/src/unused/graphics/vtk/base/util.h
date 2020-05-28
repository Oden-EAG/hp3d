#ifndef UTIL_H
#define UTIL_H

int hexa(int,int,int,int);
int nrpt_hexa(int);

int tetra(int,int,int,int);
int nrpt_tetra(int);

int prism(int,int,int,int);
int nrpt_prism(int);

int pyramid(int,int,int,int);
int nrpt_pyramid(int);

#define max(a,b) ((b<a)?(a):(b))
#define min(a,b) ((b<a)?(b):(a))
#define abs(a)   ((a<0)?(-(a)):(a))
#define mag(a,b) sqrt((a)*(a)+(b)*(b))

#endif

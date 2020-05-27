
#include "util.h"

int hexa(int n,int i,int j,int k){
  return (k*(n+1)+j)*(n+1)+i;  
}

int nrpt_hexa(int n){
  return (n+1)*(n+1)*(n+1);
}

int tetra(int n,int i,int j,int k){
  return (3*k*(n+1)*(n+2)+k*(k-1)*(k-3*n-5))/6+j*(n+1-k)-j*(j-1)/2+i;
}

int nrpt_tetra(int n){
  return (n+1)*(n+2)*(n+3)/6;
}

int prism(int n,int i,int j,int k){
  return k*(n+1)*(n+2)/2 + (2*n+3-j)*j/2 + i;
}

int nrpt_prism(int n){
  return (n+1)*(n+1)*(n+2)/2;
}

int pyramid(int n,int i,int j,int k){
  return (n+2)*(n+1-k)*k+k*(k+1)*(2*k+1)/6 + (n+1-k)*j + i;
}

int nrpt_pyramid(int n){
  return (n+1)*(n+2)*(2*n+3)/6;
}

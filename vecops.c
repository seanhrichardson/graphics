#include <stdio.h>
#include <FPT.h>

void add_vec(double ret[3], double v[3], double u[3]) {
  int i;
  for (i = 0; i < 3; i++) {
    ret[i] = v[i] + u[i];
  }
}
void scale_vec(double v[3], double u[3], double s) {
//v = s*u
  int i;
  for (i = 0; i < 3; i++) {
    v[i] = u[i]*s;
  }
}
void subtract_vec(double ret[3], double v[3], double u[3]) {
  //ret = v-u
  int i;
  for (i = 0; i < 3; i++) {
    ret[i] = v[i]-u[i];
  }
}

double d_prod(double u[3], double v[3]) {
  int i;
  double sum = 0;
  for (i = 0; i < 3; i++) {
    sum += u[i]*v[i];
  }
  return sum;
}

void projection(double ret[3], double u[3], double v[3]) {
  //projection of u onto v
  double factor = d_prod(u,v)/d_prod(v,v);
  scale_vec(ret, v, factor);
}

void orth_part(double ret[3], double u[3], double v[3]) {
  //part of u that is orthogonal to v
  double proj[3];
  projection(proj, u, v);
  subtract_vec(ret, u, proj);
}

void normalize(double v[3]) {
  int i;
  double mag = sqrt(d_prod(v,v));
  for (i = 0; i < 3; i++) {
    v[i] = v[i]/mag;
  }
}
void x_prod(double cross[3], double u[3], double v[3]) {
  cross[0] = u[1]*v[2]-u[2]*v[1];
  cross[1] = u[2]*v[0]-u[0]*v[2];
  cross[2] = u[0]*v[1]-u[1]*v[0];
}
void copy_vec(double ret[3], double u[3]) {
    ret[0] = u[0]; ret[1] = u[1]; ret[2] = u[2];
}
double mag_vec(double u[3]) {
    return sqrt(d_prod(u,u));
}
void print_vec(double u[3]) {
    printf("[%lf,%lf,%lf]",u[0],u[1],u[2]);
}

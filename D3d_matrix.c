#include <D3d_matrix.h>

int D3d_print_mat (double a[4][4]) {
  int r,c;
  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      printf(" %12.4lf ",a[r][c]);
    }
    printf("\n");
  }
}
int D3d_copy_mat (double a[4][4], double b[4][4]) {
// a = b
  int r,c;
  for (r = 0; r < 4; r++) {
      for (c = 0; c < 4; c++ ) {
           a[r][c] = b[r][c];
      }
  }
  return 1;
}
int D3d_mat_mult (double res[4][4], double a[4][4], double b[4][4]) {
// res = a * b
// this is SAFE, i.e. the user can make a call such as 
// D3d_mat_mult(p,  p,q) or D3d_mat_mult(p,  q,p) or  D3d_mat_mult(p, p,p)
  int r,c,i;
  double sum;
  double a_copy[4][4];
  double b_copy[4][4];
  D3d_copy_mat(a_copy,a);
  D3d_copy_mat(b_copy,b);
  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      sum = 0;
      for (i = 0; i < 4; i++) {
	sum += a_copy[r][i]*b_copy[i][c];
      }
      res[r][c] = sum;
    }
  }
  return 1;
}
int D3d_make_identity (double a[4][4]) {
  // a = I
  int r,c ;
  for (r = 0; r < 4; r++) {
    for (c = 0; c < 4; c++) {
      if (r == c) a[r][c] = 1.0;
      else        a[r][c] = 0.0;
    }
  }
  return 1;
}
int D3d_translate (double a[4][4], double b[4][4], 
                   double dx, double dy, double dz) {
// a = translation*a  
// b = b*translation_inverse
  double t[4][4] ;

  D3d_make_identity(t);

  t[0][3] = dx; t[1][3] = dy; t[2][3] = dz;
  D3d_mat_mult(a,t,a);

  t[0][3] = -dx; t[1][3] = -dy; t[2][3] = -dz;
  D3d_mat_mult(b,b,t);

  return 1 ;
}

int D3d_scale (double a[4][4], double b[4][4],
               double sx, double sy, double sz) {
// a = scale*a  
// b = b*scale_inverse
  double a_scale[4][4];
  D3d_make_identity(a_scale);
  double b_scale[4][4];
  D3d_make_identity(b_scale);
  a_scale[0][0] = sx;   a_scale[1][1] = sy;   a_scale[2][2] = sz;
  b_scale[0][0] = 1/sx; b_scale[1][1] = 1/sy; b_scale[2][2] = 1/sz;
  D3d_mat_mult(a,a_scale,a);
  D3d_mat_mult(b,b_scale,b);
  return 1;
}

int D3d_rotate_x (double a[4][4], double b[4][4], double radians) {
// a = rotate*a  
//  b = b*rotate_inverse
  double rotate[4][4];
  double rotate_inv[4][4];
  D3d_make_identity(rotate);
  double cs = cos(radians);
  double sn = sin(radians);
  rotate[1][1] = cs; rotate[1][2] = -sn;
  rotate[2][1] = sn; rotate[2][2] =  cs;
  D3d_make_identity(rotate_inv);
  rotate_inv[1][1] =  cs; rotate_inv[1][2] = sn;
  rotate_inv[2][1] = -sn; rotate_inv[2][2] = cs;
  D3d_mat_mult(a,rotate,a);
  D3d_mat_mult(b,b,rotate_inv);
  return 1;
  
}
int D3d_rotate_y (double a[4][4], double b[4][4], double radians) {
// a = rotate*a  
// b = b*rotate_inverse  
  double rotate[4][4];
  double rotate_inv[4][4];
  D3d_make_identity(rotate);
  double cs = cos(radians);
  double sn = sin(radians);
  rotate[0][0] =  cs; rotate[0][2] = sn;
  rotate[2][0] = -sn; rotate[2][2] = cs;
  D3d_make_identity(rotate_inv);
  rotate_inv[0][0] = cs; rotate_inv[0][2] = -sn;
  rotate_inv[2][0] = sn; rotate_inv[2][2] =  cs;
  D3d_mat_mult(a,rotate,a);
  D3d_mat_mult(b,b,rotate_inv);
  return 1;
}

int D3d_rotate_z (double a[4][4], double b[4][4], double radians) {
// a = rotate*a  
// b = b*rotate_inverse  
  double rotate[4][4];
  double rotate_inv[4][4];
  D3d_make_identity(rotate);
  double cs = cos(radians);
  double sn = sin(radians);
  rotate[0][0] = cs; rotate[0][1] = -sn;
  rotate[1][0] = sn; rotate[1][1] =  cs;
  D3d_make_identity(rotate_inv);
  rotate_inv[0][0] =  cs; rotate_inv[0][1] = sn;
  rotate_inv[1][0] = -sn; rotate_inv[1][1] = cs;
  D3d_mat_mult(a,rotate,a);
  D3d_mat_mult(b,b,rotate_inv);
  return 1;
}
int D3d_cs_rotate_x (double a[4][4], double b[4][4], double cs, double sn)
// a = rotate*a  
// b = b*rotate_inverse  
{
  double t[4][4] ;

  D3d_make_identity(t) ;

  t[1][1] =   cs ;  t[1][2] = -sn ;
  t[2][1] =   sn ;  t[2][2] =  cs ;
  D3d_mat_mult(a,  t,a) ;

  t[1][1] =   cs ;  t[1][2] =  sn ;
  t[2][1] =  -sn ;  t[2][2] =  cs ;
  D3d_mat_mult(b,  b,t) ;

  return 1 ;
}



int D3d_cs_rotate_y (double a[4][4], double b[4][4], double cs, double sn)
// a = rotate*a  
// b = b*rotate_inverse  
{
  double t[4][4] ;

  D3d_make_identity(t) ;

  t[0][0] =   cs ;  t[0][2] =  sn ;
  t[2][0] =  -sn ;  t[2][2] =  cs ;
  D3d_mat_mult(a,  t,a) ;

  t[0][0] =   cs ;  t[0][2] = -sn ;
  t[2][0] =   sn ;  t[2][2] =  cs ;
  D3d_mat_mult(b,  b,t) ;

  return 1 ;
}



int D3d_cs_rotate_z (double a[4][4], double b[4][4], double cs, double sn)
// a = rotate*a  
// b = b*rotate_inverse  
{
  double t[4][4] ;

  D3d_make_identity(t) ;

  t[0][0] =   cs ;  t[0][1] = -sn ;
  t[1][0] =   sn ;  t[1][1] =  cs ;
  D3d_mat_mult(a,  t,a) ;

  t[0][0] =   cs ;  t[0][1] =  sn ;
  t[1][0] =  -sn ;  t[1][1] =  cs ;
  D3d_mat_mult(b,  b,t) ;

  return 1 ;
}
int D3d_negate_x (double a[4][4], double b[4][4]) {
// negate the x....reflects in the yz-plane
// a = reflect*a 
// b = b*reflect_inverse  
}

int D3d_negate_y (double a[4][4], double b[4][4]) {
// negate the y....reflects in the xz-plane
// a = reflect*a 
// b = b*reflect_inverse  
}

int D3d_negate_z (double a[4][4], double b[4][4]) {
// negate the z....reflects in the xy-plane
// a = reflect*a 
// b = b*reflect_inverse
}
int D3d_mat_mult_points (double *X, double *Y, double *Z,
                         double m[4][4],
                         double *x, double *y, double *z, int numpoints) {
  // |X0 X1 X2 ...|       |x0 x1 x2 ...|
  // |Y0 Y1 Y2 ...| = m * |y0 y1 y2 ...|
  // |Z0 Z1 Z2 ...|       |z0 z1 z2 ...|
  // | 1  1  1 ...|       | 1  1  1 ...|
  
  // SAFE, user may make a call like D3d_mat_mult_points (x,y,z, m, x,y,z, n);
  int i;
  double x_copy[numpoints];
  double y_copy[numpoints];
  double z_copy[numpoints];
  for (i = 0; i < numpoints; i++) {
    x_copy[i] = x[i];
    y_copy[i] = y[i];
    z_copy[i] = z[i];
  }
  for (i = 0; i < numpoints; i++) {
    X[i] = m[0][0]*x_copy[i] + m[0][1]*y_copy[i] + m[0][2]*z_copy[i] + m[0][3]; 
    Y[i] = m[1][0]*x_copy[i] + m[1][1]*y_copy[i] + m[1][2]*z_copy[i] + m[1][3]; 
    Z[i] = m[2][0]*x_copy[i] + m[2][1]*y_copy[i] + m[2][2]*z_copy[i] + m[2][3]; 
  }
  return 1;
}

int D3d_mat_mult_pt (double P[3],
                     double m[4][4],
                     double Q[3]) {
  double q_copy[] = {Q[0],Q[1],Q[2]};
  P[0] = m[0][0]*q_copy[0] + m[0][1]*q_copy[1] + m[0][2]*q_copy[2] + m[0][3]; 
  P[1] = m[1][0]*q_copy[0] + m[1][1]*q_copy[1] + m[1][2]*q_copy[2] + m[1][3]; 
  P[2] = m[2][0]*q_copy[0] + m[2][1]*q_copy[1] + m[2][2]*q_copy[2] + m[2][3];
// multiplies a SINGLE point by a matrix
// | P[0] |       | Q[0] |
// | P[1] | = m * | Q[1] |
// | P[2] |       | Q[2] |
// |  1   |       |  1   |

// SAFE, user may make a call like 
// D3d_mat_mult_pt (q,  m,q) ;

}
int D3d_transpose(double ret[4][4], double m[4][4]) {
    int i, j;
    D3d_make_identity(ret);
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
           ret[i][j] = m[j][i]; 
        }
    } 
}
int D3d_x_product (double res[3], double a[3], double b[3]) {
  res[0] = a[1]*b[2]-a[2]*b[1];
  res[1] = a[2]*b[0]-a[0]*b[2];
  res[2] = a[0]*b[1]-a[1]*b[0];
// res = a x b  , cross product of two vectors
// SAFE: it is ok to make a call such as
// D3d_x_product (a,  a,b) or
// D3d_x_product (b,  a,b) or
// D3d_x_product (a,  a,a)
  
}

int D3d_make_movement_sequence_matrix (
                              double mat[4][4],
                              double inv[4][4],
                              int num_movements,
                              int *movement_type_list,
                              double *parameter_list ) {
  
     // create a matrix (mat) and its inverse (inv)
     // that specify a sequence of movements....

     // movement_type_list[k] is an integer that
     // specifies the type of matrix to be used in the
     // the k-th movement.  the parameter that each
     // matrix needs is supplied in parameter_list[k].

     // return 1 if successful, 0 if error


     // the codes for movement_type_list are :
     // 0 - scale x
     // 1 - scale y
     // 2 - scale z

     // 3 - rotate x
     // 4 - rotate y
     // 5 - rotate z
   
     // 6 - translate x
     // 7 - translate y
     // 8 - translate z

     // 9  - negate x...reflect in the yz plane
     // 10 - negate y...relfect in the xz plane
     // 11 - negate z...reflect in the zy plane

  D3d_make_identity(mat);
  D3d_make_identity(inv);
  
  int i;
  for (i = 0; i < num_movements; i++) {
    
    if (movement_type_list[i] == 0) {
      D3d_scale(mat,inv,parameter_list[i],1,1);
    }
    if (movement_type_list[i] == 1) {
      D3d_scale(mat,inv,1,parameter_list[i],1);
    }
    if (movement_type_list[i] == 2) {
      D3d_scale(mat,inv,1,1,parameter_list[i]);
    } //
    if (movement_type_list[i] == 3) {
      D3d_rotate_x(mat,inv,parameter_list[i]*M_PI/180);
    }
    if (movement_type_list[i] == 4) {
      D3d_rotate_y(mat,inv,parameter_list[i]*M_PI/180);
    }
    if (movement_type_list[i] == 5) {
      D3d_rotate_z(mat,inv,parameter_list[i]*M_PI/180);
    } //
    if (movement_type_list[i] == 6) {
      D3d_translate(mat,inv,parameter_list[i],0,0);
    }
    if (movement_type_list[i] == 7) {
      D3d_translate(mat,inv,0,parameter_list[i],0);
    }
    if (movement_type_list[i] == 8) {
      D3d_translate(mat,inv,0,0,parameter_list[i]);
    } //
    if (movement_type_list[i] == 9) {
      D3d_negate_x(mat,inv);
    }
    if (movement_type_list[i] == 10) {
      D3d_negate_y(mat,inv);
    }
    if (movement_type_list[i] == 11) {
      D3d_negate_z(mat,inv);
    }
  }
  return 1;
}

int D3d_viewAA (double view[4][4], double view_inverse[4][4],
                double eye[3], double coi[3], double up[3],
                double a, double b, double c, double p, double r)
// return 1 if successful, 0 otherwise
{
    double Uplen ;
    double Up[3] ;

    //printf("\nAA\n");
    
    if ((p == 0) || (r == 0)) return 0 ;

    D3d_make_identity(view) ;     
    D3d_make_identity(view_inverse) ;     
    D3d_translate(view,view_inverse, -eye[0], -eye[1], -eye[2]) ;
    D3d_cs_rotate_y (view,view_inverse, c/p, -a/p) ;
    D3d_cs_rotate_x (view,view_inverse, p/r,  b/r) ;

    D3d_mat_mult_pt (Up, view, up) ;
    Uplen = sqrt(Up[0]*Up[0] + Up[1]*Up[1]) ;

    //printf("\nUplen:%lf\n",Uplen);

    if (Uplen == 0) return 0 ;

    D3d_cs_rotate_z (view,view_inverse, Up[1]/Uplen, Up[0]/Uplen) ;

    return 1 ;
}

int D3d_viewBB (double view[4][4], double view_inverse[4][4],
                double eye[3], double coi[3], double up[3],
                double a, double b, double c, double p, double r)
// return 1 if successful, 0 otherwise
{
    double Uplen ;
    double Up[3] ;

    //printf("\nBB:\n");

    p = sqrt(b*b + c*c) ; // alter the incoming p

    if ((p == 0) || (r == 0)) return 0 ;

    D3d_make_identity(view) ;     
    D3d_make_identity(view_inverse) ;     
    D3d_translate(view,view_inverse, -eye[0], -eye[1], -eye[2]) ;
    D3d_cs_rotate_x (view,view_inverse, c/p,  b/p) ;
    D3d_cs_rotate_y (view,view_inverse, p/r, -a/r) ;

    D3d_mat_mult_pt (Up, view, up) ;
    Uplen = sqrt(Up[0]*Up[0] + Up[1]*Up[1]) ;

    if (Uplen == 0) return 0 ;

    D3d_cs_rotate_z (view,view_inverse, Up[1]/Uplen, Up[0]/Uplen) ;

    return 1 ;
}

int D3d_view (double view[4][4], double view_inverse[4][4],
              double eye[3], double coi[3], double up[3])
// Construct the view matrix and its inverse given the location
// of the eye, the center of interest, and an up point.
// return 1 if successful, 0 otherwise.
{
    double a,b,c,p,r ;
    int s ;


    // printf("entering D3d_view\n") ;
    a = coi[0] - eye[0] ;
    b = coi[1] - eye[1] ;
    c = coi[2] - eye[2] ;
    p = sqrt(a*a + c*c) ;
    r = sqrt(a*a + b*b + c*c) ;

    //printf("a:%lf - b:%lf - c:%lf - p:%lf - r:%lf",a,b,c,p,r);

    // printf("a,b,c,p,r = %lf %lf %lf %lf %lf\n\n",a,b,c,p,r) ;

    if (fabs(b) < p) {
         // printf("choose AA\n") ;
         s = D3d_viewAA (view, view_inverse,
                         eye, coi, up,
                         a,b,c,p,r) ;
    } else {
         // printf("choose BB\n") ;
         s = D3d_viewBB (view, view_inverse,
                         eye, coi, up,
                         a,b,c,p,r) ;
    }

    if (s == 0) printf("D3d_view error\n") ;

    // printf("leaving D3d_view\n") ;

    return s ;
}

/*int D3d_view (double view[4][4],  double view_inverse[4][4],
              double eye[3], double coi[3], double up[3]) {
// Construct the view matrix and its inverse given the location
// of the eye, the center of interest, and an up point.
// return 1 if successful, 0 otherwise.
  double a = coi[0]-eye[0];
  double b = coi[1]-eye[1];
  double c = coi[2]-eye[2];

  double p = sqrt(pow(a,2)+pow(c,2));
  double r = sqrt(pow(a,2)+pow(b,2)+pow(c,2));

  if (p == 0 || r == 0) return 0;

  //printf("a:%12.4lf\tb:%12.4lf\tc:%12.4lf\t\n",a,b,c);
  //printf("p:%12.4lf\tr:%12.4lf\n",p,r);

  D3d_make_identity(view);
  D3d_make_identity(view_inverse);
  
  //1. translate so eye is in origin
  D3d_translate(view, view_inverse, -eye[0],-eye[1],-eye[2]);
  //2. rotate into yz plane about y
  D3d_cs_rotate_y(view, view_inverse, c/p,-a/p);
  //3. rotate onto z axis about x
  D3d_cs_rotate_x(view, view_inverse, p/r,b/r);

  //printf("before up:\n");
  D3d_print_mat(view);
  
  //4. rotate up' into yz plane about z
  //apply previous transformations to up:
  D3d_mat_mult_pt(up,view,up);
  double u0 = up[0]; double u1 = up[1];
  double ud = sqrt(pow(u0,2)+pow(u1,2));

  //printf("u0:%12.4lf\tu1:%12.4lf\tud:%12.4lf\t\n",u0,u1,ud);
  
  D3d_cs_rotate_z(view, view_inverse, u0/ud,u1/ud);

  //printf("after up:\n");
  D3d_print_mat(view);

}*/

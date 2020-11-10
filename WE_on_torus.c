#include <stdio.h>
#include <xwd_tools.h>
#include <FPT.h>

#include <D3d_matrix.h>


//-------INSTRUCTIONS---------
//  "acom -I. D3d_matrix.c xwd_tools.c sean_final_code.c"
//  This will expect a folder called "waving_torus"
//  Variable below is the increment of u and v values. Lower "res" to increase resolution. 
double RES = 0.02;
//This will change how many waves are displayed:
double WAVES = 6;
//------------------------------

double HA = M_PI/4;
double zbuf[800][800];
double AMB = 0.2;
double DIFF_MAX = 0.5;
double SPEC_POW = 30;
double LIGHT[3] = {100,100,-50};
//double LIGHT[3] = {0,0,-50};
double EYE[3] = {0,0,0};
double IN_COLOR = 0.7;
//xwd:
int map_id = 1;

double phase[8];

//VECOPS FILE:
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

//REST
void init_zbuf() {
    int i,j;
    for (i = 0; i < 800; i++) {
        for (j = 0; j < 800; j++) {
            zbuf[i][j] = 0;
        }
    }
}
void draw_point(double p[3], double rgb[3]) {
    int x = (400/tan(HA))*(p[0]/p[2])+400;
    int y = (400/tan(HA))*(p[1]/p[2])+400;
    int inside = (x >= 0 && x < 800 && y >= 0 && y < 800 && p[2] > 0);
    if (inside && (zbuf[x][y] == 0 || zbuf[x][y] > p[2])) {
        zbuf[x][y] = p[2];
        set_xwd_map_color(map_id,x,y,rgb[0],rgb[1],rgb[2]);
    }
}
/*void draw_marker(double p[3]) {
    int x = (400/tan(HA))*(p[0]/p[2])+400;
    int y = (400/tan(HA))*(p[1]/p[2])+400;
    G_rgb(1,0,1);
    G_fill_circle(x,y,10);
}*/

//LIGHT MODEL:
double light_intensity(double p[3], double n[3]) {
    normalize(n);
    double SPEC_MAX = (1-AMB-DIFF_MAX);
    double e[3], l[3], r[3];
    //set vecs
    subtract_vec(e, EYE,p); normalize(e);
    subtract_vec(l, LIGHT,p); normalize(l);
    scale_vec(r, n, 2*d_prod(n,l));
    subtract_vec(r,r,l);
    //normal direction check
    double diffuse, specular;
    if (d_prod(e,n) < 0 && d_prod(l,n) < 0)
       scale_vec(n,n,-1);
    //calc
    diffuse = DIFF_MAX*d_prod(n,l);
    specular = SPEC_MAX*pow(d_prod(r,e),SPEC_POW);
    if (diffuse < 0) {diffuse = 0;}
    if (specular < 0) {specular = 0;}
    if (d_prod(e,n)*d_prod(l,n) < 0) { 
        specular = 0; diffuse = 0;
    }
    //printf("%lf+%lf+%lf\n",AMB,diffuse,specular);
    return AMB + diffuse + specular;
}
void color(double ret[3], double RGB[3], double inten) {
    double white[3] = {1,1,1};
    double black[3] = {0,0,0};
    double hold[3];
    double scale;
    if (inten > IN_COLOR) {
        subtract_vec(hold, white,RGB);
        scale = (inten-IN_COLOR)/(1-IN_COLOR);
    }
    else {
        subtract_vec(hold, black, RGB);
        scale = (IN_COLOR-inten)/(IN_COLOR);
    }
    scale_vec(ret,hold,scale);
    add_vec(ret, ret, RGB);
}
void color1(double RGB[3], double t) {
    RGB[0] = 1-t;
    RGB[1] = t;
    RGB[2] = 0;
}
void old_sphere_par(double p[3], double u, double v) {
    // 0 < u < 2pi; -pi/2 < v < pi/2;
    p[0] = cos(u)*cos(v);
    p[1] = sin(v);
    p[2] = sin(u)*cos(v);
}
void sphere_par(double p[3], double r, double u, double v) {
    // 0 < u < 2pi; 0 < v < pi; 0 < r < inf;
    p[0] = r*cos(u)*sin(v);
    p[1] = r*sin(u)*sin(v);
    p[2] = r*cos(v);
}
void torus_par(double p[3], double r, double u, double v) {
    double R = 2;
    p[0] = -(r*cos(u)+R)*sin(v);
    p[1] =  (r*cos(u)+R)*cos(v);
    p[2] =  (r*sin(u));
}
double get_r(double t, double u, double v) {
    return (cos(t)+1.5);
}
double r1(double t, double u, double v) {
    return sin(t)*(0.5*sqrt(3/2)*cos(u)+1);
}
double r2(double t, double u, double v) {
   return sin(t)*(0.5*sqrt(5/M_PI)*(3*pow((cos(u)),2)-1)+1); 
}
double full_path(double eye[3], double coi[3], double up[3], double T) {
    double t, t0, t1, t2, t3, t4, t5, t6, t7, ti, dt;

    //PHASE 0:
    ti = 0; t0 = ti + 2*M_PI;
    if (T >= ti && T < t0) {
        t = T-ti;
        dt = t0-ti;
        eye[0] = 0; eye[1] = 0; eye[2] = 5;
        coi[0] = coi[0] + cos(M_PI/2*t/dt); coi[1] = 0; coi[2] = eye[2] - sin(M_PI/2*t/dt); 
        up[0] = eye[0]+0; up[1] = eye[1]+1; up[2] = eye[2]+0;
    }   
    //PHASE 1:
    t1 = t0 + 6*M_PI;
    if (T >= t0 && T < t1) {
        t = T-t0;
        dt = t1-t0;
        eye[0] = 5*sin(M_PI*t/dt); eye[1] = 2*sin(M_PI*t/dt); eye[2] = 5*cos(M_PI*t/dt);
        coi[0] = 0; coi[1] = 0; coi[2] = 0;
        up[0] = eye[0]+0; up[1] = eye[1]+1; up[2] = eye[2]+0;
    }

    //PHASE 2:
    t2 = t1 + 2*M_PI;
    if (T >= t1 && T < t2) {
        t = T-t1;
        dt = t2-t1;
        eye[0] = 0; eye[1] = 0; eye[2] = -5+2.5*t/dt;
        coi[0] = 0; coi[1] = t/dt; coi[2] = 0;
        up[0] = eye[0]+0; up[1] = eye[1]+1; up[2] = eye[2]+0;
    }

    //PHASE 3:
    t3 = t2 + 4*M_PI;
    if (T >= t2 && T < t3) {
        t = T-t2;
        dt = t3-t2;
        eye[0] = 0; eye[1] = 0; eye[2] = -2.5+5*t/dt;
        coi[0] = sin(M_PI*t/dt); coi[1] = cos(M_PI*t/dt); coi[2] = 0; 
        up[0] = eye[0]+0; up[1] = eye[1]+1; up[2] = eye[2]+0;
    }

    //PHASE 4:
    t4 = t3 + 2*M_PI;
    if (T >= t3 && T < t4) {
        t = T-t3;
        dt = t2-t1;
        eye[0] = 0; eye[1] = 0; eye[2] = 2.5+2.5*t/dt;
        coi[0] = 0; coi[1] = -1+t/dt; coi[2] = 0;
        up[0] = eye[0]+0; up[1] = eye[1]+1; up[2] = eye[2]+0;
    }

    //PHASE 5:
    t5 = t4 + 6*M_PI;
    if (T >= t4 && T < t5) {
        t = T-t4;
        dt = t5-t4;
        eye[0] = -5*sin(M_PI*t/dt); eye[1] = -2*sin(M_PI*t/dt); eye[2] = 5*cos(M_PI*t/dt);
        coi[0] = 0; coi[1] = 0; coi[2] = 0;
        up[0] = eye[0]+0; up[1] = eye[1]+1; up[2] = eye[2]+0;
    }

    //PHASE 6:
    t6 = t5 + 2*M_PI;
    if (T >= t5 && T < t6) {
        t = T-t5;
        dt = t6-t5;
        eye[0] = 0; eye[1] = 0; eye[2] = -5+2.5*t/dt;
        coi[0] = 0; coi[1] = 0; coi[2] = eye[2] + 1;
         up[0] = eye[0]+0; up[1] = eye[1]+1; up[2] = eye[2]+0;
    }   

    //PHASE 7:
    t7 = t6 + 4*M_PI;
    if (T >= t6 && T < t7) {
        t = T-t6;
        dt = t7-t6;
        eye[0] = 0; eye[1] = 0; eye[2] = -2.5+5*t/dt;
        coi[0] = 0; coi[1] = 0; coi[2] = 20;
        up[0] = eye[0] + sin(2*M_PI*t/dt);
        up[1] = eye[1] + cos(2*M_PI*t/dt);
        up[2] = eye[2];
    } 
    
    //print_vec(eye); printf("\t"); print_vec(coi); printf("\t"); print_vec(up);
    //if (T >= t7) exit(0);
}
int main() {
    phase[0] = 0;
    phase[1] = phase[0] + 2*M_PI;
    phase[2] = phase[1] + 6*M_PI;
    phase[3] = phase[2] + 2*M_PI;
    phase[4] = phase[3] + 4*M_PI;
    phase[5] = phase[4] + 2*M_PI;
    phase[6] = phase[5] + 6*M_PI;
    phase[7] = phase[6] + 2*M_PI;
    phase[8] = phase[7] + 4*M_PI;

    map_id = create_new_xwd_map(800,800);
    init_zbuf();
    clear_xwd_map(map_id,0,0,0); 
    //general:
    double p[3]; 
    double RGB[3];
    double rgb[3];
    double m[4][4];
    double minv[4][4];
    D3d_make_identity(m);
    D3d_make_identity(minv);

    //ulo, uhi, ustep; v... ; t...
    double ulo = 0; double uhi = 2*M_PI; double ustep = 0.02;
    double vlo = 0; double vhi = 2*M_PI; double vstep = 0.02;
    double tlo = phase[0]; double thi = phase[8]; double tstep = 0.1;
    //int k = 6; int l = 6;
    //r given by u(t,u,v2;

    //for normal calculation:
    double du = 0.001; double dv = 0.001;
    double a[3], b[3], n[3];

    //set matrix:
    D3d_translate(m,minv,0,0,0);
    D3d_scale(m,minv,1,1,1);
    
    //set color:
    RGB[0] = 1; RGB[1] = 0; RGB[2] = 0;

    //view stuff:
    double eye[3], up[3], coi[3];
    double view[4][4], vinv[4][4];
    char fname[100]; int i;

    double u,v,t;
    for (t = tlo; t < thi; t += tstep) { //step through time
        i = t/tstep+0.1;

        if (t >= phase[0] && t < phase[1]) {
            ustep = 0.0015; vstep = 0.0015;
        }
        //circling
        if (t >= phase[1] && t < phase[2]) {
            ustep = 0.001; vstep = 0.001;
        }
        //going in
        if (t >= phase[2] && t < phase[3]) {
            ustep = 0.001; vstep = 0.001;
        }
        //looking around
        if (t >= phase[3] && t < phase[4]) {
            ustep = 0.0006; vstep = 0.0006;
        }
        //going out
        if (t >= phase[4] && t < phase[5]) {
            ustep = 0.001; vstep = 0.001;
        }
        //circle again
        if (t >= phase[5] && t < phase[6]) {
            ustep = 0.001; vstep = 0.001;
        }
        //in:
        if (t >= phase[6] && t < phase[7]) {
            ustep = 0.001; vstep = 0.001;
        }
        //spiral:
        if (t >= phase[7] && t < phase[8]) {
            ustep = 0.0006; vstep = 0.0006;
        }
        ustep = RES; vstep = RES;

        printf("frame %d\n",i);

        clear_xwd_map(map_id,0,0,0); 
        init_zbuf();

        full_path(eye,coi,up,t);

        D3d_make_identity(view);
        D3d_make_identity(vinv);
        D3d_view(view,vinv,eye,coi,up);
       
        D3d_mat_mult(m, view, m);
        D3d_mat_mult_pt(LIGHT,view,LIGHT);

        //test:
        D3d_mat_mult_pt(coi,view,coi);
        //draw_marker(coi);


        for (u = ulo; u <= uhi; u += ustep) { //step through theta
            for (v = vlo; v <= vhi; v += vstep) { //step through phi
                //set p = {x,y,z} point
                //double r_pure = 0.5*cos(0.6*t)*cos(3*u)*cos(3*v)+0.5*cos(t)*cos(5*u)*cos(5*v);
                //double r_pure = cos(t)*cos(3*u)*cos(3*v);
                double r_pure = cos(t)*cos(WAVES*u)*cos(WAVES*v);
                double r = 0.5*r_pure+1;
                torus_par(p,r,u,v);
                D3d_mat_mult_pt(p,m,p);
                //normal:
                torus_par(a,r,u+du,v); D3d_mat_mult_pt(a,m,a);
                torus_par(b,r,u,v+dv); D3d_mat_mult_pt(b,m,b);
                subtract_vec(a,a,p);
                subtract_vec(b,b,p);
                x_prod(n,a,b);
                normalize(n);
                //light and color:
                double intensity = light_intensity(p,n);
                color1(RGB,(r_pure+1)/2);
                color(rgb,RGB,intensity); 
                draw_point(p,rgb);

            }
            //G_wait_key();
        }
        
        sprintf(fname,"waving_torus/pic%04d.xwd",i);

        xwd_map_to_named_xwd_file(map_id,fname);

        D3d_mat_mult(m,vinv,m);
        D3d_mat_mult_pt(LIGHT,vinv,LIGHT);
        //G_wait_key();  
    }
}

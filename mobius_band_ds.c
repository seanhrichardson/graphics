#include <stdio.h>
#include <FPT.h>
#include <D3d_matrix.h>
#include <vecops.h>

double HA = M_PI/4;
double zbuf[800][800];
double AMB = 0.2;
double DIFF_MAX = 0.5;
double SPEC_POW = 30;
double LIGHT[3] = {-100,200,-100};
//double LIGHT[3] = {0,0,-50};
double EYE[3] = {0,0,0};
double IN_COLOR = 0.7;

void init_zbuf() {
    int i,j;
    for (i = 0; i < 800; i++) {
        for (j = 0; j < 800; j++) {
            zbuf[i][j] = 0;
        }
    }
}
void draw_point(double p[3]) {
    int x = (400/tan(HA))*(p[0]/p[2])+400;
    int y = (400/tan(HA))*(p[1]/p[2])+400;
    int inside = (x >= 0 && x < 800 && y >= 0 && y < 800);
    if (inside && (zbuf[x][y] == 0 || zbuf[x][y] > p[2])) {
        zbuf[x][y] = p[2];
        G_point(x,y);
    }
}
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

void green_blue(double rgb[3], double u, double v) {
    u = u*2*M_PI;
    v = v*2*M_PI;
    double m[4][4], minv[4][4];
    double a = 1;
    double b = 2;
    rgb[0] = (a*cos(u)+b)*cos(v);
    rgb[1] = (a*cos(u)+b)*sin(v);
    rgb[2] = a*sin(u);
    D3d_make_identity(m);
    D3d_make_identity(minv);
    D3d_scale(m,minv,1.0/6,1.0/6,0.5);
    D3d_scale(m,minv,.5,.5,.5);
    //D3d_rotate_x(m,minv,M_PI/4);
    //D3d_rotate_y(m,minv,-M_PI/4);
    D3d_translate(m,minv,0,.5,.5);
    D3d_mat_mult_pt(rgb,m,rgb);  
}
void fall(double rgb[3], double u, double v) {
    u = u*2*M_PI;
    v = v*2*M_PI;
    double m[4][4], minv[4][4];
    double a = 1;
    double b = 2;
    rgb[0] = (a*cos(u)+b)*cos(v);
    rgb[1] = (a*cos(u)+b)*sin(v);
    rgb[2] = a*sin(u);
    D3d_make_identity(m);
    D3d_make_identity(minv);
    D3d_scale(m,minv,1.0/6,1.0/6,0.5);
    D3d_scale(m,minv,.5,.5,.5);
    //D3d_rotate_x(m,minv,M_PI/4);
    //D3d_rotate_y(m,minv,-M_PI/4);
    D3d_translate(m,minv,.5,.5,0);
    D3d_mat_mult_pt(rgb,m,rgb);  
}

void fall_mobius(double rgb[3], double u, double v) {
    double m[4][4], minv[4][4];
    rgb[0] = (1 + v/2*cos(u/2))*cos(u);
    rgb[1] = v/2*sin(u/2);
    rgb[2] = (1 + v/2*cos(u/2))*sin(u);
    D3d_make_identity(m);
    D3d_make_identity(minv);
    D3d_scale(m,minv,0.2,0.2,0.2);
    D3d_scale(m,minv,.5,.5,.5);
    //D3d_rotate_x(m,minv,M_PI/4);
    //D3d_rotate_y(m,minv,-M_PI/4);
    D3d_translate(m,minv,.6,.4,0);
    D3d_mat_mult_pt(rgb,m,rgb);  
}

void mobius_green_blue(double rgb[3], double u, double v) {
    double m[4][4], minv[4][4];
    rgb[0] = (1 + v/2*cos(u/2))*cos(u);
    rgb[1] = v/2*sin(u/2);
    rgb[2] = (1 + v/2*cos(u/2))*sin(u);
    D3d_make_identity(m);
    D3d_make_identity(minv);
    D3d_scale(m,minv,1.0/6,1.0/6,0.5);
    D3d_scale(m,minv,.5,.5,.5);
    //D3d_rotate_x(m,minv,M_PI/4);
    //D3d_rotate_y(m,minv,-M_PI/4);
    D3d_translate(m,minv,0,.6,.4);
    D3d_mat_mult_pt(rgb,m,rgb);  
}

void mobius_base(double rgb[3], double u, double v) {
    double h = 0.05;
    if (v <= h && v >= -h) {
        rgb[0] = 0;
        rgb[1] = 0;
        rgb[2] = 1;
    }
}
void mobius_fiber(double rgb[3], double u, double v, double uP) {
    double h = 0.02;
    if (u-uP <= h && u-uP >= -h) {
        rgb[0] = 1;
        rgb[1] = 0;
        rgb[2] = 0;
    }
}

void mobius_point(double rgb[3], double u, double v, double uP) {
    double h = 0.07;
    if (pow((u-uP),2) + pow(v,2) < pow(h,2))
    if (u-uP <= h && u-uP >= -h && v <= h && v >= -h) {
        rgb[0] = 0;
        rgb[1] = 1;
        rgb[2] = 0;
    }
}

void local_base(double rgb[3], double u, double v, double uP) {
    double h = 0.05;
    double d = 0.3;
    if (v <= h && v >= -h && u < uP + d && u > uP - d) {
        rgb[0] = 1;
        rgb[1] = 0;
        rgb[2] = 1;
    }
}

void mobius_area(double rgb[3], double u, double v, double uP) {
    double d = 0.3;
    if (u < uP + d && u > uP - d) {
        rgb[0] += (1-rgb[0])*0.6;
    }
    mobius_fiber(rgb, u, v, uP+d);
    mobius_fiber(rgb, u, v, uP-d);
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
void sphere_par(double p[3], double u, double v) {
    // 0 < u < 2pi; -pi/2 < v < pi/2;
    p[0] = cos(u)*cos(v);
    p[1] = sin(v);
    p[2] = sin(u)*cos(v);
}
void cylinder_par(double p[3], double u, double v) {
    // 0 < u < 2pi
    // -1 < v < 1
    double r = 1;
    double h = 1; //half height
    p[0] = r*cos(u);
    p[1] = h*v;
    p[2] = r*sin(u);
}
void mobius_par(double p[3], double u, double v) {
    p[0] = (1 + v/2*cos(u/2))*cos(u);
    p[1] = v/2*sin(u/2);
    p[2] = (1 + v/2*cos(u/2))*sin(u);
}

void mobius_par2(double p[3], double u, double v) {
    p[0] = (1 + v/2*cos(u/2))*cos(u+M_PI);
    p[1] = v/2*sin(u/2);
    p[2] = (1 + v/2*cos(u/2))*sin(u+M_PI);
}

void torus_par(double p[3], double u, double v) {
    double R = 2;
    double r = 1;
    p[0] = -(r*cos(u)+R)*sin(v);
    p[1] =  (r*cos(u)+R)*cos(v);
    p[2] =  (r*sin(u));
}

int main() {
    G_init_graphics(800,800);
    init_zbuf();
    G_rgb(0,0,0);
    G_rgb(1,1,1);
    G_clear();
    G_rgb(1,1,1);
    
    //general:
    double p[3]; 
    double RGB[3];
    double rgb[3];
    double m[4][4];
    double minv[4][4];
    D3d_make_identity(m);
    D3d_make_identity(minv);

    //ulo, uhi, ustep; v... ; t...
    double ulo = 0; double uhi = 2*M_PI;   double ustep = 0.002;
    double vlo = -1; double vhi = 1; double vstep = 0.002;
    //r given by u(t,u,v);

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
    coi[0] = 0; coi[1] = 0; coi[2] = 0;
    //double r = 5;
    //eye[0] = 5*sin(t); eye[1] = 0; eye[2] = 5*cos(t);
    eye[0] = 0; eye[1] = 1.5; eye[2] = -2;
    up[0] = eye[0]+0; up[1] = eye[1]+0; up[2] = eye[2]+0;

    D3d_make_identity(view);
    D3d_make_identity(vinv);
    D3d_view(view,vinv,eye,coi,up);
   
    D3d_mat_mult(m, view, m);
    D3d_mat_mult_pt(LIGHT,view,LIGHT);

    for (u = ulo; u <= uhi; u += ustep) { //step through theta
        for (v = vlo; v <= vhi; v += vstep) { //step through phi
            //set p = {x,y,z} point
            mobius_par(p,u,v);
            D3d_mat_mult_pt(p,m,p);
            //normal:
            //test
            mobius_par(a,u+du,v); D3d_mat_mult_pt(a,m,a);
            mobius_par(b,u,v+dv); D3d_mat_mult_pt(b,m,b);
            subtract_vec(a,a,p);
            subtract_vec(b,b,p);
            x_prod(n,a,b);
            normalize(n);
            //light and color:
            double uP = (0.7)*(2*M_PI);
            fall_mobius(RGB,u,v);
            mobius_base(RGB,u,v);
            //local_base(RGB,u,v,uP);
            //mobius_area(RGB,u,v,uP);
            //mobius_point(RGB,u,v,uP);
            //color_re(RGB,u,v);
            double intensity = light_intensity(p,n);
            color(rgb,RGB,intensity); 
            G_rgb(rgb[0],rgb[1],rgb[2]);
            draw_point(p);
        }
    }

    sprintf(fname,"images/mobius_point.xwd");
    G_save_image_to_file(fname);
    
    D3d_mat_mult(m,vinv,m);
    D3d_mat_mult_pt(LIGHT,vinv,LIGHT);

    G_wait_key();  

    //ulo, uhi, ustep; v... ; t...
    ulo = 0; uhi = 2*M_PI; ustep = 0.002;
    vlo = -1; vhi = 1; vstep = 0.002;
    //r given by u(t,u,v);

    //for normal calculation:
    du = 0.01; dv = 0.01;

    //set matrix:
    D3d_translate(m,minv,0,0,0);
    D3d_scale(m,minv,1,1,1);
    
    //set color:
    RGB[0] = 0; RGB[1] = 1; RGB[2] = 0;

    //view stuff:
    coi[0] = 0; coi[1] = 0; coi[2] = 0;
    //double r = 5;
    //eye[0] = 5*sin(t); eye[1] = 0; eye[2] = 5*cos(t);
    eye[0] = 0; eye[1] = 1.5; eye[2] = -2;
    up[0] = eye[0]+0; up[1] = eye[1]+0; up[2] = eye[2]+0;

    D3d_make_identity(view);
    D3d_make_identity(vinv);
    D3d_view(view,vinv,eye,coi,up);
   
    D3d_mat_mult(m, view, m);
    D3d_mat_mult_pt(LIGHT,view,LIGHT);

    for (u = ulo; u <= uhi; u += ustep) { //step through theta
        for (v = vlo; v <= vhi; v += vstep) { //step through phi
            //set p = {x,y,z} point
            mobius_par2(p,u,v);
            D3d_mat_mult_pt(p,m,p);
            //normal:
            mobius_par2(a,u+du,v); 
            D3d_mat_mult_pt(a,m,a);
            mobius_par2(b,u,v+dv); 
            D3d_mat_mult_pt(b,m,b);
            subtract_vec(a,a,p);
            subtract_vec(b,b,p);
            x_prod(n,a,b);
            normalize(n);
            //light and color:
            mobius_green_blue(RGB,u,v);
            mobius_base(RGB,u,v);
            //color_re(RGB,u,v);
            double intensity = light_intensity(p,n);
            color(rgb,RGB,intensity); 
            G_rgb(rgb[0],rgb[1],rgb[2]);
            draw_point(p);
        }
    }

    sprintf(fname,"images/mobius_band_ds.xwd");
    G_save_image_to_file(fname);
    
    D3d_mat_mult(m,vinv,m);
    D3d_mat_mult_pt(LIGHT,vinv,LIGHT);

    G_wait_key();  
}

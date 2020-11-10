#include <stdio.h>
#include <FPT.h>

#include <D3d_matrix.h>

//feed u,v vals s.t 0<u<1 and 0<v<1
/*
void cool_color1(double RGB[3], double u, double v, struct Thing3D D) {
    double udiff = D.ur[1]-D.ur[0]; double vdiff = D.vr[1]-D.vr[0]; 
    RGB[0] = (u-D.ur[0])/(udiff); RGB[1] = (v-D.vr[0])/vdiff; RGB[2] = 1-(u-D.ur[0])*(v-D.vr[0])/(udiff*vdiff);
}
void checker(double RGB[3], double u, double v, struct Thing3D D) {
    double udiff = D.ur[1]-D.ur[0]; double vdiff = D.vr[1]-D.vr[0]; 
    int freq = 20;
    int a = (u-D.ur[0])/udiff*freq;
    int b = (v-D.vr[0])/vdiff*freq;
    if ((a+b)%2) {copy_vec(RGB, D.color1);} 
    else {copy_vec(RGB, D.color2);}
}
void plain(double RGB[3], double u, double v, struct Thing3D D) {
    copy_vec(RGB, D.color1);
}
void v_stripes(double RGB[3], double u, double v, struct Thing3D D) {
    double udiff = D.ur[1]-D.ur[0];
    int freq = 20;
    int a = (u-D.ur[0])/udiff*freq;
    if (a%2) {copy_vec(RGB, D.color1);} 
    else {copy_vec(RGB, D.color2);}
}
void u_stripes(double RGB[3], double u, double v, struct Thing3D D) {
    double vdiff = D.vr[1]-D.vr[0];
    int freq = 20;
    int b = (v-D.vr[0])/vdiff*freq;
    if (b%2) {copy_vec(RGB, D.color1);} 
    else {copy_vec(RGB, D.color2);}
}*/

void valentines(double rgb[3], double u, double v) {
    u = u*2*M_PI;
    v = u*2*M_PI;
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
    D3d_translate(m,minv,.5,0,.5);
    D3d_mat_mult_pt(rgb,m,rgb);  
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

void all(double rgb[3], double u, double v) {
    u = u*2*M_PI;
    v = u*2*M_PI;
    double m[4][4], minv[4][4];
    double a = 1;
    double b = 2;
    rgb[0] = (a*cos(u)+b)*cos(v);
    rgb[1] = (a*cos(u)+b)*sin(v);
    rgb[2] = a*sin(u);
    D3d_make_identity(m);
    D3d_make_identity(minv);
    D3d_scale(m,minv,1.0/6,1.0/6,0.5);
    D3d_translate(m,minv,.5,.5,.5);
    D3d_mat_mult_pt(rgb,m,rgb);  
} 

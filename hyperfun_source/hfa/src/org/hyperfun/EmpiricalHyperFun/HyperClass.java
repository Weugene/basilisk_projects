/*
 * Copyright (c) 1999 Richard Cartwright
 * All rights reserved.
 *
 * This program can be redistributed and/or modified under the terms
 * of the CGPL, The Common Good Public License as published by and at CGPL.org
 * (http://CGPL.org).  It is released under version 1.0 Beta of the License
 * until the 1.0 version is released after which either version 1.0 of the
 * License, or (at your option) any later version can be applied.
 *
 * THIS WORK, OR SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED (See the
 * CGPL, The Common Good Public License for more information.)
 *
 * You should have received a copy of the CGPL along with this program;
 * if not, see -  http://CGPL.org to get a copy of the License.
 *
 * $Id: HyperClass.java,v 1.4 2006/02/02 01:33:27 y7goto Exp $
 */

/* HyperFun Tetrahedral Decomposition Polygonalisation - Version 1.0
 * (c) Richard Cartwright - March 1999
 * 
 * Extend this class and overide method "f" with translated from
 * HyperFun code.  This application will execute this function a
 * numner of times and produce a VRML-2 rendering of the object
 * described, utilising techniques described in ...
 *
 * Bloomenthal et al.  "Introduction to implicit surfaces".  
 * Morgan-Kaufman, 1997.   
 *
 * Usage "java HyperFun.HyperClass <options>"
 *
 * For more details, use "java HyperFun.HyperClass -help"
 *
 * Note that output will be sent to the Standard Output stream 
 * (System.out) unless otherwise specified with the "-f" option.
 *
 * For all enquiries about this software, contact Richard Cartwright
 * by e-mail: R.Cartwright@dcs.warwick.ac.uk
 *
 * Version 1.0a - April 1999 
 * -  Fixed bug in index. 
 * -  Made inheritance from another package work correctly.
 * -  Added elliptical cylinders.
 * -  Added hyperbolic functions
 *
 */

/*
 * modified by PA
 * 2 7 2004
 * new primitives from HF Library
 * - cauchy arc with convolution surface
 * - cauchy curve with convolution surface
 * - caucy line with convolution surface
 * - cauchy Mesh with convolution surface
 * - cauchy point with convolution surface
 * - cauchy triangle with convolution surface
 */




/*
  Modified by PA
  2 11 2004
  modification of the space mapping primitives
  in order to match the original syntax/semantic
  like in the library
 */

package org.hyperfun.EmpiricalHyperFun;

import java.util.Date;
import java.io.*;

public class HyperClass {
    
    public static final double hfSphere(double x[], double center[],
					double R) {
	
	double x0 = x[0] - center[0];
	double x1 = x[1] - center[1];
	double x2 = x[2] - center[2];
	
	return (R*R) - (x0*x0) - (x1*x1) - (x2*x2);
    }
    
    public static final double hfEllipsoid(double x[], double center[],
					   double a, double b, double c) {
	
	double x0 = (x[0] - center[0]) / a;
	double x1 = (x[1] - center[1]) / b;
	double x2 = (x[2] - center[2]) / c;
	
	return 1 - (x0*x0) - (x1*x1) - (x2*x2);
    }
    
    public static final double hfCylX(double x[], double center[],
				      double R) {
	
	double x1 = x[1] - center[1];
	double x2 = x[2] - center[2];
	
	return (R*R) - (x1*x1) - (x2*x2);
    }
    
    public static final double hfCylY(double x[], double center[],
				      double R) {
	
	double x0 = x[0] - center[0];
	double x2 = x[2] - center[2];
	
	return (R*R) - (x0*x0) - (x2*x2);
    }
    
    public static final double hfCylZ(double x[], double center[],
				      double R) {
	
	double x0 = x[0] - center[0];
	double x1 = x[1] - center[1];
	
	return (R*R) - (x0*x0) - (x1*x1);
    }
    
    public static final double hfEllCylX(double x[], double center[],
					 double a, double b) {

	double y0 = (x[1] - center[1]) / a;
	double z0 = (x[2] - center[2]) / b;

	return 1 - (y0 * y0) - (z0 * z0); 
    }

    public static final double hfEllCylY(double x[], double center[],
					 double a, double b) {

	double x0 = (x[0] - center[0]) / a;
	double z0 = (x[2] - center[2]) / b;

	return 1 - (x0 * x0) - (z0 * z0);
    }

    public static final double hfEllCylZ(double x[], double center[],
					 double a, double b) {

	double x0 = (x[0] - center[0]) / a;
	double y0 = (x[1] - center[1]) / b;

	return 1 - (x0 * x0) - (y0 * y0);
    }

    public static final double hfTorusX(double x[], double center[],
					double R, double r0) {

	double x0 = x[0] - center[0];
	double x1 = x[1] - center[1];
	double x2 = x[2] - center[2];

	return (r0*r0) - (x0*x0) - (x1*x1) - (x2*x2) - (R*R) 
	    + 2*R*Math.sqrt((x1*x1) + (x2*x2));
    }

    public static final double hfTorusY(double x[], double center[],
					double R, double r0) {

	double x0 = x[0] - center[0];
	double x1 = x[1] - center[1];
	double x2 = x[2] - center[2];

	return (r0*r0) - (x0*x0) - (x1*x1) - (x2*x2) - (R*R) 
	    + 2*R*Math.sqrt((x0*x0) + (x2*x2));
    }

    public static final double hfTorusZ(double x[], double center[],
					double R, double r0) {

	double x0 = x[0] - center[0];
	double x1 = x[1] - center[1];
	double x2 = x[2] - center[2];

	return (r0*r0) - (x0*x0) - (x1*x1) - (x2*x2) - (R*R)
	    + 2*R*Math.sqrt((x0*x0) + (x1*x1));
    }

    public static final double hfBlock(double x[], double vertex[],
				       double dx, double dy, double dz) {

	double x0 = -(x[0] - vertex[0]) * (x[0] - (vertex[0] + dx));
	double y0 = -(x[1] - vertex[1]) * (x[1] - (vertex[1] + dy));
	double z0 = -(x[2] - vertex[2]) * (x[2] - (vertex[2] + dz));

	double i0 = x0 + y0 - Math.sqrt(x0 * x0 + y0 * y0);
	return i0 + z0 - Math.sqrt(i0 * i0 + z0 * z0);
    }

    public static final double hfBlobby(double x[],
					double x0[], double y0[], 
					double z0[], double a[], 
					double b[], double T) {

	double xx0, xx1, xx2, r;
	double sum = 0.0;

	int u;

	for ( u = 0 ; u < x0.length ; u++ )
	    { 
		xx0 = x[0] - x0[u];
		xx1 = x[1] - y0[u];
		xx2 = x[2] - z0[u];

		r = (xx0*xx0) + (xx1*xx1) + (xx2*xx2);
		sum += b[u] * Math.exp(-a[u] * r);
	    } 

	return sum - T;
    }

    public static final double hfMetaBall(double x[], double x0[],
					  double y0[], double z0[],
					  double b[], double d[],
					  double T) {

	double sum = 0.0;
	double r, xx0, xx1, xx2;

	int u;
    
	for ( u = 0 ; u < x0.length ; u++ )
	    {
		xx0 = x[0] - x0[u];
		xx1 = x[1] - y0[u];
		xx2 = x[2] - z0[u];
		r = (xx0 * xx0) + (xx1 * xx1) + (xx2 * xx2);

		if (r <= d[u])
		    {
			if (r <= (d[u] / 3))
			    sum += b[u] * (1 - 3 * r / d[u] * d[u]);
			else 
			    sum += 1.5*b[u] * (1 - r / d[u]) * (1 - r / d[u]);
		    }
	    }
		
	return sum - T;
    }

    public static final double hfSoft(double x[], double x0[],
				      double y0[], double z0[],
				      double d[], double T) {

	double sum = 0.0;
	double r, xx0, xx1, xx2, d2;

	int u;

	for ( u = 0 ; u < x0.length ; u++)
	    {
		xx0 = x[0] - x0[u];
		xx1 = x[1] - y0[u];
		xx2 = x[2] - z0[u];
		r = (xx0 * xx0) + (xx1 * xx1) + (xx2 * xx2);
		d2 = d[u]*d[u];

		sum += 1 - (22 * r) / (9 * d2) + (17 * r * r) / (9 * d2 * d2) -
		    (4 * r * r * r) / (9 * d2 * d2 * d2);
	    }

	return sum - T;
    }

    public static final double hfScale3D(double x[], double sx,
					   double sy, double sz) {

	double[] d = new double[3];
	d[0] = x[0] / sx;
	d[1] = x[1] / sy;
	d[2] = x[2] / sz;

	
	x[0]=d[0]; x[1]=d[1]; x[2]=d[2];

	return 1;
    }

    public static final double hfShift3D(double x[], double dx,
					   double dy, double dz) {

	double[] d = new double[3];
	d[0] = x[0] - dx;
	d[1] = x[1] - dy;
	d[2] = x[2] - dz;


	x[0]=d[0]; x[1]=d[1]; x[2]=d[2];

	return 1;
    }

    public static final double hfRotate3DX(double x[], double theta) {

	double[] d = new double[3];

	d[0] = x[0];
	d[1] = x[1] * Math.cos(theta) + x[2] * Math.sin(theta);
	d[2] = -x[1] * Math.sin(theta) + x[2] * Math.cos(theta);
	x[0]=d[0]; x[1]=d[1]; x[2]=d[2];
	return 1;
    }

    public static final double hfRotate3DY(double x[], double theta) {

	double[] d = new double[3];
	d[2] = x[2] * Math.cos(theta) + x[0] * Math.sin(theta);
	d[0] = -x[2] * Math.sin(theta) + x[0] * Math.cos(theta);
	d[1] = x[1];
x[0]=d[0]; x[1]=d[1]; x[2]=d[2];
	return 1;
    }

    public static final double hfRotate3DZ(double x[], double theta) {

	double[] d = new double[3];
    
	d[0] = x[0] * Math.cos(theta) + x[1] * Math.sin(theta);
	d[1] = -x[0] * Math.sin(theta) + x[1] * Math.cos(theta);
	d[2] = x[2];



x[0]=d[0]; x[1]=d[1]; x[2]=d[2];
	return 1;
    }

    public static final double hfBlendUni(double f1, double f2, double a0,
					  double a1, double a2) {

	double uni = f1 + f2 + Math.sqrt(f1 * f1 + f2 * f2);
	double f1a1 = f1 / a1;
	double f2a2 = f2 / a2;
	double disp = a0 / (1 + f1a1 * f1a1 + f2a2 * f2a2);

	return uni + disp;
    }

    public static final double hfBlendInt(double f1, double f2, double a0,
					  double a1, double a2) {

	double uni = f1 + f2 - Math.sqrt(f1 * f1 + f2 * f2);
	double f1a1 = f1 / a1;
	double f2a2 = f2 / a2;
	double disp = a0 / (1 + f1a1 * f1a1 + f2a2 * f2a2);

	return uni + disp;
    }

    public static final double hfTwistX(double x[], double x1, double x2,
					  double theta1, double theta2) {

	double t = (x[0] - x1) / (x2 - x1);
	double theta = (1 - t) * theta1 + t * theta2;
	double[] d = new double[3];
	d[0] = x[0];
	d[1] = x[1] * Math.cos(theta) + x[2] * Math.sin(theta);
	d[2] = -x[1] * Math.sin(theta) + x[2] * Math.cos(theta);
x[0]=d[0]; x[1]=d[1]; x[2]=d[2];
	return 1;
    }

    public static final double hfTwistY(double x[], double y1, double y2,
					  double theta1, double theta2) {

	double t = (x[1] - y1) / (y2 - y1);
	double theta = (1 - t) * theta1 + t * theta2;
	double[] d = new double[3];
	d[1] = x[1];
	d[2] = x[2] * Math.cos(theta) + x[0] * Math.sin(theta);
	d[0] = -x[2] * Math.sin(theta) + x[0] * Math.cos(theta);
x[0]=d[0]; x[1]=d[1]; x[2]=d[2];
	return 1;
    }

    public static final double hfTwistZ(double x[], double z1, double z2,
					  double theta1, double theta2) {

	double t = (x[2] - z1) / (z2 - z1);
	double theta = (1 - t) * theta1 + t * theta2;
	double[] d = new double[3];
	d[0] = x[0] * Math.cos(theta) + x[1] * Math.sin(theta);
	d[1] = -x[0] * Math.sin(theta) + x[1] * Math.cos(theta);
	d[2] = x[2];
x[0]=d[0]; x[1]=d[1]; x[2]=d[2];
	return 1;
    }

    public static final double hfStretch3D(double x[], double x0[], double sx,
					     double sy, double sz) {

	double[] d = new double[3];
	d[0] = x0[0] + (x[0] - x0[0]) / sz;
	d[1] = x0[1] + (x[1] - x0[1]) / sy;
	d[2] = x0[2] + (x[2] - x0[2]) / sz;
x[0]=d[0]; x[1]=d[1]; x[2]=d[2];
	return 1;
    }

    public static final double hfTaperX(double x[], double x1, double x2,
					  double s1, double s2) {

	double scale = s1;
	double[] d = new double[3];

	if (x[0] < x1)
	    scale = s1;
	if (x[0] > x2)
	    scale = s2;
	if ((x[0] >= x1) && (x[0] <= x2))
	    {
		double t = (x[0] - x1) / (x2 - x1);
		scale = (1 - t) * s1 + t * s2;
	    }
	d[0] = x[0];
	d[1] = x[1] / scale;
	d[2] = x[2] / scale;
x[0]=d[0]; x[1]=d[1]; x[2]=d[2];
	return 1;
    }

    public static final double hfTaperY(double x[], double y1, 
					  double y2, double s1,
					  double s2) {

	double scale = s1;
	double[] d = new double[3];

	if (x[1] < y1)
	    scale = s1;
	if (x[1] > y2)
	    scale = s2;
	if ((x[1] >= y1) && (x[1] <= y2))
	    {
		double t = (x[1] - y1) / (y2 - y1);
		scale = (1 - t) * s1 + t * s2;
	    }
	d[0] = x[0] / scale;
	d[1] = x[1];
	d[2] = x[2] / scale;
x[0]=d[0]; x[1]=d[1]; x[2]=d[2];
	return 1;
    }

    public static final double hfTaperZ(double x[], double z1,
					  double z2, double s1, 
					  double s2) {

	double scale = s1;
	double[] d = new double[3];

	if (x[2] < z1)
	    scale = s1;
	if (x[2] > z2)
	    scale = s2;
	if ((x[2] >= z1) && (x[2] <= z2))
	    {
		double t = (x[2] - z1) / (z2 - z1);
		scale = (1 - t) * s1 + t * s2;
	    }
	d[0] = x[0] / scale;
	d[1] = x[1] / scale;
	d[2] = x[2];
x[0]=d[0]; x[1]=d[1]; x[2]=d[2];
	return 1;
    }

    // Thanks to Ian Burnett for the hyperbolic function code 
    public final static double sinh(double x) {

	return (Math.exp(x) - Math.exp(-x)) / 2.0;
    }

    public final static double cosh(double x) {

	return (Math.exp(x) + Math.exp(-x)) / 2.0;
    }

    public final static double tanh(double x) {

	double plusX = Math.exp(x);
	double minusX = Math.exp(-x);

	return ((plusX - minusX) / (plusX + minusX));
    }
  
    // Dubious about the definition of sign at this time.
    public final static double sign(double x) {

	return -x;
    }

    public final static int index(double d) {

	return ((int) d) - 1;
    }


    /* ##### New primitives 2 7 2004 */
    private static double EPS=0.1e-5;
    private static double PI=3.141592;
    // also uses cond ? stat1 : stat2
    // maybe not in Java
    // only in C in that case needs to replace by if

    private final static double atanh(double x){
	// Bad code I don't check for Nan or infinity
	double top = (1+x);
	double bottom = (1-x);
	return 0.5*Math.log(top/bottom);
    }    

    // Cauchy Arc with Convolution Surface
    // orignal code Ken Yoshikawa
    public final static double hfConvArc(double x[], double center[], double radius[], double theta[], double axis[], double angle[], double S[], double T){
	int n;
	double r, th, rd=PI/180.0;
	double f1, f2, d2, b, a, p1, p2, p3, f=0.0;
	double i,j,k,c,s,ii,jj,kk,ij,jk,ki,is,js,ks,one_c,length;
	double cx,cy,cz,new_x,new_y,new_z;
	double tempx,tempy,tempz;
	double over_i=0.0,over_j=0.0,over_k=1.0,over_th,over_c,over_s,
	    over_one_c,over_ii,over_jj,over_kk,over_ij,over_jk,over_ki,
	    over_is,over_js,over_ks,over_x,over_y,over_z;
	int array_index,dim;
	int N = S.length;  /* the number of primitive */
	
	for(n=0;n<N;n++) {
	    cx = center[3*n];    /* Center of Arc */
	    cy = center[3*n+1];
	    cz = center[3*n+2];
	    
	    r = radius[n];
	    angle[n] += EPS;  /* avoid error */
	    
	    i = axis[3*n] + EPS; /* avoid error */
	    j = axis[3*n+1] + EPS; /* avoid error */
	    k = axis[3*n+2] + EPS; /* avoid error */
	    
	    length = Math.sqrt(i*i + j*j + k*k);
	    if( length < EPS ) {
		length = EPS;
	    }
	    
	    i /= length;   /* Calculate normal vector around which Arc rotates */
	    j /= length;
	    k /= length;
	    
	    c = Math.cos(rd * (-angle[n]));
	    s = Math.sin(rd * (-angle[n]));
	    
	    one_c = 1.0 - c;
	    
	    ii = i*i;  jj = j*j;  kk = k*k;
	    ij = i*j;  jk = j*k;  ki = k*i;
	    is = i*s;  js = j*s;  ks = k*s;
	    
	    if(theta[n] > 360.0)
		theta[n] = 360.0;
	    
	    /********** [Begin] over PI operation ***************************/
	    if(theta[n] > 180.0) {
		over_th = (theta[n] - 180.0)*rd;
		theta[n] = 180.0;
		
		/* rotate by -angle */
		tempx = (c + ii * one_c)*(x[0]-cx) + (-ks + ij * one_c)*(x[1]-cy) + (js + ki * one_c)*(x[2]-cz);
		tempy = (ks + ij * one_c)*(x[0]-cx) +  (c + jj * one_c)*(x[1]-cy) + (-is + jk * one_c)*(x[2]-cz);
		tempz = (-js + ki * one_c)*(x[0]-cx) + (is + jk * one_c)*(x[1]-cy) + (c + kk * one_c)*(x[2]-cz);
		
		/************* [Begin] rotate -PI operation **********************/
		over_c = Math.cos(rd * (-180.0));
		over_s = Math.sin(rd * (-180.0));
		over_one_c = 1.0 - over_c;
		
		over_ii = over_i*over_i; over_jj = over_j*over_j; over_kk = over_k*over_k;
		over_ij = over_i*over_j; over_jk = over_j*over_k; over_ki = over_k*over_i;
		over_is = over_i*over_s; over_js = over_j*over_s; over_ks = over_k*over_s;
		
		over_x = (over_c + over_ii * over_one_c)*(tempx) + (-over_ks + over_ij * over_one_c)*(tempy) + (over_js + over_ki * over_one_c)*(tempz);
		over_y = (over_ks + over_ij * over_one_c)*(tempx) + (over_c + over_jj * over_one_c)*(tempy) + (-over_is + over_jk * over_one_c)*(tempz);
		over_z = (-over_js + over_ki * over_one_c)*(tempx) + (over_is + over_jk * over_one_c)*(tempy) + (over_c + over_kk * over_one_c)*(tempz);
		/************* [End] rotate -PI operation **********************/
		
		a = 2.0*r*S[n]*S[n];
		d2 = over_x*over_x + over_y*over_y + over_z*over_z;
		b = 1.0 + r*r*S[n]*S[n] + S[n]*S[n]*d2;
		p2 = - r*r*r*r*S[n]*S[n]*S[n]*S[n] + 2.0*r*r*(S[n]*S[n])*((S[n]*S[n])*(d2 - 2.0*(over_z*over_z)) - 1.0) - (1.0 + (S[n]*S[n])*d2)*(1.0 + (S[n]*S[n])*d2);
		p1 = (p2 < 0.0) ? Math.sqrt(-p2) : Math.sqrt(p2);
		p3 = p1*p2;
		
		f1 = (b*over_y) / (over_x*p2*(a*over_x-b))
		    + (a*((over_x*over_x) + (over_y*over_y))*Math.sin(over_th) - b*over_y) / (over_x*p2*(a*(over_x*Math.cos(over_th) + over_y*Math.sin(over_th)) - b));
		
		if(p2 < 0.0)
		    f2 = 2.0*b*(Math.atan(-a*over_y/p1) + Math.atan((a*over_y - (a*over_x + b)*Math.tan(over_th/2.0)) / p1)) / p3;
		else
		    f2 = 2.0*b*(atanh(a*over_y/p1) + atanh(((a*over_x + b)*Math.tan(over_th/2.0) - a*over_y) / p1)) / p3;
		
		f += f1 + f2;
	    }
	    /********** [End] over PI operation ***************************/
	    
	    th = theta[n]*rd;
	    new_x = (c + ii * one_c)*(x[0]-cx) + (-ks + ij * one_c)*(x[1]-cy) + (js + ki * one_c)*(x[2]-cz);
	    new_y = (ks + ij * one_c)*(x[0]-cx) +  (c + jj * one_c)*(x[1]-cy) + (-is + jk * one_c)*(x[2]-cz);
	    new_z = (-js + ki * one_c)*(x[0]-cx) + (is + jk * one_c)*(x[1]-cy) + (c + kk * one_c)*(x[2]-cz);
	    
	    a = 2.0*r*S[n]*S[n];
	    d2 = (new_x*new_x) + (new_y*new_y) + (new_z*new_z);
	    b = 1.0 + (r*r)*(S[n]*S[n]) + (S[n]*S[n])*d2;
	    p2 = -(r*r*r*r)*(S[n]*S[n]*S[n]*S[n]) + 2.0*(r*r)*(S[n]*S[n])*((S[n]*S[n])*(d2 - 2.0*(new_z*new_z)) - 1.0) - (1.0 + (S[n]*S[n])*d2)*(1.0 + (S[n]*S[n])*d2);
	    p1 = (p2 < 0.0) ? Math.sqrt(-p2) : Math.sqrt(p2);
	    p3 = p1*p2;
	    
	    f1 = (b*new_y) / (new_x*p2*(a*new_x-b))
		+ (a*((new_x*new_x) + (new_y*new_y))*Math.sin(th) - b*new_y) / (new_x*p2*(a*(new_x*Math.cos(th) + new_y*Math.sin(th)) - b));
	    
	    if(p2 < 0.0)
		f2 = 2.0*b*(Math.atan(-a*new_y/p1) + Math.atan((a*new_y - (a*new_x + b)*Math.tan(th/2.0)) / p1)) / p3;
	    else
		f2 = 2.0*b*(atanh(a*new_y/p1) + atanh(((a*new_x + b)*Math.tan(th/2.0) - a*new_y) / p1)) / p3;
	    
	    f += f1 + f2;
	}
	
	return f - T;
    }

    public final static double hfConvCurve(double x[], double vect[], double S[], double T) {
	
	int n;
	double ax, ay, az;   /* normalized vector from beginnig to ending Point */
	double dx, dy, dz;   /* d = r - b */
	double l, p, q, xx, f=0.0;
	int N;  /* the number of primitive */

	N = S.length;

	for(n=0;n<N;n++) {
	    l = Math.sqrt((vect[3*(n+1)] - vect[3*n])*(vect[3*(n+1)] - vect[3*n]) + (vect[3*(n+1)+1] - vect[3*n+1])*(vect[3*(n+1)+1] - vect[3*n+1]) + (vect[3*(n+1)+2] - vect[3*n+2])*(vect[3*(n+1)+2] - vect[3*n+2]));
	    
	    if (l==0)
		{
		    return -1111111.0;
		}	    
	    ax = (vect[3*(n+1)] - vect[3*n]) / l;
	    ay = (vect[3*(n+1)+1] - vect[3*n+1]) / l;
	    az = (vect[3*(n+1)+2] - vect[3*n+2]) / l;
	    
	    dx = x[0] - vect[3*n];
	    dy = x[1] - vect[3*n+1];
	    dz = x[2] - vect[3*n+2];
	    
	    xx = dx*ax + dy*ay + dz*az;
	    p = Math.sqrt(1 + S[n]*S[n] * ( dx*dx + dy*dy + dz*dz - xx*xx));
	    q = Math.sqrt(1 + S[n]*S[n] * ( dx*dx + dy*dy + dz*dz + l*l - 2*l*xx ));
	    
	    f += xx / (2*p*p*(p*p + S[n]*S[n]*xx*xx)) + (l - xx) / (2*p*p*q*q)
		+ (Math.atan(S[n]*xx/p) + Math.atan(S[n]*(l - xx)/p)) / (2*S[n]*p*p*p);
	}
	
	return f - T;
    }

    

    public final static double hfConvLine(double x[], double begin[], double end[], double S[], double T) {
	int n;
	double ax, ay, az;   /* normalized vector from beginnig to ending Point */
	double dx, dy, dz;   /* d = r - b */
	double l, p, q, xx, f=0.0;
	
	int N=S.length;  /* the number of primitive */
	
	for(n=0;n<N;n++) {
	    l = Math.sqrt((end[3*n] - begin[3*n])*(end[3*n] - begin[3*n]) + (end[3*n+1] - begin[3*n+1])*(end[3*n+1] - begin[3*n+1]) + (end[3*n+2] - begin[3*n+2])*(end[3*n+2] - begin[3*n+2]));
	    
	    if(l == 0.0) {
		return -1111111.0;
	    }
	    
	    ax = (end[3*n] - begin[3*n]) / l;
	    ay = (end[3*n+1] - begin[3*n+1]) / l;
	    az = (end[3*n+2] - begin[3*n+2]) / l;
	    
	    dx = x[0] - begin[3*n];
	    dy = x[1] - begin[3*n+1];
	    dz = x[2] - begin[3*n+2];
	    
	    xx = dx*ax + dy*ay + dz*az;
	    p = Math.sqrt(1 + S[n]*S[n] * ( dx*dx + dy*dy + dz*dz - xx*xx));
	    q = Math.sqrt(1 + S[n]*S[n] * ( dx*dx + dy*dy + dz*dz + l*l - 2*l*xx ));
	    
	    f += xx / (2*p*p*(p*p + S[n]*S[n]*xx*xx)) + (l - xx) / (2*p*p*q*q)
		+ (Math.atan(S[n]*xx/p) + Math.atan(S[n]*(l - xx)/p)) / (2*S[n]*p*p*p);
	    
	}
	
	return f - T;
    }





    //#    Call      :  hfConvMesh(x,vect,tri,S,T);
    // is the cast allowed in Java: vect[(int)tri[n]] ?
    // if not can I replace double tri[] by int tri[] ?
    // Check some HF sources to see how the tri are inputed: 
    public final static double hfConvMesh(double x[], double vect[], double tri[], double S[], double T) {
	int n;
	double a1x, a1y, a1z; /* coordinate 1 */
	double a2x, a2y, a2z; /* coordinate 2 */
	double a3x, a3y, a3z; /* coordinate 3 */
	double [] len = new double[3];  /* len[0] means distance of coordinate 1 and 2 */
	/* len[1] means distance of coordinate 2 and 3 */
	/* len[2] means distance of coordinate 3 and 1 */
	double t, bx, by, bz;
	double ux, uy, uz, ul, u;  /* (ux,uy,uz) is the normal vector from b to a2 */
	double vx, vy, vz, vl, v;  /* (vx,vy,vz) is the normal vector from b to h */
	double dx, dy, dz, d2;     /* vector d = (x,y,z) - (bx,by,bz) */
	double a21x, a21y, a21z;   /* vector a21 = (a2x,a2y,a2z) - (a1x,a1y,a1z) */
	double a13x, a13y, a13z;   /* vector a1h = (a1x,a1y,a1z) - (a3x,a3y,a3z) */
	double a1, a2, h;
	double tempx, tempy, tempz;
	double A2, A, B2, B, C2, C, g, q, w, m, k, f=0.0, arc1, arc2, arc3;
	double n1, n2, n3, n4, n5, n6; 
	int N;  /* the number of primitive */
	
	N=S.length;
	
	for(n=0;n<N;n++) {
	    /* triangle coodinates */
	    a1x = vect[3*((int)tri[3*n]-1)];
	    a1y = vect[3*((int)tri[3*n]-1)+1];
	    a1z = vect[3*((int)tri[3*n]-1)+2];
	    a2x = vect[3*((int)tri[3*n+1]-1)];
	    a2y = vect[3*((int)tri[3*n+1]-1)+1];
	    a2z = vect[3*((int)tri[3*n+1]-1)+2];
	    a3x = vect[3*((int)tri[3*n+2]-1)];
	    a3y = vect[3*((int)tri[3*n+2]-1)+1];
	    a3z = vect[3*((int)tri[3*n+2]-1)+2];
	    
	    len[0] = Math.sqrt((a2x - a1x)*(a2x - a1x) + (a2y - a1y)*(a2y - a1y) + (a2z - a1z)*(a2z - a1z));
	    len[1] = Math.sqrt((a3x - a2x)*(a3x - a2x) + (a3y - a2y)*(a3y - a2y) + (a3z - a2z)*(a3z - a2z));
	    len[2] = Math.sqrt((a1x - a3x)*(a1x - a3x) + (a1y - a3y)*(a1y - a3y) + (a1z - a3z)*(a1z - a3z));
	    
	    if ((len[1] >= len[2]) && (len[1] > len[0])){
		tempx = a1x;    tempy = a1y;    tempz = a1z;
		a1x = a2x;      a1y = a2y;      a1z = a2z;
		a2x = a3x;       a2y = a3y;       a2z = a3z;
		a3x = tempx;     a3y = tempy;     a3z = tempz;
	    }else if ((len[2] >= len[1]) && (len[2] > len[0])){
		tempx = a1x;    tempy = a1y;    tempz = a1z;
		a1x = a3x;       a1y = a3y;       a1z = a3z;
		a3x = a2x;       a3y = a2y;       a3z = a2z;
		a2x = tempx;    a2y = tempy;    a2z = tempz;
	    }
	    len[0] = Math.sqrt((a2x - a1x)*(a2x - a1x) + (a2y - a1y)*(a2y - a1y) + (a2z - a1z)*(a2z - a1z));
	    len[1] = Math.sqrt((a3x - a2x)*(a3x - a2x) + (a3y - a2y)*(a3y - a2y) + (a3z - a2z)*(a3z - a2z));
	    len[2] = Math.sqrt((a1x - a3x)*(a1x - a3x) + (a1y - a3y)*(a1y - a3y) + (a1z - a3z)*(a1z - a3z));
	    
	    a21x = a2x - a1x;
	    a21y = a2y - a1y;
	    a21z = a2z - a1z;
	    a13x = a1x - a3x;
	    a13y = a1y - a3y;
	    a13z = a1z - a3z;
	    
	    t = -(a21x*a13x + a21y*a13y + a21z*a13z) / (len[0]*len[0]);
	    bx = a1x + t * a21x;
	    by = a1y + t * a21y;
	    bz = a1z + t * a21z;
	    
	    dx = x[0] - bx;
	    dy = x[1] - by;
	    dz = x[2] - bz;
	    
	    ux = a2x - bx;
	    uy = a2y - by;
	    uz = a2z - bz;
	    ul = Math.sqrt((ux*ux) + (uy*uy) + (uz*uz));
	    ux = ux / ul;
	    uy = uy / ul;
	    uz = uz / ul;
	    
	    vx = a3x - bx;
	    vy = a3y - by;
	    vz = a3z - bz;
	    vl = Math.sqrt((vx*vx) + (vy*vy) + (vz*vz));
	    vx = vx / vl;
	    vy = vy / vl;
	    vz = vz / vl;
	    
	    d2 = (dx*dx) + (dy*dy) + (dz*dz);
	    u = dx*ux + dy*uy + dz*uz;
	    v = dx*vx + dy*vy + dz*vz;
	    h = Math.sqrt( (a3x - bx)*(a3x - bx) + (a3y - by)*(a3y - by) + (a3z - bz)*(a3z - bz) );
	    a1 = Math.sqrt( (a1x - bx)*(a1x - bx) + (a1y - by)*(a1y - by) + (a1z - bz)*(a1z - bz) );
	    a2 = Math.sqrt( (a2x - bx)*(a2x - bx) + (a2y - by)*(a2y - by) + (a2z - bz)*(a2z - bz) );
	    
	    g = v - h;
	    m = a2*g + u*h;
	    k = u*h - a1*g;
	    C2 = 1.0/((S[n]*S[n])) + d2 - (u*u);
	    C = Math.sqrt(C2);
	    q = C2 - (v*v);
	    w = C2 - 2.0*v*h + (h*h);
	    A2 = (a1*a1)*w + (h*h)*(q + (u*u)) - 2.0*a1*h*u*g;
	    A = Math.sqrt(A2);
	    B2 = (a2*a2)*w + (h*h)*(q + (u*u)) + 2.0*a2*h*u*g;
	    B = Math.sqrt(B2);
	    
	    n1 = a1 + u;
	    n2 = a2 - u;
	    n3 = a1*n1 + v*h;
	    n4 = -a1*u - g*h;
	    n5 = -a2*n2 - v*h;
	    n6 = -a2*u + g*h;
	    
	    arc1 = k * (Math.atan(n3/A) + Math.atan(n4/A))/A;
	    arc2 = m * (Math.atan(n5/B) + Math.atan(n6/B))/B;
	    arc3 = v * (Math.atan(n1/C) + Math.atan(n2/C))/C;
	    f += (arc1 + arc2 + arc3)/(2.0*q*S[n]);
	}
	
	return f - T;
    }
    

    //#    Call      :  hfConvPoint(x,vect,S,T);
    public final static double hfConvPoint(double x[], double vect[], double S[], double T) {
	int n;
	double r2, f=0.0;
  
	int N=S.length;   /* the number of primitive */
	
	for(n=0;n<N;n++) {
	    r2 = (vect[3*n] - x[0])*(vect[3*n] - x[0]) + (vect[3*n+1] - x[1])*(vect[3*n +1] - x[1]) + (vect[3*n+2] - x[2])*(vect[3*n+2] - x[2]);
	    f += 1 / ((1 + (S[n]*S[n]) * r2)*(1 + (S[n]*S[n])*r2));
	}
	
	return f - T;
    }


    //#    Call      :  hfConvTriangle(x,vect,S,T);
    public final static double hfConvTriangle(double x[], double vect[], double S[], double T) {
	int n;
	double a1x, a1y, a1z; /* coordinate 1 */
	double a2x, a2y, a2z; /* coordinate 2 */
	double a3x, a3y, a3z; /* coordinate 3 */
	double[] len = new double[3];  /* len[0] means distance of coordinate 1 and 2 */
	/* len[1] means distance of coordinate 2 and 3 */
	/* len[2] means distance of coordinate 3 and 1 */
	double t, bx, by, bz;
	double ux, uy, uz, ul, u;  /* (ux,uy,uz) is the normal vector from b to a2 */
	double vx, vy, vz, vl, v;  /* (vx,vy,vz) is the normal vector from b to h */
	double dx, dy, dz, d2;     /* vector d = (x,y,z) - (bx,by,bz) */
	double a21x, a21y, a21z;   /* vector a21 = (a2x,a2y,a2z) - (a1x,a1y,a1z) */
	double a13x, a13y, a13z;   /* vector a1h = (a1x,a1y,a1z) - (a3x,a3y,a3z) */
	double a1, a2, h;
	double tempx, tempy, tempz;
	double A2, A, B2, B, C2, C, g, q, w, m, k, f=0.0, arc1, arc2, arc3;
	double n1, n2, n3, n4, n5, n6; 
	int N=S.length;  /* the number of primitive */
	
	for(n=0;n<N;n++) {
	    /* triangle coodinates */
	    a1x = vect[9*n];
	    a1y = vect[9*n+1];
	    a1z = vect[9*n+2];
	    a2x = vect[9*n+3];
	    a2y = vect[9*n+4];
	    a2z = vect[9*n+5];
	    a3x = vect[9*n+6];
	    a3y = vect[9*n+7];
	    a3z = vect[9*n+8];
	    
	    len[0] = Math.sqrt((a2x - a1x)*(a2x - a1x) + (a2y - a1y)*(a2y - a1y) + (a2z - a1z)*(a2z - a1z));
	    len[1] = Math.sqrt((a3x - a2x)*(a3x - a2x) + (a3y - a2y)*(a3y - a2y) + (a3z - a2z)*(a3z - a2z));
	    len[2] = Math.sqrt((a1x - a3x)*(a1x-a3x) + (a1y - a3y)*(a1y - a3y) + (a1z - a3z)*(a1z - a3z));
	    
	    if ((len[1] >= len[2]) && (len[1] > len[0])){
		tempx = a1x;    tempy = a1y;    tempz = a1z;
		a1x = a2x;      a1y = a2y;      a1z = a2z;
		a2x = a3x;       a2y = a3y;       a2z = a3z;
		a3x = tempx;     a3y = tempy;     a3z = tempz;
	    }else if ((len[2] >= len[1]) && (len[2] > len[0])){
		tempx = a1x;    tempy = a1y;    tempz = a1z;
		a1x = a3x;       a1y = a3y;       a1z = a3z;
		a3x = a2x;       a3y = a2y;       a3z = a2z;
		a2x = tempx;    a2y = tempy;    a2z = tempz;
	    }
	    len[0] = Math.sqrt((a2x - a1x)*(a2x - a1x) + (a2y - a1y)*(a2y - a1y) + (a2z - a1z)*(a2z - a1z));
	    len[1] = Math.sqrt((a3x - a2x)*(a3x - a2x) + (a3y - a2y)*(a3y - a2y) + (a3z - a2z)*(a3z - a2z));
	    len[2] = Math.sqrt((a1x - a3x)*(a1x -a3x) + (a1y - a3y)*(a1y - a3y) + (a1z - a3z)*(a1z -a3z));
	    
	    a21x = a2x - a1x;
	    a21y = a2y - a1y;
	    a21z = a2z - a1z;
	    a13x = a1x - a3x;
	    a13y = a1y - a3y;
	    a13z = a1z - a3z;
	    
	    t = -(a21x*a13x + a21y*a13y + a21z*a13z) / (len[0]*len[0]);
	    bx = a1x + t * a21x;
	    by = a1y + t * a21y;
	    bz = a1z + t * a21z;
	    
	    dx = x[0] - bx;
	    dy = x[1] - by;
	    dz = x[2] - bz;
	    
	    ux = a2x - bx;
	    uy = a2y - by;
	    uz = a2z - bz;
	    ul = Math.sqrt((ux*ux) + (uy*uy) + (uz*uz));
	    ux = ux / ul;
	    uy = uy / ul;
	    uz = uz / ul;
	    
	    vx = a3x - bx;
	    vy = a3y - by;
	    vz = a3z - bz;
	    vl = Math.sqrt((vx*vx) + (vy*vy) + (vz*vz));
	    vx = vx / vl;
	    vy = vy / vl;
	    vz = vz / vl;
	    
	    d2 = (dx*dx) + (dy*dy) + (dz*dz);
	    u = dx*ux + dy*uy + dz*uz;
	    v = dx*vx + dy*vy + dz*vz;
	    h = Math.sqrt( (a3x - bx)*(a3x - bx) + (a3y - by)*(a3y - by) + (a3z - bz)*(a3z - bz) );
	    a1 = Math.sqrt( (a1x - bx)*(a1x - bx) + (a1y - by)*(a1y - by) + (a1z - bz)*(a1z - bz) );
	    a2 = Math.sqrt( (a2x - bx)*(a2x - bx) + (a2y - by)*(a2y - by) + (a2z - bz)*(a2z - bz) );
	    
	    g = v - h;
	    m = a2*g + u*h;
	    k = u*h - a1*g;
	    C2 = 1.0/((S[n]*S[n])) + d2 - (u*u);
	    C = Math.sqrt(C2);
	    q = C2 - (v*v);
	    w = C2 - 2.0*v*h + (h*h);
	    A2 = (a1*a1)*w + (h*h)*(q + (u*u)) - 2.0*a1*h*u*g;
	    A = Math.sqrt(A2);
	    B2 = (a2*a2)*w + (h*h)*(q + (u*u)) + 2.0*a2*h*u*g;
	    B = Math.sqrt(B2);
	    
	    n1 = a1 + u;
	    n2 = a2 - u;
	    n3 = a1*n1 + v*h;
	    n4 = -a1*u - g*h;
	    n5 = -a2*n2 - v*h;
	    n6 = -a2*u + g*h;
	    
	    arc1 = k * (Math.atan(n3/A) + Math.atan(n4/A))/A;
	    arc2 = m * (Math.atan(n5/B) + Math.atan(n6/B))/B;
	    arc3 = v * (Math.atan(n1/C) + Math.atan(n2/C))/C;
	    f += (arc1 + arc2 + arc3)/(2.0*q*S[n]);
	}
	
	return f - T;
    }
    
    





    ///#    Call      :  hfConvLineR(x,begin,end,R);
    public final static double hfConvLineR(double x[], double begin[], double end[], double R) {
	int n;
	double T=0;        /* threshold value */
	double ax, ay, az;   /* normalized vector from beginnig to ending Point */
	double dx, dy, dz;   /* d = r - b */
	double pT,qT2;
	double l, p, q2, xx, /*S=1.0*/ S=1.0, f_tmp, f=0.0; 
	int N;  /* the number of primitive */

	if(R<1) S= 1/R;

	N = end.length;
	N /= 3;

	for(n=0;n<N;n++) {
	    l = Math.sqrt((end[3*n] - begin[3*n])*(end[3*n] - begin[3*n]) + (end[3*n+1] - begin[3*n+1])*(end[3*n + 1] - begin[3*n + 1]) + (end[3*n+2] - begin[3*n+2])*(end[3*n + 2] - begin[3*n + 2]));

	    ax = (end[3*n] - begin[3*n]) / l;
	    ay = (end[3*n+1] - begin[3*n+1]) / l;
	    az = (end[3*n+2] - begin[3*n+2]) / l;
    
	    dx = x[0] - begin[3*n];
	    dy = x[1] - begin[3*n+1];
	    dz = x[2] - begin[3*n+2];
    
	    pT = Math.sqrt(1 + (S*S)*(R*R));
	    qT2 = 1 + (S*S)*((R*R) + (l*l));

	    T = l / (2.0*(pT*pT)*qT2)
		+ (Math.atan(S*l/pT)) / (2.0*S*(pT*pT)*pT);

	    xx = dx*ax + dy*ay + dz*az;
	    p = Math.sqrt(1 + (S*S)*((dx*dx) + (dy*dy) + (dz*dz) - (xx*xx)));
	    q2 = 1 + (S*S)*((dx*dx) + (dy*dy) + (dz*dz) + (l*l) - 2.0*l*xx );

	    f_tmp = xx / (2.0*p*p*(p*p + S*S*xx*xx)) + (l - xx) / (2.0*p*p*q2)
		+ (Math.atan(S*xx/p) + Math.atan(S*(l - xx)/p)) / (2.0*S*p*p*p);

	    if(f < 1.0) {
		f = f + f_tmp;
		if(f > 1.0)
		    f = 1.0;
	    }
	}

	return f - T;
    }


    // 2005 7 15
    // 
    
    public final static double hfConeX(double x[], double center[], double R){
	double xx = x[0]-center[0];
	double yy = (x[1]-center[1])/R;
	double zz = (x[2]-center[2])/R;
	
	return xx*xx - yy*yy - zz*zz;
    }

    public final static double hfConeY(double x[], double center[], double R){
	double xx = (x[0]-center[0])/R;
	double yy = x[1]-center[1];
	double zz = (x[2]-center[2])/R;
	
	return yy*yy - xx*xx - zz*zz;
    }

    public final static double hfConeZ(double x[], double center[], double R){
	double xx = (x[0]-center[0])/R;
	double yy = (x[1]-center[1])/R;
	double zz = x[2]-center[2];

	return zz*zz - xx*xx - yy*yy;
    }
    
    public final static double hfEllConeX(double x[], double center[], double a, double b){
	double xx = x[0]-center[0];
	double yy = (x[1]-center[1])/a;
	double zz = (x[2]-center[2])/b;
	
	return xx*xx - yy*yy - zz*zz;
    }
    
    public final static double hfEllConeY(double x[], double center[], double a, double b){
	double xx = (x[0]-center[0])/a;
	double yy = x[1]-center[1];
	double zz = (x[2]-center[2])/b;
	
	return yy*yy - xx*xx - zz*zz;
    }
    
    public final static double hfEllConeZ(double x[], double center[], double a, double b){
	double xx = (x[0]-center[0])/a;
	double yy = (x[1]-center[1])/b;
	double zz = x[2]-center[2];
	
	return zz*zz - xx*xx - yy*yy;
    }

    public final static double hfFMapBlob(double x[], double fobj, double x0[], double y0[], double z0[], double fobj0[], double sigma[]){
	int num = z0.length;
	
	int i;
	double xt, yt, zt;
	double d2;
	
	for(i=0;i<num;i++){
	    xt = x[0] - x0[i];
	    yt = x[1] - y0[i];
	    zt = x[2] - z0[i];
	    if(sigma[i]>EPS){
		d2 = (xt*xt+yt*yt+zt*zt)/(sigma[i]*sigma[i]);
		fobj = fobj - fobj0[i]/(1+d2);
	    }
	}
	return fobj;
    }
    
    public static final double hfNoiseG(double x[], double amp, double freq, double phase){


	double xt, yt, zt;
	double a2x, a2y, a2z;
	double sx, sy, sz;
	double sx2, sy2, sz2;
	double serx, sery, serz;
	double a1, a2, ss;

	double a1d;

	a1=amp;
	a2=freq;

	xt=x[0];
	yt=x[1];
	zt=x[2];

	a2x=a2*xt;
	a2y=a2*yt;
	a2z=a2*zt;
	sx=Math.sin(a2x);
	sy=Math.sin(a2y);
	sz=Math.sin(a2z);
	a1d=a1/1.17;
	sx2=a1d*Math.sin(a2x/1.35+phase*sz);
	sy2=a1d*Math.sin(a2y/1.35+phase*sx);
	sz2=a1d*Math.sin(a2z/1.35+phase*sy);
	serx=a1*sx+sx2;
	sery=a1*sy+sy2;
	serz=a1*sz+sz2;
	ss=serx*sery*serz;
	return ss;
    }

    public final static double hfSpaceMapCubic(double xt[], double from[], double to[], double fratio[]){
    
	double dt[] = new double[3];
	double tmp=0.0;
	int N; 
	int i;
	double x, y, z;
	
	N = from.length;
	
	for(i=0;i<N;i+=3){
	    dt[0]=(to[i+0]-from[i+0]);
	    dt[1]=(to[i+1]-from[i+1]);
	    dt[2]=(to[i+2]-from[i+2]);
	    
	    x= (xt[0]-to[i+0])/(fratio[i+0]*(1.0+Math.sqrt(dt[0]*dt[0])));
	    y= (xt[1]-to[i+1])/(fratio[i+1]*(1.0+Math.sqrt(dt[1]*dt[1])));
	    z= (xt[2]-to[i+2])/(fratio[i+2]*(1.0+Math.sqrt(dt[2]*dt[2])));

	    x=x*x;
	    y=y*y;
	    z=z*z;
	    
	    tmp=Math.sqrt(x+y+z);
	    if(tmp<=1.0){
		tmp = ((1-tmp)*(1-tmp))*(1-tmp)/(1+tmp);
		
		xt[0]=xt[0]-tmp*(dt[0]);
		xt[1]=xt[1]-tmp*(dt[1]);
		xt[2]=xt[2]-tmp*(dt[2]);
	    }
	}
	return 1.0;
    } 

    public final static double hfSpaceMapExp(double xt[], double from[], double to[], double fratio[]){
    
	double dt[] = new double[3];
	double tmp=0.0;
	int N = from.length;
	int i;
	double x,y,z;

	for(i=0;i<N;i+=3){
	    dt[0]=(to[i+0]-from[i+0]);
	    dt[1]=(to[i+1]-from[i+1]);
	    dt[2]=(to[i+2]-from[i+2]);
	    
	    x= (xt[0]-to[i+0])/(fratio[i+0]*(1.0+Math.sqrt(dt[0]*dt[0])));
	    y= (xt[1]-to[i+1])/(fratio[i+1]*(1.0+Math.sqrt(dt[1]*dt[1])));
	    z= (xt[2]-to[i+2])/(fratio[i+2]*(1.0+Math.sqrt(dt[2]*dt[2])));
	    
	    x=x*x;
	    y=y*y;
	    z=z*z;

	    tmp = (x+y+z);
	    tmp = Math.exp(-tmp);
	    
	    xt[0]=xt[0]-tmp*dt[0];
	    xt[1]=xt[1]-tmp*dt[1];
	    xt[2]=xt[2]-tmp*dt[2];
	}
	return 1.0;
    } 

    public final static double hfSuperEll(double x[], double center[], double a, double b, double c, double s1, double s2){
	double xt, yt, zt;
	
	xt=(x[0]-center[0])/a;
	yt=(x[1]-center[1])/b;
	zt=(x[2]-center[2])/c;

	double p=2.0/s2;
	double xp=Math.pow(Math.abs(xt),p);
	double yp=Math.pow(Math.abs(yt),p);
	double zp=Math.pow(Math.abs(zt),2.0/s1);
	
	double xyp = Math.pow(xp+yp, s2/s1);
	return 1 - xyp -zp;
    }

}



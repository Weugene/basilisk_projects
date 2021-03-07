/*
 * Copyright (c) 2004 Pierre-Alain FAYOLLE
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
 * $Id: PovrayOut.java,v 1.2 2006/02/02 01:33:28 y7goto Exp $
 */

package org.hyperfun.gui;

import javax.vecmath.*;
import java.util.*;
import java.io.*;

public class PovrayOut {

    private static final int X = 0;
    private static final int Y = 1;
    private static final int Z = 2;
    Point3d[] vertices;
    Vector3f[] normals;
    int[] indices;
    StringBuffer out = new StringBuffer("");
    double min[];
    double max[];
    double lookAt[];
    double pos[];

    public PovrayOut(Point3d vertices[], Vector3f normals[], int indices[], double min[], double max[]){
	this.vertices = vertices;
	this.normals = normals;
	this.indices = indices;
	this.min = min;
	this.max = max;

	lookAt = new double[3];
	pos = new double[3];
	
	lookAt[0] = (min[0] + max[0])/2.0;
	lookAt[1] = (min[1] + max[1])/2.0;
	lookAt[2] = (min[2] + max[2])/2.0;

	pos[0] = (lookAt[0]);
	pos[1] = lookAt[1] + 2.0 * max[1];
	pos[2] = lookAt[2];
    }

    public void write(File file){
	try {
	    FileWriter fw = new FileWriter(file);
	    fw.write("#declare my_texture = texture {\n");
	    fw.write(" pigment { color rgb<0.2, 0.8, 0.2> }\n");
	    fw.write(" finish { ambient 0.2 diffuse 0.5 }\n");
	    fw.write("}\n");

	    fw.write("camera {\n");
	    fw.write(" location < " + pos[0] + ", " + 
		     pos[1] + ", " + pos[2] + " >\n");
	    
	    fw.write(" look_at < " + lookAt[0] + ", " + 
		     lookAt[1] + ", " + lookAt[2] + " >\n");
	    fw.write("}\n");
	    
	    fw.write("light_source {\n");
	    fw.write("< " + pos[0] + ", " + pos[1] + ", " + 
		     pos[2] + " >\n" + " color 1\n" + "}\n");
	    	    
	    fw.write("mesh {\n");

	    for (int i=0; i<indices.length / 3; i++){
		fw.write("smooth_triangle {\n");
		double x,y,z;
		double nx,ny,nz;

		x=vertices[indices[3*i]].x;
		y=vertices[indices[3*i]].y;
		z=vertices[indices[3*i]].z;
		
		if (Math.abs(x)<=1e-5) x=0;
		if (Math.abs(y)<=1e-5) y=0;
		if (Math.abs(z)<=1e-5) z=0;
		
		fw.write(" < " + x + ", " +
			 y + ", " + z + " >, ");
		
		x=normals[indices[3*i]].x;
		y=normals[indices[3*i]].y;
		z=normals[indices[3*i]].z;
		
		if (Math.abs(x)<=1e-5) x=0;
		if (Math.abs(y)<=1e-5) y=0;
		if (Math.abs(z)<=1e-5) z=0;
		
		fw.write("< " + x + ", " + y + 
			 ", " + z + " >\n");

		x=vertices[indices[3*i+1]].x;
		y=vertices[indices[3*i+1]].y;
		z=vertices[indices[3*i+1]].z;
		
		if (Math.abs(x)<=1e-5) x=0;
		if (Math.abs(y)<=1e-5) y=0;
		if (Math.abs(z)<=1e-5) z=0;

		fw.write(" < " + x + ", " + 
			 y + ", " + z + 
			 " >, ");
		
		x=normals[indices[3*i+1]].x;
		y=normals[indices[3*i+1]].y;
		z=normals[indices[3*i+1]].z;
		
		if (Math.abs(x)<=1e-5) x=0;
		if (Math.abs(y)<=1e-5) y=0;
		if (Math.abs(z)<=1e-5) z=0;

		fw.write("< " + x + ", " + 
			 y + ", " + 
			 z + " >\n"); 

		x=vertices[indices[3*i+2]].x;
		y=vertices[indices[3*i+2]].y;
		z=vertices[indices[3*i+2]].z;
		
		if (Math.abs(x)<=1e-5) x=0;
		if (Math.abs(y)<=1e-5) y=0;
		if (Math.abs(z)<=1e-5) z=0;

		fw.write(" < " + x + ", " + 
			 y + ", " + z + " >, ");

		x=normals[indices[3*i+2]].x;
		y=normals[indices[3*i+2]].y;
		z=normals[indices[3*i+2]].z;
		
		if (Math.abs(x)<=1e-5) x=0;
		if (Math.abs(y)<=1e-5) y=0;
		if (Math.abs(z)<=1e-5) z=0;

		fw.write("< " + x + ", " + 
			 y + ", " + z + 
			 " >\n");
		
		//fw.write(" texture { Green }\n");
		fw.write("}\n");
	    }

	    fw.write(" texture { my_texture }\n");
	    fw.write("}\n");
	    fw.close();
	}
	catch (IOException ioe) {
	    System.out.println("Error: " + ioe);
	}
    }
}

/*
 * Copyright (c) 2003 Mio Hiraga
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
 * $Id: StlOut.java,v 1.4 2006/03/24 01:06:56 pafzz Exp $
 */

package org.hyperfun.gui;
import javax.vecmath.*;
import java.util.*;
import java.io.*;

/**
 * @author MIO
 *
 * To change the template for this generated type comment go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */

/* modifications : pa, writes directly to file and not to a StringBuffer
*/

public class StlOut {
    private static final int X = 0;
    private static final int Y = 1;
    private static final int Z = 2;
    Point3d[] vertices;
    int[] indices;

    // StlOut Constructor
    public StlOut(Point3d vertices[], int indices[]) {
	this.vertices = vertices;	
	this.indices = indices;
    }

    public void write(File file) {
	double[] CB = new double[3];
	double[] CA = new double[3];
	double vec_length = 1;

	// Normal Vector for a facet
	double nx, ny, nz;

	try{
	    FileWriter fw = new FileWriter(file);

	    fw.write("solid\n");

	    for (int i = 0; i < indices.length/3; i++) {
		// Calculate a normal vector from the cross product
		CB[X] = vertices[indices[3*i + 1]].x - vertices[indices[3*i + 2]].x;
		CB[Y] = vertices[indices[3*i + 1]].y - vertices[indices[3*i + 2]].y;
		CB[Z] = vertices[indices[3*i + 1]].z - vertices[indices[3*i + 2]].z;
		CA[X] = vertices[indices[3*i]].x - vertices[indices[3*i + 2]].x;
		CA[Y] = vertices[indices[3*i]].y - vertices[indices[3*i + 2]].y;
		CA[Z] = vertices[indices[3*i]].z - vertices[indices[3*i + 2]].z;
		nx = CB[Y]*CA[Z] - CB[Z]*CA[Y];
		ny = CB[Z]*CA[X] - CB[X]*CA[Z];
		nz = CB[X]*CA[Y] - CB[Y]*CA[X];
		
		// Normalize the calculated normal vector
		vec_length = Math.sqrt(nx*nx + ny*ny + nz*nz);
		nx = nx / vec_length;
		ny = ny / vec_length;
		nz = nz / vec_length;
		
		fw.write(" facet normal " + (float)nx + " " + (float)ny + " " + (float)nz + "\n");
		fw.write("  outer loop" + "\n"); 
		fw.write("   vertex ");
		fw.write((float)vertices[indices[3*i]].x + " ");
		fw.write((float)vertices[indices[3*i]].y + " ");
		fw.write((float)vertices[indices[3*i]].z + "\n");
		
		fw.write("   vertex ");
		fw.write((float)vertices[indices[3*i + 1]].x + " ");
		fw.write((float)vertices[indices[3*i + 1]].y + " ");
		fw.write((float)vertices[indices[3*i + 1]].z + "\n");
		
		fw.write("   vertex ");
		fw.write((float)vertices[indices[3*i + 2]].x + " ");
		fw.write((float)vertices[indices[3*i + 2]].y + " ");
		fw.write((float)vertices[indices[3*i + 2]].z + "\n");
		
		fw.write("  endloop" + "\n"); 
		fw.write(" endfacet" + "\n"); 
	    }
	    
	    fw.write("endsolid\n");
	    fw.close();
	} catch(IOException ioe) {
	    System.out.println("Error: " + ioe);
	}
    }
}

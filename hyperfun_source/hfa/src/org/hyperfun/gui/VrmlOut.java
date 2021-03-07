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
 * $Id: VrmlOut.java,v 1.2 2006/02/02 01:33:28 y7goto Exp $
 */

package org.hyperfun.gui;

import java.awt.Color;
import javax.vecmath.*;

/**
 * @author MIO
 *
 * To change the template for this generated type comment go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
public class VrmlOut {
    // HyperFun code from TextArea
    String buf;

    // Object Color
    Color materialColor;

    Point3d[] vertices;
    Vector3f[] normals;
    int[] indices;

    // CameraPosition(0, 0, -15)
    int[] itsCPos = {0, 0, 15};

    // For outtting vrml data converter from HFcode
    StringBuffer outputtxt = new StringBuffer("");

    public VrmlOut(String buf, Color materialColor, Point3d[] vertices, Vector3f[] normals, int[] indices) {
	this.buf = buf;
	this.materialColor = materialColor;
	this.vertices = vertices;	
	this.normals = normals;
	this.indices = indices;
	//vrmlout();
    }

    public String vrmlout() {
	//String itsModelName = getModelName(buf);
	String itsModelName = "my_model";

	// COMMENTS
	outputtxt.append("#VRML V2.0 utf8\n");
	outputtxt.append("\n" );
	outputtxt.append("#Produced by HyperFun Applet Viewer\n");
	outputtxt.append("\n");
	outputtxt.append("#Object: " + itsModelName + "\n");
	outputtxt.append("\n");

	// CAMERA
	outputtxt.append("DEF Start Viewpoint {\n");
	outputtxt.append("  position " + itsCPos[0]+ " " + itsCPos[1] + " " + 2*itsCPos[2] + "\n");
	outputtxt.append("  orientation " + -itsCPos[0] + " " + -itsCPos[1] + " " + -itsCPos[2] + " 0\n");
	outputtxt.append("  fieldOfView " + (60/360)*Math.PI + "\n");
	outputtxt.append("  description \"Start\"\n");
	outputtxt.append("  }\n");

	// OBJECT
	outputtxt.append("DEF " + itsModelName + " Transform {\n");
	outputtxt.append("  translation 0 0 0\n");
	outputtxt.append("  children [\n");
	outputtxt.append("    Shape {\n");
	outputtxt.append("      appearance Appearance {\n");
	outputtxt.append("        material Material {\n");
	outputtxt.append("          diffuseColor " +(float)materialColor.getRed()/255 + " " + (float)materialColor.getGreen()/255 + " " + (double)materialColor.getBlue()/255 + "\n");
	outputtxt.append("          specularColor 1 1 1\n");
	outputtxt.append("          shininess 0.64\n");
	outputtxt.append("          }\n");
	outputtxt.append("        }\n");
	outputtxt.append("      geometry DEF " + itsModelName + "-FACES IndexedFaceSet {\n");
	outputtxt.append("        ccw TRUE\n");
	outputtxt.append("        solid FALSE\n");

	// Vertices
	outputtxt.append("        coord DEF " + itsModelName + "-COORD Coordinate {\n");
	outputtxt.append("          point [\n");
	for(int A=0; A<vertices.length; A++){
	outputtxt.append("            ");
	outputtxt.append((float)vertices[A].x + " ");
	outputtxt.append((float)vertices[A].y + " ");
	outputtxt.append((float)vertices[A].z + "\n");
	}
	outputtxt.append("            ]\n");
	outputtxt.append("          }\n");

	// Normals
	outputtxt.append("        normal Normal {\n");
	outputtxt.append("          vector [\n");
	for (int B = 0; B < normals.length; B++) {
	    outputtxt.append("            ");
	    outputtxt.append(normals[B].x + " ");
	    outputtxt.append(normals[B].y + " ");
	    outputtxt.append(normals[B].z + "\n");
	}
	outputtxt.append("            ]\n");
	outputtxt.append("          }\n");
	outputtxt.append("        normalPerVertex TRUE\n");

	// Vert Index
	outputtxt.append("        coordIndex [\n");
	for (int C = 0; C < indices.length/3; C++) {
	    outputtxt.append("          ");
	    outputtxt.append(indices[3*C] + ", ");
	    outputtxt.append(indices[3*C+1] + ", ");
	    outputtxt.append(indices[3*C+2] + ", -1, \n");
	}
	outputtxt.append("          ]\n");

	//Norm Index
	outputtxt.append("        normalIndex [\n");
	for(int D = 0; D < indices.length/3; D++) {
	    outputtxt.append("          ");
	    outputtxt.append(indices[3*D] + ", ");
	    outputtxt.append(indices[3*D+1] + ", ");
	    outputtxt.append(indices[3*D+2] + ", -1, \n");
	}
	outputtxt.append("          ]\n");

	//Color Attribute
	/*	pw.print( "        colorPerVertex TRUE\n");
	pw.print( "        color Color { color [\n");//1 0 0, 0 1 0, 0 0 1, 1,1,0, 0,1,1]}
	for(int E=0); E<itsPolyMesh_->VertexNum()); E++){
	pw.print( "          ");
	pw.print( itsPolyMesh_->itsData.getAttributes(E, 0) + ", ");
	pw.print( itsPolyMesh_->itsData.getAttributes(E, 1) + ", ");
	pw.print( itsPolyMesh_->itsData.getAttributes(E, 2) + ", \n");
	}
	pw.print( "          ]}\n");
	*/	

	outputtxt.append( "        }\n");
	outputtxt.append( "      }\n");
	outputtxt.append( "    ]\n");
	outputtxt.append( "  }\n");

	return outputtxt.toString();
    }

    public String getModelName(String buf) {
	char c_hfcode[] = buf.toCharArray();
	char dummy[] = new char[36];
	int cnt=0;

	for (int i = 0; i < buf.length(); i++ ) {
	    if (c_hfcode[i] == '(' ||  c_hfcode[i] == ' ') {
		break;
	    } else {
		dummy[i] = c_hfcode[i];
		cnt++;
	    }
	}

	char modelName[] = new char[cnt];
	for (int i = 0; i < cnt; i++) {
	    modelName[i] = dummy[i];
	}
	String itsModelName = new String(modelName);

	return itsModelName;
    }
}

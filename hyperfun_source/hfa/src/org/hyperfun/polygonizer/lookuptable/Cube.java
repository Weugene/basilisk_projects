/*
 * Copyright (c) 2006 HyperFun Project
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
 * $Id: Cube.java,v 1.1 2006/11/08 09:13:15 y7goto Exp $
 */

package org.hyperfun.polygonizer.lookuptable;

/**
 * @author Yuichiro Goto
 * @version $Revision: 1.1 $
 */
public final class Cube {
    private int index;
    private Edge[] edges;
    private Face[] faces;

    Cube(int index) {
	if (index < 0 || index > 255) {
	    throw new IllegalArgumentException(
		"index must be in the range between 0 and 255");
	}
	this.index = index;

	edges = new Edge[12];
	for (int edgeIndex = 0; edgeIndex < 12; edgeIndex++) {
	    edges[edgeIndex] = new Edge(edgeIndex);
	}
	faces = new Face[6];
	for (int faceIndex = 0; faceIndex < 6; faceIndex++) {
	    faces[faceIndex] = FaceFactory.createFace(faceIndex, index, edges);
	}
    }

    public int getIndex() {
	return index;
    }

    public Edge getEdge(int index) {
	return edges[index];
    }

    public int getEdgeCount() {
	return edges.length;
    }

    public Face getFace(int index) {
	return faces[index];
    }

    public int getFaceCount() {
	return faces.length;
    }

    private static String indexToString(int index) {
	return (((index & (1<<7)) != 0) ? "1" : "0") +
	       (((index & (1<<6)) != 0) ? "1" : "0") +
	       (((index & (1<<5)) != 0) ? "1" : "0") +
	       (((index & (1<<4)) != 0) ? "1" : "0") +
	       (((index & (1<<3)) != 0) ? "1" : "0") +
	       (((index & (1<<2)) != 0) ? "1" : "0") +
	       (((index & (1<<1)) != 0) ? "1" : "0") +
	       (((index & (1<<0)) != 0) ? "1" : "0");
    }

    public int[] getEdgeConnectivity(int[] connectionSwitches) throws FaceNotAmbiguousException {
	int[] connectivity = {
	    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
	};
	for (int faceIndex = 0; faceIndex < 6; faceIndex++) {
	    Face face = faces[faceIndex];
	    if (face.isAmbiguous() == false && connectionSwitches[faceIndex] != 0) {
		throw new FaceNotAmbiguousException(faceIndex +
		    ((faceIndex == 1) ? "st" :
		    ((faceIndex == 2) ? "nd" :
		    ((faceIndex == 3) ? "rd" : "th"))) + " face of the cube " +
		    index + " (" + indexToString(index) + ") is not ambiguous");
	    }
	    for (int edgeIndex = 0; edgeIndex < 4; edgeIndex++) {
		Edge edge = face.getEdge(edgeIndex);
		if (edge.getConnectedEdge(0) != null && face.contains(edge.getConnectedEdge(0))) {
		    connectivity[edge.getIndex()] =
			edge.getConnectedEdge(connectionSwitches[faceIndex]).getIndex();
		}
	    }
	}
	return connectivity;
    }

    public String toString() {
	return "Cube" + index + "[" + faces[0] + "," + faces[1] + "," + faces[2] + "," +
				      faces[3] + "," + faces[4] + "," + faces[5] + "]";
    }
}

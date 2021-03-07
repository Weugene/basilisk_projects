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
 * $Id: FaceFactory.java,v 1.1 2006/11/08 09:13:15 y7goto Exp $
 */

package org.hyperfun.polygonizer.lookuptable;

/**
 * @author Yuichiro Goto
 * @version $Revision: 1.1 $
 */
class FaceFactory {
    private static final int[][] FACE_VERTICES = {
	{0, 1, 2, 3},
	{0, 1, 5, 4},
	{0, 3, 7, 4},
	{4, 5, 6, 7},
	{3, 2, 6, 7},
	{1, 2, 6, 5}
    };

    private static final int[][] FACE_EDGES = {
	{0,  1,  2,  3},
	{0,  9,  4,  8},
	{3, 11,  7,  8},
	{4,  5,  6,  7},
	{2, 10,  6, 11},
	{1, 10,  5,  9}
    };

    private static final int[][][] EDGE_CONNECTIVITY_ON_FACE = {
	{{-1,-1,-1,-1}, null},
	{{-1,-1,-1, 0}, null},
	{{ 1,-1,-1,-1}, null},
	{{-1,-1,-1, 1}, null},
	{{-1, 2,-1,-1}, null},
	{{-1, 0,-1, 2}, {-1, 2,-1, 0}},
	{{ 2,-1,-1,-1}, null},
	{{-1,-1,-1, 2}, null},
	{{-1,-1, 3,-1}, null},
	{{-1,-1, 0,-1}, null},
	{{ 1,-1, 3,-1}, { 3,-1, 1,-1}},
	{{-1,-1, 1,-1}, null},
	{{-1, 3,-1,-1}, null},
	{{-1, 0,-1,-1}, null},
	{{ 3,-1,-1,-1}, null},
	{{-1,-1,-1,-1}, null} 
    };

    // Clockwise
    private static final int CW  = 1;
    // Counter clockwise
    private static final int CCW = 0;

    private static final int[] FACE_ORIENTATION = {CW, CCW, CW, CCW, CW, CCW};

    private static boolean isAmbiguousBitPattern(int bitPatternOnFace) {
	return (bitPatternOnFace == 5 || bitPatternOnFace == 10) ? true : false;
    }

    private static boolean isBitOn(int bitPatternOnCube, int vertexIndex) {
	return ((bitPatternOnCube & (1 << vertexIndex)) != 0) ? true : false;
    }

    private static int buildBitPatternOnFace(int bitPatternOnCube, int faceIndex) {
	int bitPatternOnFace = 0;
	for (int vertexIndex = 0; vertexIndex < 4; vertexIndex++) {
	    if (isBitOn(bitPatternOnCube, FACE_VERTICES[faceIndex][vertexIndex])) {
		bitPatternOnFace |= 1 << vertexIndex;
	    }
	}
	return bitPatternOnFace;
    }

    static Face createFace(int faceIndex, int bitPatternOnCube, Edge[] edges) {
	if (faceIndex < 0 || faceIndex > 5) {
	    throw new IllegalArgumentException(
		"faceIndex must be in the range between 0 and 5");
	}
	if (bitPatternOnCube < 0 || bitPatternOnCube > 255) {
	    throw new IllegalArgumentException(
		"bitPatternOnCube must be in the range between 0 and 255");
	}
	if (edges.length != 12) {
	    throw new IllegalArgumentException(
		"length of edges must be 12");
	}
	int bitPatternOnFace = buildBitPatternOnFace(bitPatternOnCube, faceIndex);

	Face face = new Face(faceIndex, new Edge[] {
			     edges[FACE_EDGES[faceIndex][0]],
			     edges[FACE_EDGES[faceIndex][1]],
			     edges[FACE_EDGES[faceIndex][2]],
			     edges[FACE_EDGES[faceIndex][3]]},
			     isAmbiguousBitPattern(bitPatternOnFace));

	int[][] connectivity = EDGE_CONNECTIVITY_ON_FACE[bitPatternOnFace];
	for (int i = 0; i < 2; i++) {
	    if (connectivity[i] != null) {
		for (int vertexIndex = 0; vertexIndex < 4; vertexIndex++) {
		    if (connectivity[i][vertexIndex] != -1) {
			if (FACE_ORIENTATION[faceIndex] == CW) {
			    Edge edge = face.getEdge(vertexIndex);
			    edge.setConnectedEdge(i, face.getEdge(connectivity[i][vertexIndex]));
			} else {
			    Edge edge = face.getEdge(connectivity[i][vertexIndex]);
			    edge.setConnectedEdge(i, face.getEdge(vertexIndex));
			}
		    }
		}
	    }
	}
	return face;
    }

    /** Don't allow instantiating this class */
    private FaceFactory() {
    }
}

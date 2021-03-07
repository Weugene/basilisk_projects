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
 * $Id: Polygonizer.java,v 1.7 2006/11/28 12:41:41 y7goto Exp $
 */

package org.hyperfun.polygonizer;

import java.util.Collection;
import java.util.Hashtable;

import org.hyperfun.polygonizer.lookuptable.Cube;
import org.hyperfun.polygonizer.lookuptable.Edge;
import org.hyperfun.polygonizer.lookuptable.Face;
import org.hyperfun.polygonizer.lookuptable.FaceNotAmbiguousException;
import org.hyperfun.polygonizer.lookuptable.LookupTable;

/**
 * Isosurface polygonizer based on the following thesis and the Marching Cubes
 * algorithm:
 * <pre>
 * Pasko A.A., Pilyugin V.V., Pokrovskiy V.N. "Geometric modeling in the
 * analysis of trivariate functions", Communications of Joint Insititute of
 * Nuclear Research, P10-86-310, Dubna, Russia, 1986 (in Russian).
 * </pre>
 *
 * @author Yuichiro Goto
 * @version $Revision: 1.7 $
 */
public final class Polygonizer {
    // Lower left front corner of the bounding box
    private double xMin, yMin, zMin;
    // Upper right back corner of the bounding box
    private double xMax, yMax, zMax;
    // Number of divisions along each axis of the bounding box
    private int xDiv, yDiv, zDiv;
    // Isovalue of the isosurface
    private double isovalue;
    // Function defining the isosurface
    private Function function;
    private double dx, dy, dz;
    private double ndx, ndy, ndz;

    public Polygonizer() {
	this(new double[] {-1.0, -1.0, -1.0}, new double[] {1.0, 1.0, 1.0},
	     new int[] {1, 1, 1}, 0.0, null);
    }

    public Polygonizer(Function function) {
	this(new double[] {-1.0, -1.0, -1.0}, new double[] {1.0, 1.0, 1.0},
	     new int[] {1, 1, 1}, 0.0, function);
    }

    public Polygonizer(double[] min, double[] max,
		       int[] div, double isovalue, Function function) {
	setBoundingBox(min, max);
	setDivisions(div);
	setIsovalue(isovalue);
	setFunction(function);
    }

    public void getBoundingBox(double[] min, double[] max) {
	min[0] = xMin; min[1] = yMin; min[2] = zMin;
	max[0] = xMax; max[1] = yMax; max[2] = zMax;
    }

    public void setBoundingBox(double xMin, double yMin, double zMin,
			       double xMax, double yMax, double zMax) {
	if (xMin < xMax) {
	    this.xMin = xMin;
	    this.xMax = xMax;
	} else if (xMin > xMax) {
	    this.xMin = xMax;
	    this.xMax = xMin;
	} else {
	    this.xMin = -1.0;
	    this.xMax =  1.0;
	}
	if (yMin < yMax) {
	    this.yMin = yMin;
	    this.yMax = yMax;
	} else if (yMin > yMax) {
	    this.yMin = yMax;
	    this.yMax = yMin;
	} else {
	    this.yMin = -1.0;
	    this.yMax =  1.0;
	}
	if (zMin < zMax) {
	    this.zMin = zMin;
	    this.zMax = zMax;
	} else if (zMin > zMax) {
	    this.zMin = zMax;
	    this.zMax = zMin;
	} else {
	    this.zMin = -1.0;
	    this.zMax =  1.0;
	}
    }

    public void setBoundingBox(double min[], double max[]) {
	setBoundingBox(min[0], min[1], min[2],
		       max[0], max[1], max[2]);
    }

    public void getDivisions(int[] div) {
	div[0] = xDiv;
	div[1] = yDiv;
	div[2] = zDiv;
    }

    public void setDivisions(int xDiv, int yDiv, int zDiv) {
	this.xDiv = (xDiv > 0) ? xDiv : 1;
	this.yDiv = (yDiv > 0) ? yDiv : 1;
	this.zDiv = (zDiv > 0) ? zDiv : 1;
    }

    public void setDivisions(int div[]) {
	setDivisions(div[0], div[1], div[2]);
    }

    public double getIsovalue() {
	return isovalue;
    }

    public void setIsovalue(double isovalue) {
	this.isovalue = isovalue;
    }

    public Function getFunction() {
	return function;
    }

    public void setFunction(Function function) {
	this.function = function;
    }

    public void polygonize(Collection vertices,
			   Collection normals,
			   Collection indices) {
	if (function == null) {
	    return;
	}
	double[] values = new double[8];
	double[][] positionsD = new double[8][3];
	int[][] positionsI = new int[8][3];
	int[] connectionSwitches = new int[6];
	int[] edgeToIndex = new int[12];
	Hashtable indexTable = new Hashtable();
	double[][] upperPlane = new double[yDiv+1][xDiv+1];
	double[][] lowerPlane = new double[yDiv+1][xDiv+1];
	double eps = (isovalue == 0.0) ? 1.0E-5 : isovalue * 1.0E-5;

	dx = (xMax - xMin) / xDiv;
	dy = (yMax - yMin) / yDiv;
	dz = (zMax - zMin) / zDiv;

	ndx = 0.001 * dx;
	ndy = 0.001 * dy;
	ndz = 0.001 * dz;

	sample(lowerPlane, zMin);
	for (int k = 0; k < zDiv; k++) {
	    double z1 = zMin + k * dz;
	    double z2 = zMin + (k + 1) * dz;
	    sample(upperPlane, z2);
	    for (int j = 0; j < yDiv; j++) {
		double y1 = yMin + j * dy;
		double y2 = yMin + (j + 1) * dy;
		for (int i = 0; i < xDiv; i++) {
		    double x1 = xMin + i * dx;
		    double x2 = xMin + (i + 1) * dx;
		    // Set sampled function values on each corner of the
		    // cube
		    values[0] = lowerPlane[j  ][i  ];
		    values[1] = lowerPlane[j+1][i  ];
		    values[2] = lowerPlane[j+1][i+1];
		    values[3] = lowerPlane[j  ][i+1];
		    values[4] = upperPlane[j  ][i  ];
		    values[5] = upperPlane[j+1][i  ];
		    values[6] = upperPlane[j+1][i+1];
		    values[7] = upperPlane[j  ][i+1];

		    // Adjust the function values which are almost same as the
		    // isovalue
		    if (Math.abs(values[0] - isovalue) < eps) {values[0] += 10.0 * eps;}
		    if (Math.abs(values[1] - isovalue) < eps) {values[1] += 10.0 * eps;}
		    if (Math.abs(values[2] - isovalue) < eps) {values[2] += 10.0 * eps;}
		    if (Math.abs(values[3] - isovalue) < eps) {values[3] += 10.0 * eps;}
		    if (Math.abs(values[4] - isovalue) < eps) {values[4] += 10.0 * eps;}
		    if (Math.abs(values[5] - isovalue) < eps) {values[5] += 10.0 * eps;}
		    if (Math.abs(values[6] - isovalue) < eps) {values[6] += 10.0 * eps;}
		    if (Math.abs(values[7] - isovalue) < eps) {values[7] += 10.0 * eps;}

		    // Calculate index into the lookup table
		    int cubeIndex = 0;
		    if (values[0] > isovalue) {cubeIndex +=   1;}
		    if (values[1] > isovalue) {cubeIndex +=   2;}
		    if (values[2] > isovalue) {cubeIndex +=   4;}
		    if (values[3] > isovalue) {cubeIndex +=   8;}
		    if (values[4] > isovalue) {cubeIndex +=  16;}
		    if (values[5] > isovalue) {cubeIndex +=  32;}
		    if (values[6] > isovalue) {cubeIndex +=  64;}
		    if (values[7] > isovalue) {cubeIndex += 128;}

		    // Skip the empty cube
		    if (cubeIndex == 0 || cubeIndex == 255) {
			continue;
		    }
		    Cube cube = LookupTable.getCube(cubeIndex);

		    // Set up corner positions of the cube
		    positionsD[0][0] = x1;  positionsD[0][1] = y1;  positionsD[0][2] = z1;
		    positionsD[1][0] = x1;  positionsD[1][1] = y2;  positionsD[1][2] = z1;
		    positionsD[2][0] = x2;  positionsD[2][1] = y2;  positionsD[2][2] = z1;
		    positionsD[3][0] = x2;  positionsD[3][1] = y1;  positionsD[3][2] = z1;
		    positionsD[4][0] = x1;  positionsD[4][1] = y1;  positionsD[4][2] = z2;
		    positionsD[5][0] = x1;  positionsD[5][1] = y2;  positionsD[5][2] = z2;
		    positionsD[6][0] = x2;  positionsD[6][1] = y2;  positionsD[6][2] = z2;
		    positionsD[7][0] = x2;  positionsD[7][1] = y1;  positionsD[7][2] = z2;

		    positionsI[0][0] = i;   positionsI[0][1] = j;   positionsI[0][2] = k;
		    positionsI[1][0] = i;   positionsI[1][1] = j+1; positionsI[1][2] = k;
		    positionsI[2][0] = i+1; positionsI[2][1] = j+1; positionsI[2][2] = k;
		    positionsI[3][0] = i+1; positionsI[3][1] = j;   positionsI[3][2] = k;
		    positionsI[4][0] = i;   positionsI[4][1] = j;   positionsI[4][2] = k+1;
		    positionsI[5][0] = i;   positionsI[5][1] = j+1; positionsI[5][2] = k+1;
		    positionsI[6][0] = i+1; positionsI[6][1] = j+1; positionsI[6][2] = k+1;
		    positionsI[7][0] = i+1; positionsI[7][1] = j;   positionsI[7][2] = k+1;

		    // Find the cube edges which have intersection points with
		    // the isosurface
		    for (int edgeIndex = 0; edgeIndex < 12; edgeIndex++) {
			Edge edge = cube.getEdge(edgeIndex);
			if (edge.getConnectedEdge(0) != null) {
			    EdgeKey key = new EdgeKey(positionsI[edge.getStartVertexIndex()],
						      positionsI[edge.getEndVertexIndex()]);
			    if (indexTable.containsKey(key)) {
				edgeToIndex[edgeIndex] = ((Integer)indexTable.get(key)).intValue();
			    } else {
				double t = (isovalue - values[edge.getStartVertexIndex()]) /
					   (values[edge.getEndVertexIndex()] - values[edge.getStartVertexIndex()]);
				double[] v = lerp(t, positionsD[edge.getStartVertexIndex()],
						     positionsD[edge.getEndVertexIndex()]);
				vertices.add(v);
				if (normals != null) {
				    normals.add(calcNormal(v));
				}
				indexTable.put(key, new Integer(edgeToIndex[edgeIndex] = vertices.size()-1));
			    }
			}
		    }

		    // Resolve topological ambiguity on cube faces
		    for (int faceIndex = 0; faceIndex < 6; faceIndex++) {
			Face face = cube.getFace(faceIndex);
			if (face.isAmbiguous()) {
			    double d0 = values[face.getEdge(0).getEndVertexIndex()] -
					values[face.getEdge(0).getStartVertexIndex()];
			    double d1 = values[face.getEdge(2).getEndVertexIndex()] -
					values[face.getEdge(2).getStartVertexIndex()];
			    double t = (isovalue - values[face.getEdge(1).getStartVertexIndex()]) /
				       (values[face.getEdge(1).getEndVertexIndex()] -
					values[face.getEdge(1).getStartVertexIndex()]);
			    connectionSwitches[faceIndex] = (t > -d0 / (d1 - d0)) ? 1 : 0;
			} else {
			    connectionSwitches[faceIndex] = 0;
			}
		    }

		    // Get the connectivity graph of the cube edges and trace
		    // it to generate triangles
		    int[] connectivity;
		    try {
			connectivity = cube.getEdgeConnectivity(connectionSwitches);
		    } catch (FaceNotAmbiguousException e) {
			// This should not happen
			throw new RuntimeException(e);
		    }
		    for (int edgeIndex = 0; edgeIndex < 12;) {
			if (connectivity[edgeIndex] != -1) {
			    int index0 = edgeIndex;
			    int index1 = connectivity[index0];
			    int index2 = connectivity[index1];

			    indices.add(new Integer(edgeToIndex[index0]));
			    indices.add(new Integer(edgeToIndex[index1]));
			    indices.add(new Integer(edgeToIndex[index2]));

			    connectivity[index0] = -1;
			    connectivity[index1] = -1;
			    if (connectivity[index2] != index0) {
				connectivity[index0] = index2;
				continue;
			    }
			    connectivity[index2] = -1;
			}
			edgeIndex++;
		    }
		}
	    }
	    // Swap the lower and upper plane
	    double[][] tmp = lowerPlane;
	    lowerPlane = upperPlane;
	    upperPlane = tmp;
	}
    }

    private void sample(double[][] plane, double z) {
	for (int j = 0; j <= yDiv; j++) {
	    double y = yMin + j * dy;
	    for (int i = 0; i <= xDiv; i++) {
		double x = xMin + i * dx;
		plane[j][i] = function.evaluate(x, y, z);
	    }
	}
    }

    private static double[] lerp(double t, double[] v0, double[] v1) {
	return new double[] {v0[0] + t * (v1[0] - v0[0]),
			     v0[1] + t * (v1[1] - v0[1]),
			     v0[2] + t * (v1[2] - v0[2])};
    }

    private float[] calcNormal(double[] v) {
	double x = v[0];
	double y = v[1];
	double z = v[2];

	double f = function.evaluate(x, y, z);
	double nx = -(function.evaluate(x+ndx, y, z) - f) / ndx;
	double ny = -(function.evaluate(x, y+ndy, z) - f) / ndy;
	double nz = -(function.evaluate(x, y, z+ndz) - f) / ndz;

	double len = Math.sqrt(nx*nx + ny*ny + nz*nz);
	if (len > 0.0) {
	    nx /= len;
	    ny /= len;
	    nz /= len;
	}
	return new float[] {(float)nx, (float)ny, (float)nz};
    }

    /**
     * EdgeKey class is based on the code found in Jules Bloomenthal's
     * implicit.c:
     * <pre>
     * implicit.c
     *     an implicit surface polygonizer, translated from Mesa
     *     applications should call polygonize()
     *
     * to compile a test program for ASCII output:
     *     cc implicit.c -o implicit -lm
     *
     * to compile a test program for display on an SGI workstation:
     *     cc -DSGIGFX implicit.c -o implicit -lgl_s -lm
     *
     * Authored by Jules Bloomenthal, Xerox PARC.
     * Copyright (c) Xerox Corporation, 1991.  All rights reserved.
     * Permission is granted to reproduce, use and distribute this code for
     * any and all purposes, provided that this notice appears in all copies.
     * </pre>
     */
    private static final class EdgeKey {
	private static final int BIT_SHIFT = 10;
	private static final int BIT_MASK = (1<<BIT_SHIFT)-1;

	private int i0, j0, k0;
	private int i1, j1, k1;

	public EdgeKey(int[] p0, int[] p1) {
	    this(p0[0], p0[1], p0[2],
		 p1[0], p1[1], p1[2]);
	}
	public EdgeKey(int i0, int j0, int k0,
		       int i1, int j1, int k1) {
	    if (i0 < i1 || (i0 == i1 && (j0 < j1 || (j0 == j1 && k0 < k1)))) {
		this.i0 = i0; this.j0 = j0; this.k0 = k0;
		this.i1 = i1; this.j1 = j1; this.k1 = k1;
	    } else {
		this.i0 = i1; this.j0 = j1; this.k0 = k1;
		this.i1 = i0; this.j1 = j0; this.k1 = k0;
	    }
	}
	public boolean equals(Object obj) {
	    if (this == obj) {
		return true;
	    }
	    if (obj instanceof EdgeKey) {
		EdgeKey key = (EdgeKey)obj;
		if (i0 == key.i0 && j0 == key.j0 && k0 == key.k0 &&
		    i1 == key.i1 && j1 == key.j1 && k1 == key.k1) {
		    return true;
		}
	    }
	    return false;
	}
	public int hashCode() {
	    return (((((i0&BIT_MASK)<<BIT_SHIFT)|(j0&BIT_MASK))<<BIT_SHIFT)|(k0&BIT_MASK)) +
		   (((((i1&BIT_MASK)<<BIT_SHIFT)|(j1&BIT_MASK))<<BIT_SHIFT)|(k1&BIT_MASK));
	}
    }
}

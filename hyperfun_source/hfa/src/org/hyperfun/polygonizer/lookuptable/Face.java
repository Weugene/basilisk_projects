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
 * $Id: Face.java,v 1.1 2006/11/08 09:13:15 y7goto Exp $
 */

package org.hyperfun.polygonizer.lookuptable;

/**
 * @author Yuichiro Goto
 * @version $Revision: 1.1 $
 */
public final class Face {
    private int index;
    private Edge[] edges;
    private boolean ambiguous;

    Face(int index, Edge[] edges, boolean ambiguous) {
	if (index < 0 || index > 5) {
	    throw new IllegalArgumentException(
		"index must be in the range between 0 and 5");
	}
	this.index = index;

	if (edges.length != 4) {
	    throw new IllegalArgumentException(
		"length of edges must be 4");
	}
	this.edges = edges;

	this.ambiguous = ambiguous;
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

    public boolean isAmbiguous() {
	return ambiguous;
    }

    public boolean contains(Edge edge) {
	return (edge == edges[0] || edge == edges[1] ||
		edge == edges[2] || edge == edges[3]) ? true : false;
    }

    public String toString() {
	return "Face" + index + "[" + edges[0] + "," + edges[1] + "," +
				      edges[2] + "," + edges[3] + "]" + (isAmbiguous() ? "*" : "");
    }
}

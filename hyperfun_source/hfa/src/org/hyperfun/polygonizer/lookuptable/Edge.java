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
 * $Id: Edge.java,v 1.1 2006/11/08 09:13:15 y7goto Exp $
 */

package org.hyperfun.polygonizer.lookuptable;

/**
 * @author Yuichiro Goto
 * @version $Revision: 1.1 $
 */
public final class Edge {
    private static final int[][] EDGE_VERTICES = {
	{0, 1}, {1, 2}, {3, 2}, {0, 3},
	{4, 5}, {5, 6}, {7, 6}, {4, 7},
	{0, 4}, {1, 5}, {2, 6}, {3, 7}
    };
    private int index;
    private int startVertexIndex;
    private int endVertexIndex;
    private Edge connectedEdge0;
    private Edge connectedEdge1;

    Edge(int index) {
	if (index < 0 || index > 11) {
	    throw new IllegalArgumentException(
		"index must be in the range between 0 and 11");
	}
	this.index = index;
	this.startVertexIndex = EDGE_VERTICES[index][0];
	this.endVertexIndex = EDGE_VERTICES[index][1];
	this.connectedEdge0 = null;
	this.connectedEdge1 = null;
    }

    public int getIndex() {
	return index;
    }

    public int getStartVertexIndex() {
	return startVertexIndex;
    }

    public int getEndVertexIndex() {
	return endVertexIndex;
    }

    void setConnectedEdge(int index, Edge edge) {
	if (index != 0 && index != 1) {
	    throw new IndexOutOfBoundsException();
	}
	if (index == 0) {
	    connectedEdge0 = edge;
	} else {
	    connectedEdge1 = edge;
	}
    }

    public Edge getConnectedEdge(int index) {
	if (index != 0 && index != 1) {
	    throw new IndexOutOfBoundsException();
	}
	return (index == 0) ? connectedEdge0 : connectedEdge1;
    }

    public String toString() {
	return "Edge" + index + "[" + startVertexIndex + "," + endVertexIndex + "]";
    }
}

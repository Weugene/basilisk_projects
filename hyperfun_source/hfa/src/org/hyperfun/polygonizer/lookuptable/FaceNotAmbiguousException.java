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
 * $Id: FaceNotAmbiguousException.java,v 1.1 2006/11/08 09:13:15 y7goto Exp $
 */

package org.hyperfun.polygonizer.lookuptable;

/**
 * @author Yuichiro Goto
 * @version $Revision: 1.1 $
 */
public class FaceNotAmbiguousException extends Exception {
    public FaceNotAmbiguousException() {
	super();
    }

    public FaceNotAmbiguousException(String s) {
	super(s);
    }
}

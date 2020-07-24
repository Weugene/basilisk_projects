/*
 * Copyright (c) 1999 Richard Cartwright  All rights reserved.
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
 * $Id: Loader.java,v 1.2 2006/02/02 01:33:28 y7goto Exp $
 */

package org.hyperfun.EmpiricalHyperFun;

import java.io.*;
import java.security.SecureClassLoader;

/** <P></P>

    @author Richard Cartwright
    @version 2.0 alpha
*/

public class Loader extends SecureClassLoader {

    /*
    final static Loader loader = new Loader();

    static Package ehfPackage = null;

    static {
	ehfPackage = loader.definePackage("EmpiricalHyperFun.Generated",
					  "EHF Generated Classes",
					  "0.1",
					  "HyperFun Consortium",
					  "EHF Generated Classes",
					  "0.1",
					  "HyperFun Consortium", null);
    }
    */

    /** <P>Loader is user as a utility and its constructor is therefore
	trivial.  The superclass constructor is delegated with the
	responsibility of setting up a working class loader based on
	the system class loaded
	<code>Class.getSystemClassLoader()</code>.</P>
     */
    public Loader() {
    }

    public Loader(ClassLoader parent) {
      super(parent);
    }

    /** <P>Create a class from a byte array containing class data.  The
	name of the class must match the name of the class encoded in
	the class data.  The processing of the class data relies on
	the methods <code>ClassLoader.defineClass()</code> and
	<code>ClassLoader.resolveClass()<code>.</P>

	@param name Name of the class to be loaded.
	@param classData Representation of the class as an array of bytes.
	@return A <code>java.land.Class</code> representing the linked
	class described in the passed byte data.
     */

    public Class classFromBytes(String name, byte[] classData)
	throws IndexOutOfBoundsException, SecurityException,
	       ClassFormatError, NullPointerException {

	Class cl = null;

	cl = defineClass(name, classData, 0, classData.length);

        resolveClass(cl);

	return cl;
    }
}

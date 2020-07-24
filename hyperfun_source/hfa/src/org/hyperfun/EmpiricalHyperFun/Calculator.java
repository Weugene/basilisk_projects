/*
 * Copyright (c) 2001 Kazuhiro Mochizuki
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
 * $Id: Calculator.java,v 1.3 2006/02/02 01:33:27 y7goto Exp $
 */

package org.hyperfun.EmpiricalHyperFun;

import java.io.*;
import java.lang.*;
import java.lang.reflect.Method;
import java.util.Hashtable;
import java_cup.runtime.*;
import javax.vecmath.Vector3f;
import org.apache.bcel.generic.*;
import org.apache.bcel.*;
import java.util.*;

import org.hyperfun.polygonizer.Function;

public class Calculator implements Constants, Function {
    Method method;
    Object obj = new Object();

    public Calculator(FileInputStream file) throws Exception {
	//try {
	//  Parser p = new Parser(new Yylex(file));
	//  cal_Init( p );
	//} catch ( Exception e ) {
	//  System.err.println( e );
	//}
	Parser p = new Parser(new Yylex(file));
	cal_Init(p);
    }

    public Calculator(String in) throws Exception {
	//try {
	//  Parser p = new Parser(new Yylex(new StringReader(in)));
	//  cal_Init( p );
	//} catch ( Exception e ) {
	//  System.err.println( e );
	//}
	Parser p = new Parser(new Yylex(new StringReader(in)));
	cal_Init(p);
    }

    private void cal_Init(Parser p) throws Exception {
	Symbol   top;
	Vector   bytes;
	ClassGen cg;

	int avg = 0;

	cg = new ClassGen("GeneratedClass", "org.hyperfun.EmpiricalHyperFun.FrepType", "<generated>",
			  ACC_PUBLIC | ACC_SUPER, null);

	p.setClassGen(cg, "f");
	p.setVarTable(new Hashtable());

	top = p.parse();
	if (((Vector)top.value).size() != 0) {
	    StringBuffer buf = new StringBuffer();

	    for (Enumeration e = ((Vector)top.value).elements(); e.hasMoreElements();) {
		buf.append((String)e.nextElement() + "\n");
	    }
	    throw new Exception(buf.toString());
	}

	FieldGen test_field = new FieldGen(ACC_PUBLIC, Type.DOUBLE, "test_field",
					   cg.getConstantPool());
	cg.addField(test_field.getField());

	cg.addEmptyConstructor(ACC_PUBLIC);

	Loader load = new Loader(this.getClass().getClassLoader());

	byte[] input = cg.getJavaClass().getBytes();
	// For DEBUG
	//cg.getJavaClass().dump(new FileOutputStream("dump.class"));
	Class tmp = load.classFromBytes("GeneratedClass", input);
	obj = tmp.newInstance();

	java.lang.reflect.Method[] M_test = tmp.getDeclaredMethods();

	method = M_test[0];
    }

    public double evaluate(double x, double y, double z) throws RuntimeException {
	double result = 0.0;

	double[] X = {x, y, z};
	double[] S = {};
	double[] A = {};

	Object[] IN = {X, S, A};

	try {
	    result = ((Double)method.invoke(obj, IN)).doubleValue();
	} catch (Exception e) {
	    throw new RuntimeException("Calculation error");
	}
	return result;
    }
}

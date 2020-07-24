/*
 * Copyright (c) 1999 Richard Cartwright
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
 * $Id: HFParameter.java,v 1.2 2006/02/02 01:33:27 y7goto Exp $
 */

package org.hyperfun.EmpiricalHyperFun;

import org.apache.bcel.generic.*;
import org.apache.bcel.classfile.*;
import org.apache.bcel.*;

public class HFParameter implements Constants {

    public static final int NUMBER = 1;
    public static final int ARRAY = 2;
    public static final int FUNCTION = 3;

    public String name;
    public int type;
    private int index = -1;
    Field field;
    boolean local = false;

    public HFParameter() {
    }

    HFParameter(String name, int type, int index) {
	local = true;
	this.name = name;
	this.type = type;
	this.index = index;
	field = null;
    }

    int getIndex() {
	return index;
    }

    void setIndex(int idx) {
	index = idx;
    }

    void addField(ConstantPoolGen cp, ClassGen cg) {
	switch (type) {
	case NUMBER:
	    field = (new FieldGen(ACC_PUBLIC, Type.DOUBLE,
				  name, cp)).getField();
	    cg.addField(field);
	    index = cp.addFieldref(cg.getClassName(), name, "D");
	    break;

	case ARRAY:
	    field = (new FieldGen(ACC_PUBLIC, new ArrayType(Type.DOUBLE, 1),
				  name, cp)).getField();
	    cg.addField(field);
	    index = cp.addFieldref(cg.getClassName(), name, "[D");
	    break;

	case FUNCTION:
	    field = (new FieldGen(ACC_PUBLIC, new ObjectType("org.hyperfun.EmpiricalHyperFun.FrepType"),
			 name, cp)).getField();
	    cg.addField(field);
	    index = cp.addFieldref(cg.getClassName(), name, "Lorg.hyperfun.EmpiricalHyperFun/FrepType;");
	    break;

	default:
	    break;
	}
    }
}

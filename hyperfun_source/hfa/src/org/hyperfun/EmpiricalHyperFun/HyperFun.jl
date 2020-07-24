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
 * $Id: HyperFun.jl,v 1.4 2006/10/07 15:48:14 y7goto Exp $
 */

package org.hyperfun.EmpiricalHyperFun;

import java.io.*;
import java_cup.runtime.*;

%%
%type Symbol
%eofval{
  return new Symbol(Sym.EOF);
%eofval}
%state COMMENT
%cup
%public
%unicode
%notunix
%line
digit = [0-9]
letter = [a-zA-Z]
unary1 = sqrt|exp|log|sin|cos|tan|asin|acos|atan|abs
unary2 = logd
unary3 = sinh|cosh|tanh|sign
binary1 = max|min|atan2
binary2 = mod
%{
	public int getCurrentLine() {
           return yyline;
        }
%}
%%
<YYINITIAL>	--	{ yybegin(COMMENT); }
<COMMENT>	.	{ }
<COMMENT>	[\n\r]	{ yybegin(YYINITIAL); }
<YYINITIAL>	\(	{ return new Symbol(Sym.LPAREN); }
<YYINITIAL>	\)	{ return new Symbol(Sym.RPAREN); }
<YYINITIAL>	,	{ return new Symbol(Sym.COMMA); }
<YYINITIAL>	\{	{ return new Symbol(Sym.LCURL); }
<YYINITIAL>	\}	{ return new Symbol(Sym.RCURL); }
<YYINITIAL>	\;	{ return new Symbol(Sym.SEMI); }
<YYINITIAL>	\=	{ return new Symbol(Sym.EQ); }
<YYINITIAL>	\.	{ return new Symbol(Sym.DOT); }
<YYINITIAL>	\[	{ return new Symbol(Sym.LSQ); }
<YYINITIAL>	\]	{ return new Symbol(Sym.RSQ); }
<YYINITIAL>	"/="	{ return new Symbol(Sym.NEQ); }
<YYINITIAL>	"<"	{ return new Symbol(Sym.LT); }
<YYINITIAL>	">"	{ return new Symbol(Sym.GT); }
<YYINITIAL>	"<="	{ return new Symbol(Sym.LEQ); }
<YYINITIAL>	">="	{ return new Symbol(Sym.GEQ); }
<YYINITIAL>	\+	{ return new Symbol(Sym.PLUS); }
<YYINITIAL>	-	{ return new Symbol(Sym.MINUS); }
<YYINITIAL>	\*	{ return new Symbol(Sym.STAR); }
<YYINITIAL>	"/"	{ return new Symbol(Sym.DIV); }
<YYINITIAL>	"^"	{ return new Symbol(Sym.POW); }
<YYINITIAL>	\|	{ return new Symbol(Sym.UNION); }
<YYINITIAL>	"&"	{ return new Symbol(Sym.INTER); }
<YYINITIAL>	\\	{ return new Symbol(Sym.SUB); }
<YYINITIAL>	\~	{ return new Symbol(Sym.NEG); }
<YYINITIAL>	array	{ return new Symbol(Sym.ARRAY); }
<YYINITIAL>	if	{ return new Symbol(Sym.IF); }
<YYINITIAL>	then	{ return new Symbol(Sym.THEN); }
<YYINITIAL>	else	{ return new Symbol(Sym.ELSE); }
<YYINITIAL>	endif	{ return new Symbol(Sym.ENDIF); }
<YYINITIAL>	while	{ return new Symbol(Sym.WHILE); }
<YYINITIAL>	loop	{ return new Symbol(Sym.LOOP); }
<YYINITIAL>	endloop	{ return new Symbol(Sym.ELOOP); }
<YYINITIAL>	and	{ return new Symbol(Sym.AND); }
<YYINITIAL>	or	{ return new Symbol(Sym.OR); }
<YYINITIAL>	not	{ return new Symbol(Sym.NOT); }
<YYINITIAL>	{unary1}	{ return new Symbol(Sym.UNARY1, yytext()); }
<YYINITIAL>	{unary2}	{ return new Symbol(Sym.UNARY2, yytext()); }
<YYINITIAL>	{unary3}	{ return new Symbol(Sym.UNARY3, yytext()); }
<YYINITIAL>	{binary1} 	{ return new Symbol(Sym.BIN1, yytext()); }
<YYINITIAL>	{binary2}	{ return new Symbol(Sym.BIN2, yytext()); }
<YYINITIAL>	{digit}+	{ return new Symbol(Sym.INT, 
					new Integer(yytext())); }
<YYINITIAL>	{digit}+(\.{digit}*)?([Ee][\+\-]?{digit}+)?
				{ return new Symbol(Sym.NUM,
					new Double(yytext())); }
<YYINITIAL>	{digit}*(\.{digit}+)([Ee][\+\-]?{digit}+)? {
				return new Symbol(Sym.NUM, 
					new Double(yytext())); }
<YYINITIAL>	hfsphere        { return new Symbol(Sym.HFSPHERE); }
<YYINITIAL>	hfellipsoid     { return new Symbol(Sym.HFELLIPSOID); }
<YYINITIAL>	hfcylinderx     { return new Symbol(Sym.HFCYLX); }
<YYINITIAL>	hfcylindery	{ return new Symbol(Sym.HFCYLY); }
<YYINITIAL>	hfcylinderz	{ return new Symbol(Sym.HFCYLZ); }
<YYINITIAL>	hfellcylx	{ return new Symbol(Sym.HFELLCYLX); }
<YYINITIAL>	hfellcyly	{ return new Symbol(Sym.HFELLCYLY); }
<YYINITIAL>	hfellcylz	{ return new Symbol(Sym.HFELLCYLZ); }
<YYINITIAL>	hftorusx	{ return new Symbol(Sym.HFTORUSX); }
<YYINITIAL>	hftorusy	{ return new Symbol(Sym.HFTORUSY); }
<YYINITIAL>	hftorusz	{ return new Symbol(Sym.HFTORUSZ); }
<YYINITIAL>	hfblock		{ return new Symbol(Sym.HFBLOCK); }
<YYINITIAL>	hfblobby	{ return new Symbol(Sym.HFBLOBBY); }
<YYINITIAL>	hfmetaball	{ return new Symbol(Sym.HFMETABALL); }
<YYINITIAL>	hfsoft		{ return new Symbol(Sym.HFSOFT); }
<YYINITIAL>	hfblenduni	{ return new Symbol(Sym.HFBLENDUNI); }
<YYINITIAL>	hfblendint	{ return new Symbol(Sym.HFBLENDINT); }
<YYINITIAL>	hfscale3d	{ return new Symbol(Sym.HFSCALE3D); /* [D */ }
<YYINITIAL>	hfshift3d	{ return new Symbol(Sym.HFSHIFT3D); /* [D */ }
<YYINITIAL>	hfrotate3dx	{ return new Symbol(Sym.HFROTATE3DX); /* [D */ }
<YYINITIAL>	hfrotate3dy	{ return new Symbol(Sym.HFROTATE3DY); /* [D */ }
<YYINITIAL>	hfrotate3dz	{ return new Symbol(Sym.HFROTATE3DZ); /* [D */ }
<YYINITIAL>	hftwistx	{ return new Symbol(Sym.HFTWISTX); /* [D */ }
<YYINITIAL>	hftwisty	{ return new Symbol(Sym.HFTWISTY); /* [D */ }
<YYINITIAL>	hftwistz	{ return new Symbol(Sym.HFTWISTZ); /* [D */ }
<YYINITIAL>	hfstretch3d	{ return new Symbol(Sym.HFSTRETCH3D); /* [D */ }
<YYINITIAL>	hftaperx	{ return new Symbol(Sym.HFTAPERX); /* [D */ }
<YYINITIAL>	hftapery	{ return new Symbol(Sym.HFTAPERY); /* [D */ }
<YYINITIAL>	hftaperz	{ return new Symbol(Sym.HFTAPERZ); /* [D */ }
<YYINITIAL>	hfconvarc	{ return new Symbol(Sym.HFCONVARC); }
<YYINITIAL>	hfconvcurve	{ return new Symbol(Sym.HFCONVCURVE); }
<YYINITIAL>	hfconvline	{ return new Symbol(Sym.HFCONVLINE); }
<YYINITIAL>	hfconvmesh	{ return new Symbol(Sym.HFCONVMESH); }
<YYINITIAL>	hfconvpoint	{ return new Symbol(Sym.HFCONVPOINT); }
<YYINITIAL>	hfconvtriangle	{ return new Symbol(Sym.HFCONVTRIANGLE); }
<YYINITIAL>	hfconvliner	{ return new Symbol(Sym.HFCONVLINER); }
<YYINITIAL>	hfconex		{ return new Symbol(Sym.HFCONEX); }
<YYINITIAL>	hfconey		{ return new Symbol(Sym.HFCONEY); }
<YYINITIAL>	hfconez		{ return new Symbol(Sym.HFCONEZ);}
<YYINITIAL>	hfellconex	{ return new Symbol(Sym.HFELLCONEX); }
<YYINITIAL>	hfellconey	{ return new Symbol(Sym.HFELLCONEY); }
<YYINITIAL>	hfellconez	{ return new Symbol(Sym.HFELLCONEZ); }
<YYINITIAL>	hffmapblob	{ return new Symbol(Sym.HFFMAPBLOB); }
<YYINITIAL>	hfnoiseg	{ return new Symbol(Sym.HFNOISEG); }
<YYINITIAL>	hfspacemapcubic { return new Symbol(Sym.HFSPACEMAPCUBIC); }
<YYINITIAL>	hfspacemapexp	{ return new Symbol(Sym.HFSPACEMAPEXP); }
<YYINITIAL>	hfsuperell	{ return new Symbol(Sym.HFSUPERELL); }
<YYINITIAL>	{letter}({letter}|{digit}|_)*	
				{ return new Symbol(Sym.ID, yytext()); }
<YYINITIAL>	" "|\n|\r|\t	{ }
	

#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*uoldxi=1;*)
(*uoldxg=1;*)
(*dt=0.01;*)
(*CC=1;*)
(*Delta=4/2^9;*)
(*uxi=1;*)
(*Print[(uoldxi+uoldxg+(2.0*(dt*CC/Delta)-1)*uxi)/(1.0+2.0*dt*CC/Delta)]*)


(* ::Input:: *)
(*Print[(2.0*(dt*CC/Delta)-1)*(((uoldxi+uoldxg)/(2.0*dt*CC-Delta))*Delta+ uxi)/(1+2.0*dt*CC/Delta)]*)
(**)


(* ::Input:: *)
(*Print[(uoldxi+uoldxg+(2.0*(dt*CC/Delta)-1)*uxi) - (((uoldxi+uoldxg)/(2.0*dt*CC-Delta))*Delta+ uxi)]*)


(* ::Input:: *)
(*u[x_,y_]:=0.01*Sin[2*Pi*x]*Sin[2*Pi*y];*)
(*v[x_,y_]:=0.01*Cos[2*Pi*x]*Cos[2*Pi*y];*)
(**)
(**)

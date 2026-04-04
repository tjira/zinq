#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"];$Assumptions={rho\[Element]Reals,rho>0,rs\[Element]Reals,rs>0};


(* ::Subsection:: *)
(*Helper Functions*)


FormatExpression[expr_,format_]:=Fold[FullSimplify[#1/.#2]&,expr,format]


(* ::Subsection:: *)
(*Slater Exchange*)


Ex=-Cx*rho^(4/3);


eps=Ex/rho
v=D[Ex,rho];
f=D[v,rho];


FullSimplify[FormatExpression[v,{}]-eps]+"eps"
FormatExpression[f,{}]


(* ::Subsection:: *)
(*Chachiyo Correlation*)


(* ::Subsubsection:: *)
(*rs=cbrt(3/(4 Pi rho))*)


eps=a*Log[1+b/rs+b/rs^2]
v=eps-rs*D[eps,rs]/3;
f=-4*Pi*rs^4*D[v,rs]/9;


FullSimplify[FormatExpression[v,{}]-eps]+"eps"
FormatExpression[Together@f,{}]


(* ::Subsection:: *)
(*VWN5 Correlation*)


(* ::Subsubsection:: *)
(*x=sqrt(x)*)
(*Q=sqrt(4c-b^2)*)
(*X(x)=x^2+bx+c*)


eps=A*(Log[x^2/X[x]]+2*b/Q*ArcTan[Q/(2*x+b)]-b*x0/X[x0]*(Log[(x-x0)^2/X[x]]+2*(b+2*x0)/Q*ArcTan[Q/(2x+b)]))
v=eps-x*D[eps,x]/6;
f=-2*Pi*x^7*D[v,x]/9;


format={
	Q^2->4*c-b^2,
	c+x*(b+x)->X[x],
	b+2*x0->X'[x0],
	X''[x]->2
};


FullSimplify[FormatExpression[v,format]-eps]+"eps"
FormatExpression[f,format]




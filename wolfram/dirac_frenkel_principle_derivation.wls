#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


(* ::Section:: *)
(*Assumptions*)


$Assumptions={k\[Element]Reals,m\[Element]Reals,x\[Element]Reals,q[t]\[Element]Reals,p[t]\[Element]Reals,\[Gamma]\[Element]Reals,\[Gamma]>0};


(* ::Section:: *)
(*Potential*)


V=k*x^2/2;


(* ::Section:: *)
(*Wavefunction and Parameters*)


psi=Exp[-\[Gamma]/2*(x-q[t])^2+I*p[t]*(x-q[t])];


psi=psi/Sqrt@Integrate[Conjugate@psi*psi,{x,-Infinity,Infinity}];


phi={q[t],p[t]};


(* ::Section:: *)
(*Dirac-Frenkel Principle System*)


S=Table[Integrate[FullSimplify[Conjugate@D[psi,phi[[i]]]*D[psi,phi[[j]]]],{x,-Infinity,Infinity}],{i,1,Length@phi},{j,1,Length@phi}];


dEkin=Transpose@{Grad[Integrate[FullSimplify[Conjugate@psi*(-D[psi,{x,2}]/(2*m))],{x,-Infinity,Infinity}],phi]};


dEpot=Transpose@{Grad[Integrate[FullSimplify[Conjugate@psi*V*psi],{x,-Infinity,Infinity}],phi]};


dE=dEkin+dEpot;dE//MatrixForm


F=-2*FullSimplify@Im@S;F//MatrixForm


dPhi=Dot[Inverse@F,dE];dPhi//MatrixForm

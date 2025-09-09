#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]


(* ::Section:: *)
(*Assumptions*)


$Assumptions={k\[Element]Reals,m\[Element]Reals,x\[Element]Reals,q[t]\[Element]Reals,p[t]\[Element]Reals,gamma\[Element]Reals,gamma>0};


(* ::Section:: *)
(*Potential*)


V=k*x^2/2;


(* ::Section:: *)
(*Wavefunction and Parameters*)


psi=Exp[-gamma/2*(x-q[t])^2+I*p[t]*x];


psi=psi/Sqrt@Integrate[Conjugate@psi*psi,{x,-Infinity,Infinity}]


phi={q[t],p[t]};


(* ::Section:: *)
(*Dirac-Frenkel Principle System*)


S=Table[Integrate[FullSimplify[Conjugate@D[psi,phi[[i]]]*D[psi,phi[[j]]]],{x,-Infinity,Infinity}],{i,1,Length@phi},{j,1,Length@phi}];


dE=Transpose@{Grad[Integrate[FullSimplify[Conjugate@psi*(V*psi-D[psi,{x,2}]/(2*m))],{x,-Infinity,Infinity}],phi]};dE//MatrixForm


F=-2*FullSimplify@Im@S;F//MatrixForm


dPhi=Inverse@F . dE;dPhi//MatrixForm

#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"];


(* ::Subsection:: *)
(*Assumptions*)


$Assumptions={Subscript[x,1]\[Element]Reals,Subscript[x,2]\[Element]Reals,Subscript[p,1]\[Element]Reals,Subscript[p,2]\[Element]Reals,Re@Subscript[g,1]>0,Re@Subscript[g,2]>0};


(* ::Subsection:: *)
(*Initial Conditions*)


IC={m->5,Subscript[x,1]->1,Subscript[x,2]->2,Subscript[p,1]->2,Subscript[p,2]->1,Subscript[g,1]->1+2I,Subscript[g,2]->3+4I};


(* ::Subsection:: *)
(*Wavefunction Definition*)


G1=Exp[-Subscript[g,1]*(q-Subscript[x,1])^2/2+I*Subscript[p,1]*(q-Subscript[x,1])];
G2=Exp[-Subscript[g,2]*(q-Subscript[x,2])^2/2+I*Subscript[p,2]*(q-Subscript[x,2])];


G1=G1/Sqrt@Integrate[FullSimplify[Conjugate@G1*G1,q\[Element]Reals],{q,-Infinity,Infinity}];
G2=G2/Sqrt@Integrate[FullSimplify[Conjugate@G2*G2,q\[Element]Reals],{q,-Infinity,Infinity}];


(* ::Subsection:: *)
(*Overlap Integral*)


S=Integrate[FullSimplify[Conjugate@G1*G2,{q\[Element]Reals}],{q,-Infinity,Infinity}]


(* ::Subsection:: *)
(*Kinetic Integral*)


T=Integrate[FullSimplify[-Conjugate@G1*D[G2,{q,2}]/(2*m),{m\[Element]Reals,q\[Element]Reals}],{q,-Infinity,Infinity}]


(* ::Subsection:: *)
(*Time Derivative Integral*)


taux=Integrate[FullSimplify[Conjugate@G1*D[G2,Subscript[x,2]],{q\[Element]Reals}],{q,-Infinity,Infinity}]
taup=Integrate[FullSimplify[Conjugate@G1*D[G2,Subscript[p,2]],{q\[Element]Reals}],{q,-Infinity,Infinity}]


(* ::Subsection:: *)
(*Integral Values*)


N[S/.IC]
N[T/.IC]


dx=N[taux/.IC]
dp=N[taup/.IC]




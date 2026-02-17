#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"];Off[General::munfl];


(* ::Subsection:: *)
(*Variables*)


kB=1.380649*^-23;
Eh=4.359744722206*^-18;
Na=6.02214076*^23;
T=1000000000;
m=1;
k=1;


kBT=kB*T/Eh;
omega=Sqrt[k/m];
sigma=kBT/m/omega^2;


(* ::Subsection:: *)
(*Data*)


data=RandomVariate[NormalDistribution[0,Sqrt@sigma],10000000];


(* ::Subsection:: *)
(*Exact Entropy*)


Sexact=kB*Na*(1/(kBT*(Exp[1/kBT]-1))-Log[1-Exp[-1/kBT]])


(* ::Subsection:: *)
(*Schlitter Entropy*)


Sschlitter=0.5*kB*Na*Log[1+kBT*Exp[2]*m*Covariance@data]

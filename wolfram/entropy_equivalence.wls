#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"];Off[General::munfl];


kB=1.380649*^-23;
Eh=4.359744722206*^-18;
T=1000000000;
m=1;
k=1;


kBT=kB*T/Eh;
omega=Sqrt[k/m];
sigma=kBT/m/omega^2;


Sschlitter=0.5*kB/Eh*Log[1+kBT*Exp[2]*m*Covariance@RandomVariate[NormalDistribution[0,Sqrt@sigma],10000000]]


Sexact=kB/Eh*(1/(kBT*(Exp[1/kBT]-1))-Log[1-Exp[-1/kBT]])

(* ::Package:: *)

(* ::Section:: *)
(*Package for Linear thermo-viscoelasticity*)


(* ::Text:: *)
(*This first package carries with tensor calculus operations in  both harmonic and Hill's decomposition*)


(* ::Input:: *)
(*(* Name: LinThViscoElast *)*)
(*(* Author: Camilo Andr\[EAcute]s SUAREZ AFANADOR *)*)
(*(* Date : May 4,2021 *)*)
(* *)


BeginPackage["LinThViscoElast`"]

RelaxTimeFun::usage = 
"RelaxTimeFun[\"mo\",\"\[Tau]\",\"m\",\"var\"] Gives a discrete viscoelastic function in time-domain given: \
 Instantaneous modulus \"\!\(\*SubscriptBox[\(m\), \(o\)]\)\",  the vectors \"\!\(\*SubscriptBox[\(\[Tau]\), \
 \(\(i\)\(\\\ \)\)]\)\" and \"\!\(\*SubscriptBox[\(m\), \(i\)]\)\" of the discrete distribution and \"var\" the symbol of the time variable"
RelaxFreqFun::usage = 
"RelaxFreqFun[\"mo\",\"\[Tau]\",\"m\",\"var\"] Gives a discrete viscoelastic function in the frequency domain given: \ 
Instantaneous \"\!\(\*SubscriptBox[\(m\), \(o\)]\)\", the vectors \"\!\(\*SubscriptBox[\(\[Tau]\), \(\(i\)\(\\\ \)\)]\)\" and \"\!\(\*SubscriptBox[\(m\), \(i\)]\)\"
of the discrete distribution and \"var\" the symbol of the time variable"
PlotRelax::usage = 
"PlotRelax[\"log10start\",\"log10end\",\"RelaxFun\",\"var\",\"name\" ] \
plot a viscoelastic function in time domain given: \
the base10 logarithm of the initial and end time (\"log10start\",\"log10end\"), \
the function \"RelaxFun\", the argument's function \"var\", and an optional Title \"name\" "
PlotDMA::usage = 
"PlotDMA[\"log10start\",\"log10end\",\"FourFun\",\"var\",\"name\",\"scaleFun\"] \
plot a viscoelastic function in frequency domain given: \
the base10 logarithm of the initial and ending frequency (\"log10start\",\"log10end\"), \
 the function \"FourFun\", the argument's function \"var\", \
 an two optionals, Title \"name\" and \"scalefuns\" specifying the \
axis scale (by default = {\"Log10\",\"Linear\"})"
Epsthvec::usage= 
"\[Epsilon]thvec[\"Tlist\", \"\[Alpha]tg\", \"\[Alpha]ti\", \"\[Tau]\[Alpha]ti\", \"\[Alpha]lg\", \"\[Alpha]li\", \"\[Tau]\[Alpha]li\", \"\[CapitalDelta]timeList\"]
Gives a list of the evolution of the thermal-strain {\[Epsilon]11,\[Epsilon]22,\[Epsilon]33,\[Epsilon]12,\[Epsilon]13,\[Epsilon]23} as a function of: \
Tlist, the list of temperature evolution \
\[Alpha]tg, the glassy transverse expansion coefficient\
\[Alpha]ti, the list of the discrete distriution of transverse expansion coefficient\
\[Tau]\[Alpha]ti, the list of the discrete distriution of transverse expansion times\
\[Alpha]lg, the glassy axial expansion coefficient\
\[Alpha]li, the list of the discrete distriution of axial expansion coefficient\
\[Tau]\[Alpha]li, the list of the discrete distriution of axial expansion times\
\[CapitalDelta]timeList, the list of time steps starting with zero, {0,\[CapitalDelta]t2,...,\[CapitalDelta]tN} "
Epsthvec2::usage= 
"The same as Epstvec, adding a Taylor expansion of 2nd order whe \[CapitalDelta]\[Xi]/\[Tau] -> 0"
EpsthvecEulerImp::usage= 
"\[Epsilon]thvec[\"Tlist\", \"\[Alpha]tg\", \"\[Alpha]ti\", \"\[Tau]\[Alpha]ti\", \"\[Alpha]lg\", \"\[Alpha]li\", \"\[Tau]\[Alpha]li\", \"\[CapitalDelta]timeList\"]
Gives a list of the evolution of the thermal-strain {\[Epsilon]11,\[Epsilon]22,\[Epsilon]33,\[Epsilon]12,\[Epsilon]13,\[Epsilon]23} using Implicit Euler's scheme as a function of: \
Tlist, the list of temperature evolution \
\[Alpha]tg, the glassy transverse expansion coefficient\
\[Alpha]ti, the list of the discrete distriution of transverse expansion coefficient\
\[Tau]\[Alpha]ti, the list of the discrete distriution of transverse expansion times\
\[Alpha]lg, the glassy axial expansion coefficient\
\[Alpha]li, the list of the discrete distriution of axial expansion coefficient\
\[Tau]\[Alpha]li, the list of the discrete distriution of axial expansion times\
\[CapitalDelta]timeList, the list of time steps starting with zero, {0,\[CapitalDelta]t2,...,\[CapitalDelta]tN} "
TIsoTensorStress::usage = 
"TensorStress[\"Hillg\", \"Hilli\", \"Hill\[Tau]i\", \"\[Epsilon]mec\", \"dir\", \"\[CapitalDelta]timeList\"] \
Gives a list of the evolution of the stress tensor as a function of: \
Hillg, the 6 length list of instantaneous coefs of the materials in Hill's basis {2\!\(\*SubscriptBox[\(K\), \(g\)]\),\!\(\*SubscriptBox[\(l\), \(g\)]\),\!\(\*SubscriptBox[\(l\), \(g\)]\),\!\(\*SubscriptBox[\(n\), \(g\)]\),2\!\(\*SubscriptBox[\(m\), \(g\)]\),2\!\(\*SubscriptBox[\(g\), \(g\)]\)}. \
Hilli, the 6 length list of lists of the coefs distributions of the material in Hill's basis. \
Hill\[Tau]i, the 6 length list of lists of the time distributions of the material in Hill's basis.  \
\[Epsilon]mec, the list of time evolution of mechanical strains as {\!\(\*SubscriptBox[\(\[Epsilon]\), \(\(11\)\(,\)\)]\) \!\(\*SubscriptBox[\(\[Epsilon]\), \(\(22\)\(,\)\)]\) \!\(\*SubscriptBox[\(\[Epsilon]\), \(33\)]\), \!\(\*SubscriptBox[\(\[Epsilon]\), \(12\)]\), \!\(\*SubscriptBox[\(\[Epsilon]\), \(13\)]\), \!\(\*SubscriptBox[\(\[Epsilon]\), \(23\)]\)}for each step. \
dir, the principal direction of the material symmetry (1 or 2 or 3). \
\[CapitalDelta]timeList, the list of time steps starting with zero, {0,\[CapitalDelta]t2,...,\[CapitalDelta]tN}. "
TIsoTensorStressEulerImp::usage = 
"TensorStress[\"Hillg\", \"Hilli\", \"Hill\[Tau]i\", \"\[Epsilon]mec\", \"dir\", \"\[CapitalDelta]timeList\"] \
Gives a list of the evolution of the stress tensor  using Implicit Euler's scheme as a function of: \
Hillg, the 6 length list of instantaneous coefs of the materials in Hill's basis {2\!\(\*SubscriptBox[\(K\), \(g\)]\),\!\(\*SubscriptBox[\(l\), \(g\)]\),\!\(\*SubscriptBox[\(l\), \(g\)]\),\!\(\*SubscriptBox[\(n\), \(g\)]\),2\!\(\*SubscriptBox[\(m\), \(g\)]\),2\!\(\*SubscriptBox[\(g\), \(g\)]\)}. \
Hilli, the 6 length list of lists of the coefs distributions of the material in Hill's basis. \
Hill\[Tau]i, the 6 length list of lists of the time distributions of the material in Hill's basis.  \
\[Epsilon]mec, the list of time evolution of mechanical strains as {\!\(\*SubscriptBox[\(\[Epsilon]\), \(\(11\)\(,\)\)]\) \!\(\*SubscriptBox[\(\[Epsilon]\), \(\(22\)\(,\)\)]\) \!\(\*SubscriptBox[\(\[Epsilon]\), \(33\)]\), \!\(\*SubscriptBox[\(\[Epsilon]\), \(12\)]\), \!\(\*SubscriptBox[\(\[Epsilon]\), \(13\)]\), \!\(\*SubscriptBox[\(\[Epsilon]\), \(23\)]\)}for each step. \
dir, the principal direction of the material symmetry (1 or 2 or 3). \
\[CapitalDelta]timeList, the list of time steps starting with zero, {0,\[CapitalDelta]t2,...,\[CapitalDelta]tN}. "
IsotropicStress::usage = 
"IsotropicStress[\"ag\", \"ai\", \"\[Tau]i\", \"\[Epsilon]mec\", \"\[CapitalDelta]timeList\"] \
Gives a list of the evolution of the stress tensor as a function of: \
ag, The Isotropic instantaneous coefficients {3\!\(\*SubscriptBox[\(\[Kappa]\), \(g\)]\), 2\!\(\*SubscriptBox[\(\[Mu]\), \(g\)]\)}. \
ai, The lis of isotropic coefficients distribution {{3\!\(\*SubscriptBox[\(\[Kappa]\), \(1\)]\),...,3\!\(\*SubscriptBox[\(\[Kappa]\), \(N\[Kappa]\)]\)},{2\!\(\*SubscriptBox[\(\[Mu]\), \(1\)]\),...,2\!\(\*SubscriptBox[\(\[Mu]\), \(N\[Mu]\)]\)}}. \
\[Tau]i, The lis of times-distribution of coefficients  {{\!\(\*SubscriptBox[\(\[Tau]\[Kappa]\), \(1\)]\),...,\!\(\*SubscriptBox[\(\[Tau]\[Kappa]\), \(N\[Kappa]\)]\)},{\!\(\*SubscriptBox[\(\[Tau]\[Mu]\), \(1\)]\),...,\!\(\*SubscriptBox[\(\[Tau]\[Mu]\), \(N\[Mu]\)]\)}}. \
\[Epsilon]mec, the list of time evolution of mechanical strains as {\!\(\*SubscriptBox[\(\[Epsilon]\), \(\(11\)\(,\)\)]\) \!\(\*SubscriptBox[\(\[Epsilon]\), \(\(22\)\(,\)\)]\) \!\(\*SubscriptBox[\(\[Epsilon]\), \(33\)]\), \!\(\*SubscriptBox[\(\[Epsilon]\), \(12\)]\), \!\(\*SubscriptBox[\(\[Epsilon]\), \(13\)]\), \!\(\*SubscriptBox[\(\[Epsilon]\), \(23\)]\)} for each step. \
\[CapitalDelta]timeList, the list of time steps starting with zero, {0,\[CapitalDelta]t2,...,\[CapitalDelta]tN}. "
IsotropicStress2::usage = 
"The same as IsotropicStress, adding a Taylor expansion of 2nd order whe \[CapitalDelta]\[Xi]/\[Tau] -> 0"
IsotropicStressEulerImp::usage = 
"IsotropicStress[\"ag\", \"ai\", \"\[Tau]i\", \"\[Epsilon]mec\", \"\[CapitalDelta]timeList\"] \
Gives a list of the evolution of the stress tensor using Implicit Euler's scheme as a function of: \
ag, The Isotropic instantaneous coefficients {3\!\(\*SubscriptBox[\(\[Kappa]\), \(g\)]\), 2\!\(\*SubscriptBox[\(\[Mu]\), \(g\)]\)}. \
ai, The lis of isotropic coefficients distribution {{3\!\(\*SubscriptBox[\(\[Kappa]\), \(1\)]\),...,3\!\(\*SubscriptBox[\(\[Kappa]\), \(N\[Kappa]\)]\)},{2\!\(\*SubscriptBox[\(\[Mu]\), \(1\)]\),...,2\!\(\*SubscriptBox[\(\[Mu]\), \(N\[Mu]\)]\)}}. \
\[Tau]i, The lis of times-distribution of coefficients  {{\!\(\*SubscriptBox[\(\[Tau]\[Kappa]\), \(1\)]\),...,\!\(\*SubscriptBox[\(\[Tau]\[Kappa]\), \(N\[Kappa]\)]\)},{\!\(\*SubscriptBox[\(\[Tau]\[Mu]\), \(1\)]\),...,\!\(\*SubscriptBox[\(\[Tau]\[Mu]\), \(N\[Mu]\)]\)}}. \
\[Epsilon]mec, the list of time evolution of mechanical strains as {\!\(\*SubscriptBox[\(\[Epsilon]\), \(\(11\)\(,\)\)]\) \!\(\*SubscriptBox[\(\[Epsilon]\), \(\(22\)\(,\)\)]\) \!\(\*SubscriptBox[\(\[Epsilon]\), \(33\)]\), \!\(\*SubscriptBox[\(\[Epsilon]\), \(12\)]\), \!\(\*SubscriptBox[\(\[Epsilon]\), \(13\)]\), \!\(\*SubscriptBox[\(\[Epsilon]\), \(23\)]\)} for each step. \
\[CapitalDelta]timeList, the list of time steps starting with zero, {0,\[CapitalDelta]t2,...,\[CapitalDelta]tN}. "
DeltaXiListSimp2n::usage =
"DeltaXiListSimp2n[\"afun\", \"var\", \"timelist\", \"Templist\"]
Gives a list of \[CapitalDelta]\[Xi] and \[Xi] values using Simpsons 2n integration rule for a given: \
\"afun\", the shift function. \
\"var\", the argument of the shift function. \
\"timelist\", the list of time values. \
\"Templist\", the list of temperature values. "
DeltaXiIntInt::usage = 
"DeltaXiIntInt[\"afun\", \"var\", \"timelist\", \"Templist\"] 
Gives a list of \[CapitalDelta]\[Xi] values for a given: \
\"afun\", the shift function. \
\"var\", the argument of the shift function. \
\"timelist\", the list of time values. \
\"Templist\", the list of temperature values. "
DeltaXiBooleList::usage = 
"DeltaXiBooleList[\"afun\", \"var\", \"timelist\", \"Templist\"]
Gives a list of \[CapitalDelta]\[Xi] and \[Xi] values using Boole's integration for a given: \
\"afun\", the shift function. \
\"var\", the argument of the shift function. \
\"timelist\", the list of time values. \
\"Templist\", the list of temperature values. "


Begin["`Private`"]

(* 1.Functions to generate viscoelastic functions, given: Instantaneous modulus Subscript[m, o], and the vectors Subscript[m, i] and Subscript[\[Tau], i ]
of the discrete distribution *)
(* Time domain viscoelastic function *)
RelaxTimeFun[mo_,\[Tau]_,m_,var_Symbol]:= Block[{Funt},Funt = mo - Sum[m[[i]](1-E^(-var/\[Tau][[i]])),{i,1,Length[m]}];Funt] 
(* Laplace-Carson domain viscoelastic function *)
RelaxFreqFun[mo_,\[Tau]_,m_,var_Symbol]:= Block[{Funp},Funp = mo - Sum[m[[i]]/\[Tau][[i]] (1/(var+\[Tau][[i]]^-1)),{i,1,Length[m]}]; Funp]

(* 2.Plot moduli *) 
PlotRelax[log10start_,log10end_,RelaxFun_,var_Symbol,name_:"module" ]:=Block[{timelist,Plot,PlotFun},
timelist = Rationalize[Table[10^i,{i,log10start,log10end,0.1}],2*^-16]; 
PlotFun=ListLinePlot[{Thread[{timelist,N[RelaxFun/.{var-> timelist},16]}]},ScalingFunctions->{"Log10","Linear"},PlotRange->All,PlotLabel->name <>" vs t(s)",GridLines->Automatic,Frame->True, Axes->False];
PlotFun]
PlotDMA[log10start_,log10end_,FourFun_,var_Symbol,name_:"module" ,scaleFun_:{"Log10","Linear"}]:=Block[{freqlist,Plot,PlotFun},
freqlist = Rationalize[Table[10^i,{i,log10start,log10end,0.1}],2*^-16]; PlotFun=ListLinePlot[{Thread[{freqlist,N[Re[FourFun /.{var-> I 2\[Pi] freqlist}],16]}],Thread[{freqlist,N[Im[FourFun /.{var-> I 2\[Pi] freqlist}],16]}]},
ScalingFunctions->scaleFun,PlotRange->All,PlotLabel->name <>" vs F(hz)",PlotLegends->{"Storage","Loss"},GridLines->Automatic,Frame->True, Axes->False];
PlotFun]

(* 3.Linear Thermo-Viscoelasticity functions *)

(* 3.1Taylor's integration scheme for solution of ODEs  *)

(* Scheme for a Stress-like Internal variable *)
Taylor[preq_,\[Tau]i_,Mi_,preRhs_,currRhs_,\[CapitalDelta]\[Xi]_]:= Block[{qival},
qival=E^(-\[CapitalDelta]\[Xi]/\[Tau]i)*preq+ Mi (preRhs*(1-E^(-\[CapitalDelta]\[Xi]/\[Tau]i))+ (currRhs-preRhs) *(1-\[Tau]i/\[CapitalDelta]\[Xi] (1-E^(-\[CapitalDelta]\[Xi]/\[Tau]i))));qival]

TaylorT[preq_,\[Tau]i_,Mi_,preRhs_,currRhs_,\[CapitalDelta]\[Xi]_]:= Block[{qival,X,A,B,C},
X = \[CapitalDelta]\[Xi]/\[Tau]i;
Which[
Abs[Log[$MinMachineNumber]*0.95]>X>1.0*^-8,A= E^(-\[CapitalDelta]\[Xi]/\[Tau]i);B=(1-E^(-\[CapitalDelta]\[Xi]/\[Tau]i));C=(1-\[Tau]i/\[CapitalDelta]\[Xi] (1-E^(-\[CapitalDelta]\[Xi]/\[Tau]i))),
Abs[Log[$MinMachineNumber]*0.95]<X,A= 0.0;B=(1.0);C=(1.0),
X<1.0*^-8,A=1-X+X^2/2;B=X-X^2/2;C=X/2-X^2/6];
qival=A*preq+ Mi(preRhs*B + (currRhs-preRhs)*C);qival]

EulerImp2[preq_,\[Tau]i_,Mi_,preRhs_,currRhs_,\[CapitalDelta]\[Xi]_]:= Block[{qival},
qival=((\[Tau]i)*preq+ Mi(currRhs*\[CapitalDelta]\[Xi]))/(\[Tau]i+\[CapitalDelta]\[Xi]);qival]

(* Scheme for a Strain-like Internal variable *)
Taylor2[preq_,\[Tau]i_,preRhs_,currRhs_,\[CapitalDelta]\[Xi]_]:= Block[{qival},
qival=E^(-\[CapitalDelta]\[Xi]/\[Tau]i)*preq+ (preRhs*(1-E^(-\[CapitalDelta]\[Xi]/\[Tau]i))+(currRhs-preRhs)*(1-\[Tau]i/\[CapitalDelta]\[Xi] (1-E^(-\[CapitalDelta]\[Xi]/\[Tau]i))));qival]

TaylorT2[preq_,\[Tau]i_,preRhs_,currRhs_,\[CapitalDelta]\[Xi]_]:= Block[{qival,X,A,B,C},
X = \[CapitalDelta]\[Xi]/\[Tau]i;
Which[
Abs[Log[$MinMachineNumber]*0.95]>X>1.0*^-8,A= E^(-\[CapitalDelta]\[Xi]/\[Tau]i);B=(1-E^(-\[CapitalDelta]\[Xi]/\[Tau]i));C=(1-\[Tau]i/\[CapitalDelta]\[Xi] (1-E^(-\[CapitalDelta]\[Xi]/\[Tau]i))),
Abs[Log[$MinMachineNumber]*0.95]<X,A= 0.0;B=(1.0);C=(1.0),
X<1.0*^-8,A=1-X+X^2/2;B=X-X^2/2;C=X/2-X^2/6];
qival=A*preq+ (preRhs*B + (currRhs-preRhs)*C);qival]

EulerImp[preq_,\[Tau]i_,preRhs_,currRhs_,\[CapitalDelta]\[Xi]_]:= Block[{qival},
qival=((\[Tau]i)*preq+ (currRhs)*\[CapitalDelta]\[Xi])/(\[CapitalDelta]\[Xi]+\[Tau]i) ;qival]

(* 3.2Thermal Deformations functions *)
(* thermal deformations for a principal direction *)
\[Epsilon]th[Tlist_,\[Alpha]g_, \[Alpha]i_,\[Tau]i_, \[CapitalDelta]timeList_]:= Block[{To,q\[Alpha],q\[Alpha]vals,\[Epsilon]th},
To = Tlist[[1]];
q\[Alpha] = Table[0,{i,Length[Tlist]},{j,Length[\[Alpha]i]}];
Do[q\[Alpha][[i,j]]=Taylor[q\[Alpha][[i-1,j]], \[Tau]i[[j]], \[Alpha]i[[j]],(Tlist[[i-1]]-To),(Tlist[[i]]-To),\[CapitalDelta]timeList[[i]]],{i,2,Length[Tlist]},{j,Length[\[Alpha]i]}];
q\[Alpha]vals=Total[q\[Alpha],{2}];
\[Epsilon]th = Table[0,{i,Length[Tlist]}];
Do[\[Epsilon]th[[i]]=\[Alpha]g (Tlist[[i]]-To)-q\[Alpha]vals[[i]],{i,Length[Tlist]}];
\[Epsilon]th]
(* thermal deformations for a transverse-isotropic material *)
Epsthvec[Tlist_,\[Alpha]tg_, \[Alpha]ti_,\[Tau]\[Alpha]ti_,\[Alpha]lg_, \[Alpha]li_,\[Tau]\[Alpha]li_, \[CapitalDelta]timeList_]:=Block[{vec, \[Epsilon]thplan,\[Epsilon]thaxial},
\[Epsilon]thplan=\[Epsilon]th[Tlist,\[Alpha]tg, \[Alpha]ti,\[Tau]\[Alpha]ti, \[CapitalDelta]timeList];\[Epsilon]thaxial=\[Epsilon]th[Tlist,\[Alpha]lg, \[Alpha]li,\[Tau]\[Alpha]li, \[CapitalDelta]timeList];
vec =Table[0,{i,Length[\[Epsilon]thplan]}];
Do[vec[[i]]={\[Epsilon]thplan[[i]],\[Epsilon]thplan[[i]],\[Epsilon]thaxial[[i]],0,0,0},{i,Length[\[Epsilon]thplan]}];
vec]

\[Epsilon]th2[Tlist_,\[Alpha]g_, \[Alpha]i_,\[Tau]i_, \[CapitalDelta]timeList_]:= Block[{To,q\[Alpha],q\[Alpha]vals,\[Epsilon]th},
To = Tlist[[1]];
q\[Alpha] = Table[0,{i,Length[Tlist]},{j,Length[\[Alpha]i]}];
Do[q\[Alpha][[i,j]]=TaylorT[q\[Alpha][[i-1,j]], \[Tau]i[[j]], \[Alpha]i[[j]],(Tlist[[i-1]]-To),(Tlist[[i]]-To),\[CapitalDelta]timeList[[i]]],{i,2,Length[Tlist]},{j,Length[\[Alpha]i]}];
q\[Alpha]vals=Total[q\[Alpha],{2}];
\[Epsilon]th = Table[0,{i,Length[Tlist]}];
Do[\[Epsilon]th[[i]]=\[Alpha]g (Tlist[[i]]-To)-q\[Alpha]vals[[i]],{i,Length[Tlist]}];
\[Epsilon]th]
(* thermal deformations for a transverse-isotropic material *)
Epsthvec2[Tlist_,\[Alpha]tg_, \[Alpha]ti_,\[Tau]\[Alpha]ti_,\[Alpha]lg_, \[Alpha]li_,\[Tau]\[Alpha]li_, \[CapitalDelta]timeList_]:=Block[{vec, \[Epsilon]thplan,\[Epsilon]thaxial},
\[Epsilon]thplan=\[Epsilon]th2[Tlist,\[Alpha]tg, \[Alpha]ti,\[Tau]\[Alpha]ti, \[CapitalDelta]timeList];\[Epsilon]thaxial=\[Epsilon]th2[Tlist,\[Alpha]lg, \[Alpha]li,\[Tau]\[Alpha]li, \[CapitalDelta]timeList];
vec =Table[0,{i,Length[\[Epsilon]thplan]}];
Do[vec[[i]]={\[Epsilon]thplan[[i]],\[Epsilon]thplan[[i]],\[Epsilon]thaxial[[i]],0,0,0},{i,Length[\[Epsilon]thplan]}];
vec]

\[Epsilon]thEulerImp[Tlist_,\[Alpha]g_, \[Alpha]i_,\[Tau]i_, \[CapitalDelta]timeList_]:= Block[{To,q\[Alpha],q\[Alpha]vals,\[Epsilon]th},
To = Tlist[[1]];
q\[Alpha] = Table[0,{i,Length[Tlist]},{j,Length[\[Alpha]i]}];
Do[q\[Alpha][[i,j]]=EulerImp2[q\[Alpha][[i-1,j]], \[Tau]i[[j]], \[Alpha]i[[j]],(Tlist[[i-1]]-To),(Tlist[[i]]-To),\[CapitalDelta]timeList[[i]]],{i,2,Length[Tlist]},{j,Length[\[Alpha]i]}];
q\[Alpha]vals=Total[q\[Alpha],{2}];
\[Epsilon]th = Table[0,{i,Length[Tlist]}];
Do[\[Epsilon]th[[i]]=\[Alpha]g (Tlist[[i]]-To)-q\[Alpha]vals[[i]],{i,Length[Tlist]}];
\[Epsilon]th]
(* thermal deformations for a transverse-isotropic material *)
EpsthvecEulerImp[Tlist_,\[Alpha]tg_, \[Alpha]ti_,\[Tau]\[Alpha]ti_,\[Alpha]lg_, \[Alpha]li_,\[Tau]\[Alpha]li_, \[CapitalDelta]timeList_]:=Block[{vec, \[Epsilon]thplan,\[Epsilon]thaxial},
\[Epsilon]thplan=\[Epsilon]thEulerImp[Tlist,\[Alpha]tg, \[Alpha]ti,\[Tau]\[Alpha]ti, \[CapitalDelta]timeList];\[Epsilon]thaxial=\[Epsilon]thEulerImp[Tlist,\[Alpha]lg, \[Alpha]li,\[Tau]\[Alpha]li, \[CapitalDelta]timeList];
vec =Table[0,{i,Length[\[Epsilon]thplan]}];
Do[vec[[i]]={\[Epsilon]thplan[[i]],\[Epsilon]thplan[[i]],\[Epsilon]thaxial[[i]],0,0,0},{i,Length[\[Epsilon]thplan]}];
vec]




(* thermal deformations for a constant expansion coef consideration (elastic case) *)
\[Epsilon]thvecGlass[Tlist_,\[Alpha]tg_,\[Alpha]lg_,  \[CapitalDelta]timeList_]:=Block[{vec, \[Epsilon]thplan,\[Epsilon]thaxial},\[Epsilon]thplan=\[Epsilon]th[Tlist,\[Alpha]tg, {0},{1}, \[CapitalDelta]timeList];\[Epsilon]thaxial=\[Epsilon]th[Tlist,\[Alpha]lg, {0},{1}, \[CapitalDelta]timeList];
vec =Table[0,{i,Length[\[Epsilon]thplan]}];
Do[vec[[i]]={\[Epsilon]thplan[[i]],\[Epsilon]thplan[[i]],\[Epsilon]thaxial[[i]],0,0,0},{i,Length[\[Epsilon]thplan]}];
vec]

(* 3.3Stress Computations functions *)
SetDirectory[NotebookDirectory[]];
Needs["TensorCalc`"]
(* stress function for isotropic viscoelastic materials *)
(* Entries format:
ag = {3Subscript[\[Kappa], g], 2Subscript[\[Mu], g]}
ai = {{3Subscript[\[Kappa], 1],...,3Subscript[\[Kappa], N\[Kappa]]},{2Subscript[\[Mu], 1],...,2Subscript[\[Mu], N\[Mu]]}}
\[Tau]i = {{Subscript[\[Tau]\[Kappa], 1],...,Subscript[\[Tau]\[Kappa], N\[Kappa]]},{Subscript[\[Tau]\[Mu], 1],...,Subscript[\[Tau]\[Mu], N\[Mu]]}}
\[CapitalDelta]timeList = {0,\[CapitalDelta]t2,...,\[CapitalDelta]tN}/ \[CapitalDelta]ti = Subscript[t, i]-Subscript[t, i-1], \[CapitalDelta]t1 = t1 - t0; t0=t1=0
\[Epsilon]mec = {\[Epsilon]11,\[Epsilon]22,\[Epsilon]33,\[Epsilon]12,\[Epsilon]13,\[Epsilon]23}
*)
IsotropicStress[ag_,ai_, \[Tau]i_,\[Epsilon]mec_, \[CapitalDelta]timeList_]:= Block[
{\[Sigma]Inst,\[Epsilon]Ten,\[Epsilon]TenSph,\[Epsilon]TenDev,\[Epsilon]vSph,\[Epsilon]vDev,\[Sigma]vals,Jijkl,Kijkl},
{Jijkl,Kijkl} = IsoDecomp[];
\[Epsilon]Ten= Table[0,{i,Length[\[Epsilon]mec]}];
Do[\[Epsilon]Ten[[i]]={{\[Epsilon]mec[[i,1]],\[Epsilon]mec[[i,4]],\[Epsilon]mec[[i,5]]},{\[Epsilon]mec[[i,4]],\[Epsilon]mec[[i,2]],\[Epsilon]mec[[i,6]]},{\[Epsilon]mec[[i,5]],\[Epsilon]mec[[i,6]],\[Epsilon]mec[[i,3]]}},{i,Length[\[Epsilon]mec]}];
\[Epsilon]TenSph =Table[0,{i,Length[\[Epsilon]Ten]}]; 
\[Epsilon]TenDev =Table[0,{i,Length[\[Epsilon]Ten]}]; 
Do[\[Epsilon]TenSph[[i]] = Contr4thAND2ndTen[Jijkl,\[Epsilon]Ten[[i]]],{i,Length[\[Epsilon]mec]}];
Do[\[Epsilon]TenDev[[i]] = Contr4thAND2ndTen[Kijkl,\[Epsilon]Ten[[i]]],{i,Length[\[Epsilon]mec]}];
\[Sigma]Inst = Table[0,{i,Length[\[Epsilon]Ten]}];
Do[\[Sigma]Inst[[i]] = ag[[1]]\[Epsilon]TenSph[[i]] + ag[[2]]\[Epsilon]TenDev[[i]],{i,Length[\[Epsilon]mec]}];
\[Epsilon]vSph = Table[{{0,0,0},{0,0,0},{0,0,0}},{i,Length[\[Epsilon]Ten]},{j,Length[ai[[1]]]}];
\[Epsilon]vDev = Table[{{0,0,0},{0,0,0},{0,0,0}},{i,Length[\[Epsilon]Ten]},{j,Length[ai[[2]]]}];
Do[\[Epsilon]vSph[[i,j]]=Taylor2[\[Epsilon]vSph[[i-1,j]],\[Tau]i[[1,j]],\[Epsilon]TenSph[[i-1]],\[Epsilon]TenSph[[i]],\[CapitalDelta]timeList[[i]]],{i,2,Length[\[Epsilon]Ten]},{j,Length[ai[[1]]]}];
Do[\[Epsilon]vDev[[i,j]]=Taylor2[\[Epsilon]vDev[[i-1,j]],\[Tau]i[[2,j]],\[Epsilon]TenDev[[i-1]],\[Epsilon]TenDev[[i]],\[CapitalDelta]timeList[[i]]],{i,2,Length[\[Epsilon]Ten]},{j,Length[ai[[2]]]}];
\[Sigma]vals =  Table[0,{i,Length[\[Epsilon]Ten]}];
Do[\[Sigma]vals[[i]] =\[Sigma]Inst[[i]] - Sum[ai[[1,j]]*\[Epsilon]vSph[[i,j]],{j,Length[ai[[1]]]}]- Sum[ai[[2,k]]*\[Epsilon]vDev[[i,k]],{k,Length[ai[[2]]]}] ,{i,Length[\[Epsilon]Ten]}];
\[Sigma]vals
]

IsotropicStress2[ag_,ai_, \[Tau]i_,\[Epsilon]mec_, \[CapitalDelta]timeList_]:= Block[
{\[Sigma]Inst,\[Epsilon]Ten,\[Epsilon]TenSph,\[Epsilon]TenDev,\[Epsilon]vSph,\[Epsilon]vDev,\[Sigma]vals,Jijkl,Kijkl},
{Jijkl,Kijkl} = IsoDecomp[];
\[Epsilon]Ten= Table[0,{i,Length[\[Epsilon]mec]}];
Do[\[Epsilon]Ten[[i]]={{\[Epsilon]mec[[i,1]],\[Epsilon]mec[[i,4]],\[Epsilon]mec[[i,5]]},{\[Epsilon]mec[[i,4]],\[Epsilon]mec[[i,2]],\[Epsilon]mec[[i,6]]},{\[Epsilon]mec[[i,5]],\[Epsilon]mec[[i,6]],\[Epsilon]mec[[i,3]]}},{i,Length[\[Epsilon]mec]}];
\[Epsilon]TenSph =Table[0,{i,Length[\[Epsilon]Ten]}]; 
\[Epsilon]TenDev =Table[0,{i,Length[\[Epsilon]Ten]}]; 
Do[\[Epsilon]TenSph[[i]] = Contr4thAND2ndTen[Jijkl,\[Epsilon]Ten[[i]]],{i,Length[\[Epsilon]mec]}];
Do[\[Epsilon]TenDev[[i]] = Contr4thAND2ndTen[Kijkl,\[Epsilon]Ten[[i]]],{i,Length[\[Epsilon]mec]}];
\[Sigma]Inst = Table[0,{i,Length[\[Epsilon]Ten]}];
Do[\[Sigma]Inst[[i]] = ag[[1]]\[Epsilon]TenSph[[i]] + ag[[2]]\[Epsilon]TenDev[[i]],{i,Length[\[Epsilon]mec]}];
\[Epsilon]vSph = Table[{{0,0,0},{0,0,0},{0,0,0}},{i,Length[\[Epsilon]Ten]},{j,Length[ai[[1]]]}];
\[Epsilon]vDev = Table[{{0,0,0},{0,0,0},{0,0,0}},{i,Length[\[Epsilon]Ten]},{j,Length[ai[[2]]]}];
Do[\[Epsilon]vSph[[i,j]]=TaylorT2[\[Epsilon]vSph[[i-1,j]],\[Tau]i[[1,j]],\[Epsilon]TenSph[[i-1]],\[Epsilon]TenSph[[i]],\[CapitalDelta]timeList[[i]]],{i,2,Length[\[Epsilon]Ten]},{j,Length[ai[[1]]]}];
Do[\[Epsilon]vDev[[i,j]]=TaylorT2[\[Epsilon]vDev[[i-1,j]],\[Tau]i[[2,j]],\[Epsilon]TenDev[[i-1]],\[Epsilon]TenDev[[i]],\[CapitalDelta]timeList[[i]]],{i,2,Length[\[Epsilon]Ten]},{j,Length[ai[[2]]]}];
\[Sigma]vals =  Table[0,{i,Length[\[Epsilon]Ten]}];
Do[\[Sigma]vals[[i]] =\[Sigma]Inst[[i]] - Sum[ai[[1,j]]*\[Epsilon]vSph[[i,j]],{j,Length[ai[[1]]]}]- Sum[ai[[2,k]]*\[Epsilon]vDev[[i,k]],{k,Length[ai[[2]]]}] ,{i,Length[\[Epsilon]Ten]}];
\[Sigma]vals
]
(* *)
IsotropicStressEulerImp[ag_,ai_, \[Tau]i_,\[Epsilon]mec_, \[CapitalDelta]timeList_]:= Block[
{\[Sigma]Inst,\[Epsilon]Ten,\[Epsilon]TenSph,\[Epsilon]TenDev,\[Epsilon]vSph,\[Epsilon]vDev,\[Sigma]vals,Jijkl,Kijkl},
{Jijkl,Kijkl} = IsoDecomp[];
\[Epsilon]Ten= Table[0,{i,Length[\[Epsilon]mec]}];
Do[\[Epsilon]Ten[[i]]={{\[Epsilon]mec[[i,1]],\[Epsilon]mec[[i,4]],\[Epsilon]mec[[i,5]]},{\[Epsilon]mec[[i,4]],\[Epsilon]mec[[i,2]],\[Epsilon]mec[[i,6]]},{\[Epsilon]mec[[i,5]],\[Epsilon]mec[[i,6]],\[Epsilon]mec[[i,3]]}},{i,Length[\[Epsilon]mec]}];
\[Epsilon]TenSph =Table[0,{i,Length[\[Epsilon]Ten]}]; 
\[Epsilon]TenDev =Table[0,{i,Length[\[Epsilon]Ten]}]; 
Do[\[Epsilon]TenSph[[i]] = Contr4thAND2ndTen[Jijkl,\[Epsilon]Ten[[i]]],{i,Length[\[Epsilon]mec]}];
Do[\[Epsilon]TenDev[[i]] = Contr4thAND2ndTen[Kijkl,\[Epsilon]Ten[[i]]],{i,Length[\[Epsilon]mec]}];
\[Sigma]Inst = Table[0,{i,Length[\[Epsilon]Ten]}];
Do[\[Sigma]Inst[[i]] = ag[[1]]\[Epsilon]TenSph[[i]] + ag[[2]]\[Epsilon]TenDev[[i]],{i,Length[\[Epsilon]mec]}];
\[Epsilon]vSph = Table[{{0,0,0},{0,0,0},{0,0,0}},{i,Length[\[Epsilon]Ten]},{j,Length[ai[[1]]]}];
\[Epsilon]vDev = Table[{{0,0,0},{0,0,0},{0,0,0}},{i,Length[\[Epsilon]Ten]},{j,Length[ai[[2]]]}];
Do[\[Epsilon]vSph[[i,j]]=EulerImp[\[Epsilon]vSph[[i-1,j]],\[Tau]i[[1,j]],\[Epsilon]TenSph[[i-1]],\[Epsilon]TenSph[[i]],\[CapitalDelta]timeList[[i]]],{i,2,Length[\[Epsilon]Ten]},{j,Length[ai[[1]]]}];
Do[\[Epsilon]vDev[[i,j]]=EulerImp[\[Epsilon]vDev[[i-1,j]],\[Tau]i[[2,j]],\[Epsilon]TenDev[[i-1]],\[Epsilon]TenDev[[i]],\[CapitalDelta]timeList[[i]]],{i,2,Length[\[Epsilon]Ten]},{j,Length[ai[[2]]]}];
\[Sigma]vals =  Table[0,{i,Length[\[Epsilon]Ten]}];
Do[\[Sigma]vals[[i]] =\[Sigma]Inst[[i]] - Sum[ai[[1,j]]*\[Epsilon]vSph[[i,j]],{j,Length[ai[[1]]]}]- Sum[ai[[2,k]]*\[Epsilon]vDev[[i,k]],{k,Length[ai[[2]]]}] ,{i,Length[\[Epsilon]Ten]}];
\[Sigma]vals
]


(* stress function for transverse-isotropic viscoelastic materials *)
(* format of entries 
 n = {{nx,ny,nz}}  being n the unit vector of the symetry axis of the material function
Hill coefs = {2K,l,l,n,2m,2g}
\[CapitalDelta]timeList = {0,\[CapitalDelta]t2,...,\[CapitalDelta]tN}/ \[CapitalDelta]ti = Subscript[t, i]-Subscript[t, i-1], \[CapitalDelta]t1 = t1 - t0; t0=t1=0
\[Epsilon]mec = {\[Epsilon]11,\[Epsilon]22,\[Epsilon]33,\[Epsilon]12,\[Epsilon]13,\[Epsilon]23}
*)

(* Projection of \[Epsilon] onto the Hill's basis *)
\[Epsilon]tilde1[\[Epsilon]Tensor_,n_]:=Block[{nij,Id ,\[CapitalTheta],eps1},nij = Transpose[n] . n;Id = IdentityMatrix[3];
\[CapitalTheta]= Id-nij; eps1 =1/2*\[CapitalTheta]*(Tr[\[Epsilon]Tensor]-(n . \[Epsilon]Tensor . Transpose[n])[[1,1]]);eps1]
\[Epsilon]tilde2[\[Epsilon]Tensor_,n_]:=Block[{nij,Id ,\[CapitalTheta],eps2},nij = Transpose[n] . n;Id = IdentityMatrix[3];
\[CapitalTheta]= Id-nij; eps2 =\[CapitalTheta]*((n . \[Epsilon]Tensor . Transpose[n])[[1,1]]);eps2]
\[Epsilon]tilde3[\[Epsilon]Tensor_,n_]:=Block[{nij,eps3},nij = Transpose[n] . n; eps3 =nij*(Tr[\[Epsilon]Tensor]-(n . \[Epsilon]Tensor . Transpose[n])[[1,1]]);
eps3]
\[Epsilon]tilde4[\[Epsilon]Tensor_,n_]:=Block[{nij,eps4},nij = Transpose[n] . n; eps4 =nij*((n . \[Epsilon]Tensor . Transpose[n])[[1,1]]);
eps4]
\[Epsilon]tilde5[\[Epsilon]Tensor_,n_]:=Block[{nij,Id ,\[CapitalTheta],eps5},nij = Transpose[n] . n;
 eps5 =1/2*(2\[Epsilon]Tensor -2\[Epsilon]Tensor . nij+ 2Transpose[n] . (n*((n . \[Epsilon]Tensor . Transpose[n])[[1,1]])-n . \[Epsilon]Tensor)-2 \[Epsilon]tilde1[\[Epsilon]Tensor,n]);
 eps5]
\[Epsilon]tilde6[\[Epsilon]Tensor_,n_]:=Block[{nij,Id ,\[CapitalTheta],eps6},nij = Transpose[n] . n;
 eps6 = (\[Epsilon]Tensor . nij -((n . \[Epsilon]Tensor . Transpose[n])[[1,1]])*nij)+(nij . \[Epsilon]Tensor-((n . \[Epsilon]Tensor . Transpose[n])[[1,1]])*nij);
 eps6]

(* Stress Function in tensor form *)
TIsoTensorStress[Hillg_, Hilli_,Hill\[Tau]i_,\[Epsilon]mec_,dir_, \[CapitalDelta]timeList_]:= Block[{\[Sigma]Inst,\[Epsilon]Ten,ev,\[Sigma]vals,n},
\[Epsilon]Ten= Table[0,{i,Length[\[Epsilon]mec]}];n = {Table[KroneckerDelta[i,dir],{i,3}]};
Do[\[Epsilon]Ten[[i]]={{\[Epsilon]mec[[i,1]],\[Epsilon]mec[[i,4]],\[Epsilon]mec[[i,5]]},{\[Epsilon]mec[[i,4]],\[Epsilon]mec[[i,2]],\[Epsilon]mec[[i,6]]},{\[Epsilon]mec[[i,5]],\[Epsilon]mec[[i,6]],\[Epsilon]mec[[i,3]]}},{i,Length[\[Epsilon]mec]}];
\[Sigma]Inst = Table[0,{i,Length[\[Epsilon]Ten]}];
Do[\[Sigma]Inst[[i]]=Hillg[[1]]*\[Epsilon]tilde1[\[Epsilon]Ten[[i]],n]+Hillg[[2]]*\[Epsilon]tilde2[\[Epsilon]Ten[[i]],n]+Hillg[[3]]*\[Epsilon]tilde3[\[Epsilon]Ten[[i]],n]+Hillg[[4]]*\[Epsilon]tilde4[\[Epsilon]Ten[[i]],n]+Hillg[[5]]*\[Epsilon]tilde5[\[Epsilon]Ten[[i]],n]+Hillg[[6]]*\[Epsilon]tilde6[\[Epsilon]Ten[[i]],n],{i,Length[\[Epsilon]Ten]}];
ev=  Table[Table[{{0,0,0},{0,0,0},{0,0,0}},{i,Length[Hilli]},{j,Length[Hilli[[i]]]}],{k,Length[\[Epsilon]Ten]}];
Do[ev[[k,1,j]]=Taylor2[ev[[k-1,1,j]],Hill\[Tau]i[[1,j]],\[Epsilon]tilde1[\[Epsilon]Ten[[k-1]],n],\[Epsilon]tilde1[\[Epsilon]Ten[[k]],n],\[CapitalDelta]timeList[[k]]],{j,Length[Hilli[[1]]]},{k,2,Length[\[Epsilon]Ten]}];
Do[ev[[k,2,j]]=Taylor2[ev[[k-1,2,j]],Hill\[Tau]i[[2,j]],\[Epsilon]tilde2[\[Epsilon]Ten[[k-1]],n],\[Epsilon]tilde2[\[Epsilon]Ten[[k]],n],\[CapitalDelta]timeList[[k]]],{j,Length[Hilli[[2]]]},{k,2,Length[\[Epsilon]Ten]}];
Do[ev[[k,3,j]]=Taylor2[ev[[k-1,3,j]],Hill\[Tau]i[[3,j]],\[Epsilon]tilde3[\[Epsilon]Ten[[k-1]],n],\[Epsilon]tilde3[\[Epsilon]Ten[[k]],n],\[CapitalDelta]timeList[[k]]],{j,Length[Hilli[[3]]]},{k,2,Length[\[Epsilon]Ten]}];
Do[ev[[k,4,j]]=Taylor2[ev[[k-1,4,j]],Hill\[Tau]i[[4,j]],\[Epsilon]tilde4[\[Epsilon]Ten[[k-1]],n],\[Epsilon]tilde4[\[Epsilon]Ten[[k]],n],\[CapitalDelta]timeList[[k]]],{j,Length[Hilli[[4]]]},{k,2,Length[\[Epsilon]Ten]}];
Do[ev[[k,5,j]]=Taylor2[ev[[k-1,5,j]],Hill\[Tau]i[[5,j]],\[Epsilon]tilde5[\[Epsilon]Ten[[k-1]],n],\[Epsilon]tilde5[\[Epsilon]Ten[[k]],n],\[CapitalDelta]timeList[[k]]],{j,Length[Hilli[[5]]]},{k,2,Length[\[Epsilon]Ten]}];
Do[ev[[k,6,j]]=Taylor2[ev[[k-1,6,j]],Hill\[Tau]i[[6,j]],\[Epsilon]tilde6[\[Epsilon]Ten[[k-1]],n],\[Epsilon]tilde6[\[Epsilon]Ten[[k]],n],\[CapitalDelta]timeList[[k]]],{j,Length[Hilli[[6]]]},{k,2,Length[\[Epsilon]Ten]}];
\[Sigma]vals =  Table[0,{i,Length[\[Epsilon]Ten]}];
Do[\[Sigma]vals[[k]]=\[Sigma]Inst[[k]] -Sum[Sum[Hilli[[i,j]]*ev[[k,i,j]],{j,Length[Hilli[[i]]]}],{i,Length[Hilli]}],{k,1,Length[\[Epsilon]Ten]}];
\[Sigma]vals
]

TIsoTensorStressEulerImp[Hillg_, Hilli_,Hill\[Tau]i_,\[Epsilon]mec_,dir_, \[CapitalDelta]timeList_]:= Block[{\[Sigma]Inst,\[Epsilon]Ten,ev,\[Sigma]vals,n},
\[Epsilon]Ten= Table[0,{i,Length[\[Epsilon]mec]}];n = {Table[KroneckerDelta[i,dir],{i,3}]};
Do[\[Epsilon]Ten[[i]]={{\[Epsilon]mec[[i,1]],\[Epsilon]mec[[i,4]],\[Epsilon]mec[[i,5]]},{\[Epsilon]mec[[i,4]],\[Epsilon]mec[[i,2]],\[Epsilon]mec[[i,6]]},{\[Epsilon]mec[[i,5]],\[Epsilon]mec[[i,6]],\[Epsilon]mec[[i,3]]}},{i,Length[\[Epsilon]mec]}];
\[Sigma]Inst = Table[0,{i,Length[\[Epsilon]Ten]}];
Do[\[Sigma]Inst[[i]]=Hillg[[1]]*\[Epsilon]tilde1[\[Epsilon]Ten[[i]],n]+Hillg[[2]]*\[Epsilon]tilde2[\[Epsilon]Ten[[i]],n]+Hillg[[3]]*\[Epsilon]tilde3[\[Epsilon]Ten[[i]],n]+Hillg[[4]]*\[Epsilon]tilde4[\[Epsilon]Ten[[i]],n]+Hillg[[5]]*\[Epsilon]tilde5[\[Epsilon]Ten[[i]],n]+Hillg[[6]]*\[Epsilon]tilde6[\[Epsilon]Ten[[i]],n],{i,Length[\[Epsilon]Ten]}];
ev=  Table[Table[{{0,0,0},{0,0,0},{0,0,0}},{i,Length[Hilli]},{j,Length[Hilli[[i]]]}],{k,Length[\[Epsilon]Ten]}];
Do[ev[[k,1,j]]=EulerImp[ev[[k-1,1,j]],Hill\[Tau]i[[1,j]],\[Epsilon]tilde1[\[Epsilon]Ten[[k-1]],n],\[Epsilon]tilde1[\[Epsilon]Ten[[k]],n],\[CapitalDelta]timeList[[k]]],{j,Length[Hilli[[1]]]},{k,2,Length[\[Epsilon]Ten]}];
Do[ev[[k,2,j]]=EulerImp[ev[[k-1,2,j]],Hill\[Tau]i[[2,j]],\[Epsilon]tilde2[\[Epsilon]Ten[[k-1]],n],\[Epsilon]tilde2[\[Epsilon]Ten[[k]],n],\[CapitalDelta]timeList[[k]]],{j,Length[Hilli[[2]]]},{k,2,Length[\[Epsilon]Ten]}];
Do[ev[[k,3,j]]=EulerImp[ev[[k-1,3,j]],Hill\[Tau]i[[3,j]],\[Epsilon]tilde3[\[Epsilon]Ten[[k-1]],n],\[Epsilon]tilde3[\[Epsilon]Ten[[k]],n],\[CapitalDelta]timeList[[k]]],{j,Length[Hilli[[3]]]},{k,2,Length[\[Epsilon]Ten]}];
Do[ev[[k,4,j]]=EulerImp[ev[[k-1,4,j]],Hill\[Tau]i[[4,j]],\[Epsilon]tilde4[\[Epsilon]Ten[[k-1]],n],\[Epsilon]tilde4[\[Epsilon]Ten[[k]],n],\[CapitalDelta]timeList[[k]]],{j,Length[Hilli[[4]]]},{k,2,Length[\[Epsilon]Ten]}];
Do[ev[[k,5,j]]=EulerImp[ev[[k-1,5,j]],Hill\[Tau]i[[5,j]],\[Epsilon]tilde5[\[Epsilon]Ten[[k-1]],n],\[Epsilon]tilde5[\[Epsilon]Ten[[k]],n],\[CapitalDelta]timeList[[k]]],{j,Length[Hilli[[5]]]},{k,2,Length[\[Epsilon]Ten]}];
Do[ev[[k,6,j]]=EulerImp[ev[[k-1,6,j]],Hill\[Tau]i[[6,j]],\[Epsilon]tilde6[\[Epsilon]Ten[[k-1]],n],\[Epsilon]tilde6[\[Epsilon]Ten[[k]],n],\[CapitalDelta]timeList[[k]]],{j,Length[Hilli[[6]]]},{k,2,Length[\[Epsilon]Ten]}];
\[Sigma]vals =  Table[0,{i,Length[\[Epsilon]Ten]}];
Do[\[Sigma]vals[[k]]=\[Sigma]Inst[[k]] -Sum[Sum[Hilli[[i,j]]*ev[[k,i,j]],{j,Length[Hilli[[i]]]}],{i,Length[Hilli]}],{k,1,Length[\[Epsilon]Ten]}];
\[Sigma]vals
]

(*  Thermo-viscoelastic modelling  *)
(* Internal-time technique in incremental deformations scheme *)
(* 1. Trapezoids integration scheme  *)
DeltaXiListTrap[afun_,var_,timelist_,Templist_]:= Block[{\[CapitalDelta]\[Xi]list},
\[CapitalDelta]\[Xi]list = Table[0,{i,Length[timelist]}];
Do[\[CapitalDelta]\[Xi]list[[i]]=(timelist[[i]]-timelist[[i-1]])/2 (1/(afun /. var-> Templist[[i]])  +1/(afun /. var-> Templist[[i-1]]) ),
{i,2,Length[\[CapitalDelta]\[Xi]list]}];\[CapitalDelta]\[Xi]list]
(* 2. midpoint integration scheme  *)
DeltaXiListMid[afun_,var_,timelist_,Templist_]:= Block[{\[CapitalDelta]\[Xi]list,Tempfun,aTmidl},
\[CapitalDelta]\[Xi]list = Table[0,{i,Length[timelist]}];
aTmidl = 0*\[CapitalDelta]\[Xi]list;
Do[aTmidl[[i]] = afun/.{var-> (Templist[[i-1]]+Templist[[i]])*1/2} ;\[CapitalDelta]\[Xi]list[[i]]=
((timelist[[i]]-timelist[[i-1]])/aTmidl[[i]]),
{i,2,Length[\[CapitalDelta]\[Xi]list]}];      
\[CapitalDelta]\[Xi]list]
(* 3. Simpsons 2n integration scheme  *)
DeltaXiListSimp2n[afun_,var_,timelist_,Templist_]:= Block[{\[CapitalDelta]\[Xi]Trlist,\[CapitalDelta]\[Xi]Midlist,\[CapitalDelta]\[Xi]list,\[Xi]list},
\[CapitalDelta]\[Xi]Trlist = DeltaXiListTrap[afun,var,timelist,Templist];
\[CapitalDelta]\[Xi]Midlist =  DeltaXiListMid[afun,var,timelist,Templist]; 
\[CapitalDelta]\[Xi]list =2/3 \[CapitalDelta]\[Xi]Midlist+1/3 \[CapitalDelta]\[Xi]Trlist;
\[Xi]list = 0*\[CapitalDelta]\[Xi]list;
Do[\[Xi]list[[i]] = \[Xi]list[[i-1]]+\[CapitalDelta]\[Xi]list[[i]],{i,2,Length[\[CapitalDelta]\[Xi]list]}];
{\[CapitalDelta]\[Xi]list,\[Xi]list}]
(* 3. Simpsons 2n integration scheme  *)
DeltaXiBoole[afun_,var_,ti_,tip1_,Ti_,Tip1_]:= Block[{\[CapitalDelta]\[Xi],h,tvec,Tvec,Funvec,BooleCoefs}, h =(tip1 - ti) ;
tvec = Table[i,{i,ti,tip1,h/4}];
LinIntFunT[time_]:= Ti + ((Tip1-Ti)/(tip1-ti))*(time-ti);
Tvec = LinIntFunT[tvec];
Funvec=0*Tvec;Do[Funvec[[i]] = (afun^-1) /. var -> Tvec[[i]],{i,Length[Tvec]}];
BooleCoefs = (1/90)*{7,32,12,32,7};
\[CapitalDelta]\[Xi] = h*Total[(Funvec*BooleCoefs)];
\[CapitalDelta]\[Xi]]
DeltaXiBooleList[afun_,var_,timelist_,Templist_]:= Block[{\[CapitalDelta]\[Xi]list,\[Xi]list},
\[CapitalDelta]\[Xi]list = 0*timelist; \[Xi]list = 0*timelist;
Do[\[CapitalDelta]\[Xi]list[[i]] = DeltaXiBoole[afun,var,timelist[[i-1]],timelist[[i]],Templist[[i-1]],Templist[[i]]];\[Xi]list[[i]] = \[Xi]list[[i-1]]+\[CapitalDelta]\[Xi]list[[i]],{i,2,Length[timelist]}]
;{\[CapitalDelta]\[Xi]list,\[Xi]list}]
(* 4. Integration interpolation *)
DeltaXiIntInt[afun_,var_,timelist_,Templist_]:=Block[{afuninvvals,funInt,\[Xi]vals},
afuninvvals = (afun^-1) /. var -> Templist;
\[Xi]vals = Table[0,{i,Length[timelist]}];
funInt = Interpolation[Thread[{timelist,afuninvvals}]];
Do[\[Xi]vals[[i]]=NIntegrate[funInt[u],{u,0,timelist[[i]]}],{i,2,Length[\[Xi]vals]}]
;]


End[]
EndPackage[]




(* ::Package:: *)

(* ::Section:: *)
(*Package for Analytical micro-mechanics*)


(* ::Text:: *)
(*This  package carries with Mean-field computations in thermo-viscoelasticity, focusing in transverse-isotropic media*)


(* ::Input:: *)
(*(* Name: MicroMeca *)*)
(*(* Author: Camilo Andr\[EAcute]s SUAREZ AFANADOR *)*)
(*(* Date : May 3,2021 *)*)
(* *)


BeginPackage["MicroMeca`"]

MeanFieldSchemes::usage = 
"MeanFieldSchemes[\"IdScheme\", \"MatrixProps\", \"InclusionProps\", \"InclusionForm\", \"VolFract\", \"AspectRatSpatDist\"]
Gives \"\!\(\*SubscriptBox[\(aHill\), \(eff\)]\)\" vector, given: \
\"IdScheme\" (1:Dil, 2:M&T, 3:UpH&S, 4:Liel, 5:IDD, 6:LielensModif),\"MatrixProps\" \
({\!\(\*SubscriptBox[\(\[Kappa]\), \(m\)]\),\!\(\*SubscriptBox[\(\[Mu]\), \(m\)]\)}),\"InclusionProps\" \
({\!\(\*SubscriptBox[\(\[Kappa]\), \(i\)]\),\!\(\*SubscriptBox[\(\[Mu]\), \(i\)]\)}),
\"InclusionForm\" (-3:Sphere, -2:Cylinder, -1:Penny, w:L/D), and the \
optional \"AspectRatSpatDist\" (same as inclusion by default)"

MeanFieldSchemes\[Alpha]::usage = 
"MeanFieldSchemes[\"IdScheme\", \"MatrixProps\", \"InclusionProps\", \"InclusionForm\", \"VolFract\", \"AspectRatSpatDist\"]
Gives \"\!\(\*SubscriptBox[\(aHill\), \(eff\)]\)\" vector and the 3x3 matrix of \!\(\*SubscriptBox[\(\[Alpha]\), \(eff\)]\), given: \
\"IdScheme\" (1:Dil, 2:M&T, 3:UpH&S, 4:Liel, 5:IDD),\"MatrixProps\" \
({\!\(\*SubscriptBox[\(\[Kappa]\), \(m\)]\),\!\(\*SubscriptBox[\(\[Mu]\), \(m\)]\),\!\(\*SubscriptBox[\(\[Alpha]\), \(m\)]\)}),\"InclusionProps\" \
({\!\(\*SubscriptBox[\(\[Kappa]\), \(i\)]\),\!\(\*SubscriptBox[\(\[Mu]\), \(i\)]\),\!\(\*SubscriptBox[\(\[Alpha]\), \(i\)]\)}),
\"InclusionForm\" (-3:Sphere, -2:Cylinder, -1:Penny, w:L/D), and the \
optional \"AspectRatSpatDist\" (same as inclusion by default)"

AorW::usage=
"[\[Kappa]o,\[Mu]o,\[Kappa]r,\[Mu]r,wo,cw,wmin,wmax,Nw]"

MeanFieldSchemesATW::usage = 
"MeanFieldSchemes[\"IdScheme\", \"MatrixProps\", \"InclusionProps\", \"formD\", \"MicroPars\", \"Nw\", \"Mode\"]
Gives \"\!\(\*SubscriptBox[\(L\), \(eff\)]\)\" vector, given: \
\"IdScheme\" (1:Dil, 2:M&T, 3:UpH&S, 4:Liel, 5:IDD, 6:LielensModif),\"MatrixProps\" \
({\!\(\*SubscriptBox[\(\[Kappa]\), \(m\)]\),\!\(\*SubscriptBox[\(\[Mu]\), \(m\)]\)}),\"InclusionProps\" \
({\!\(\*SubscriptBox[\(\[Kappa]\), \(i\)]\),\!\(\*SubscriptBox[\(\[Mu]\), \(i\)]\)}),
\"InclusionForm\" (-3:Sphere, -2:Cylinder, -1:Penny, w:L/D), and the \
optional \"AspectRatSpatDist\" (same as inclusion by default)"

MeanFieldSchemes\[Alpha]Freqs::usage = 
"MeanFieldSchemesFreqs[IdScheme_, MatrixPropsIn_, InclusionPropsIn_,InclusionForm_,VolFract_,AspectRatSpatDist_:0,var_Symbol,
expminMec_,expmaxMec_,expminTh_,expmaxTh_, di_, nPartsMec_, nPartsTh_, Tolfit_:0.5]"

PHillCoef::usage = 
"PHillCoef[\[Kappa]o_,\[Mu]o_,form_]"

Aor::usage = 
"Aor [\[Kappa]o_,\[Mu]o_,\[Kappa]r_,\[Mu]r_,form_]"

ArHS::usage = 
"ArHS[\[Kappa]o_,\[Mu]o_,\[Kappa]m_,\[Mu]m_,\[Kappa]f_,\[Mu]f_,form_,c_]"

EshCyl::usage = 
"EshCyl[\[Kappa]_,\[Mu]_]"

PHillCoefExplForms::usage = 
"PHillCoefExplForms[\[Kappa]o_,\[Mu]o_,form_]"

HillfunMods::usage = 
"HillfunMods[\"ListOfFunc\", \"var\", \"RelaxQ\"]
Gives a list of Prony parameters and approximation errors ({mg,mi,\[Tau]i,errors}) given: \
ListOfFunc, a list of viscoelastic functions in Laplace-Carson Space. \
var, The argument of the Laplace-Carson transform. \
RelaxQ, boolean variable for Relaxation or Creep functions identification"

MeanFieldSchemes\[Alpha]ATW::usage = 
"MeanFieldSchemes\[Alpha]2PATW[IdScheme, MatrixProps, InclusionProps,InclusionForm,VolFract,AspectRatSpatDist, MicroPars, Mode]
gives the couple {CHill,\[Alpha]Tildevec} for a given:
IdScheme = 1,2,3,4 or 5
MatrixProps = {Subscript[\[Kappa], m], Subscript[\[Mu], m], Subscript[\[Alpha], m]}
InclusionProps = {Subscript[\[Kappa], i], Subscript[\[Mu], i], Subscript[\[Alpha], i]}
InclusionForm = -3,-2,-1, w>0
AspectRatSpatDist = 0 if not IDD, \!\(\*SubscriptBox[\(w\), \(\(D\)\(\\\ \)\)]\)>0 if IDD
MicroPars = {cfTh, cfcomp, Ninc, ws, ns, wmArr, cWTh, woWTh, mATTh, wmFit, cFit, woFit,mATFit}
Mode 1 if Theoretical, 2 if Fitted, 3 if Deterministic"

MeanFieldSchemes\[Alpha]FreqsATW::usage = 
"MeanFieldSchemesFreqsATW[IdScheme, MatrixPropsIn, InclusionPropsIn,InclusionForm,VolFract,AspectRatSpatDist, MicroPars,
Mode,var, expminMec,expmaxMec,expminTh,expmaxTh, di, nPartsMec, nPartsTh, Tolfit==0.5] with output: {Mg,Mi,\[Tau]Mi,ErrsFit}"

EshCyl::usage = 
"EshCyl[\[Kappa]_,\[Mu]_]"


Begin["`Private`"]
SetDirectory[NotebookDirectory[]];

Needs["TensorCalc`"]
Needs["MicroStats`"]
Needs["ViscoIdentification`"]
(* --------------------------------------------------------------------------------------------*)

(* 1.Eshelby's tensor for prolate isotropic spheroids (Eshelby 1957, Mura 1987) *)
\[CapitalDelta] = Product[Sqrt[a[i]^2+u],{i,3}];
Coef1 = 1/(8\[Pi](1-\[Nu]));
Coef2 = (1-2\[Nu])/(8\[Pi](1-\[Nu]));
(* Prolate spheroid definitions w=aspect ratio with symmetry axis = Overscript[Subscript[e, 3], \[RightVector]] *)
a[2]=a[1];
a[3]= w a[1];
(* Values to integrate in I1 and I2 *)
integrand1 = Table[(2\[Pi] a[1] a[2] a[3])/((a[i]^2 +u)\[CapitalDelta]),{i,3}];
integrand2 = Table[(2\[Pi] a[1] a[2] a[3])/((a[i]^2 +u)(a[j]^2 +u)\[CapitalDelta]),{i,3},{j,3}];
Ini = Assuming[a[1]>0&&w>0,Integrate[integrand1,{u,0,\[Infinity]}]];
Inij = Assuming[a[1]>0&&w>0,Integrate[integrand2,{u,0,\[Infinity]}]];
(* Components of S^E *)
Sijkl = Table[0,{i,1,3},{j,1,3},{k,1,3},{l,1,3}];
(* Eshelby's tensor components *)
Do[Sijkl[[i,i,i,i]]=3Coef1 a[i]^2 Inij[[i,i]] + Coef2 Ini[[i]],{i,3}]
Do[ If[i!=j,Sijkl[[i,i,j,j]] = Coef1 a[j]^2 Inij[[i,j]] - Coef2 Ini[[i]]],{i,3},{j,3}]
Do[ If[i!=j,Sijkl[[i,j,i,j]] = Coef1/2 (a[i]^2+a[j]^2) Inij[[i,j]] + Coef2/2 (Ini[[i]]+ Ini[[j]])],{i,3},{j,3}]
Do[ If[i!=j,Sijkl[[i,j,j,i]] =Sijkl[[i,j,i,j]]//Simplify],{i,3},{j,3}]
(* Eshelby's Tensors for an isotropic media *)
EshSphere[\[Kappa]_,\[Mu]_]:=Block[{nu,S},nu=(3\[Kappa]-2\[Mu])/(2(3\[Kappa] + \[Mu])); S=Limit[Sijkl,w-> 1]/. {\[Nu]-> nu};S]
EshCyl[\[Kappa]_,\[Mu]_]:=Block[{nu,S},nu=(3\[Kappa]-2\[Mu])/(2(3\[Kappa] + \[Mu])); S=Limit[Sijkl,w-> \[Infinity]]/. {\[Nu]-> nu};S]
EshPenny[\[Kappa]_,\[Mu]_]:=Block[{nu,S},nu=(3\[Kappa]-2\[Mu])/(2(3\[Kappa] + \[Mu])); S=Limit[Sijkl,w-> 0]/. {\[Nu]-> nu};S]
EshSpheroid[\[Kappa]_,\[Mu]_,wval_]:=Block[{nu,S},nu=(3\[Kappa]-2\[Mu])/(2(3\[Kappa] + \[Mu])); S=Limit[Sijkl,w-> wval]/. {\[Nu]-> nu};S]

(* 2.Hill's Tensor Subscript[P, 0]=S^E(Subscript[L, 0])^-1.
Identificator for inclusion geometry:
form=Sphere(-3)
form=Cylinder(-2)
form=Penny(-1)
form=w=Spheroid(w>0) *)
(*PHillCoef[\[Kappa]o_,\[Mu]o_,form_] := Block[{SeCoef ,CoCoefInv,PHill},
If[form==-3,SeCoef  = From4thTenToHillCoef[EshSphere[\[Kappa]o,\[Mu]o],3]];
If[form==-2,SeCoef = From4thTenToHillCoef[EshCyl[\[Kappa]o,\[Mu]o],3]];
If[form==-1,SeCoef  = From4thTenToHillCoef[EshPenny[\[Kappa]o,\[Mu]o],3]];
If[form>0,SeCoef  = From4thTenToHillCoef[EshSpheroid[\[Kappa]o,\[Mu]o,form],3]];
CoCoefInv=Inv4thHillCoef[FromIsoCoef2HillCoef[3\[Kappa]o,2\[Mu]o]];
PHill = Contr4thTenHillBasis[SeCoef,CoCoefInv];
PHill]*)

PHillCoef[\[Kappa]o_,\[Mu]o_,form_] := Block[{SeCoef ,CoCoefInv,PHill,w},
If[form==-3,SeCoef  = From4thTenToHillCoef[EshSphere[\[Kappa]o,\[Mu]o],3];CoCoefInv=Inv4thHillCoef[FromIsoCoef2HillCoef[3\[Kappa]o,2\[Mu]o]];PHill = Contr4thTenHillBasis[SeCoef,CoCoefInv];];
If[form==-2,SeCoef = From4thTenToHillCoef[EshCyl[\[Kappa]o,\[Mu]o],3];CoCoefInv=Inv4thHillCoef[FromIsoCoef2HillCoef[3\[Kappa]o,2\[Mu]o]];PHill = Contr4thTenHillBasis[SeCoef,CoCoefInv];];
If[form==-1,SeCoef  = From4thTenToHillCoef[EshPenny[\[Kappa]o,\[Mu]o],3];CoCoefInv=Inv4thHillCoef[FromIsoCoef2HillCoef[3\[Kappa]o,2\[Mu]o]];PHill = Contr4thTenHillBasis[SeCoef,CoCoefInv];];
If[form>1,w = form;PHill= N[{(w (3 w Sqrt[-1+w^2] (-3 \[Kappa]o+(-3+2 w^2) \[Mu]o)+((3+6 w^2) \[Kappa]o+(7-4 w^2) \[Mu]o) ArcCosh[w]))/(4 (-1+w^2)^(5/2) \[Mu]o (3 \[Kappa]o+4 \[Mu]o)),-((w (3 \[Kappa]o+\[Mu]o) (-3 w Sqrt[-1+w^2]+ArcCosh[w]+2 w^2 ArcCosh[w]))/(4 (-1+w^2)^(5/2) \[Mu]o (3 \[Kappa]o+4 \[Mu]o))),-((w (3 \[Kappa]o+\[Mu]o) (-3 w Sqrt[-1+w^2]+ArcCosh[w]+2 w^2 ArcCosh[w]))/(4 (-1+w^2)^(5/2) \[Mu]o (3 \[Kappa]o+4 \[Mu]o))),(-3 Sqrt[-1+w^2] (-2 \[Mu]o+3 w^2 (\[Kappa]o+\[Mu]o))+w ((3+6 w^2) \[Kappa]o+(-5+8 w^2) \[Mu]o) ArcCosh[w])/(2 (-1+w^2)^(5/2) \[Mu]o (3 \[Kappa]o+4 \[Mu]o)),(w (w Sqrt[-1+w^2] (3 (-5+2 w^2) \[Kappa]o+(-17+14 w^2) \[Mu]o)+3 (3 \[Kappa]o+(5-4 w^2) \[Mu]o) ArcCosh[w]))/(8 (-1+w^2)^(5/2) \[Mu]o (3 \[Kappa]o+4 \[Mu]o)),(Sqrt[-1+w^2] (3 (2+3 w^2+w^4) \[Kappa]o+2 (4-3 w^2+2 w^4) \[Mu]o)-3 w (3 (1+w^2) \[Kappa]o+2 \[Mu]o) ArcCosh[w])/(4 (-1+w^2)^(5/2) \[Mu]o (3 \[Kappa]o+4 \[Mu]o))}]];
If[0<form<1,w = form;PHill=N[{(w (3 w Sqrt[1-w^2] (-3 \[Kappa]o+(-3+2 w^2) \[Mu]o)+((3+6 w^2) \[Kappa]o+(7-4 w^2) \[Mu]o) ArcCos[w]))/(4 (1-w^2)^(5/2) \[Mu]o (3 \[Kappa]o+4 \[Mu]o)),-((w (3 \[Kappa]o+\[Mu]o) (-3 w Sqrt[1-w^2]+ArcCos[w]+2 w^2 ArcCos[w]))/(4 (1-w^2)^(5/2) \[Mu]o (3 \[Kappa]o+4 \[Mu]o))),-((w (3 \[Kappa]o+\[Mu]o) (-3 w Sqrt[1-w^2]+ArcCos[w]+2 w^2 ArcCos[w]))/(4 (1-w^2)^(5/2) \[Mu]o (3 \[Kappa]o+4 \[Mu]o))),(-3 Sqrt[1-w^2] (-2 \[Mu]o+3 w^2 (\[Kappa]o+\[Mu]o))+w ((3+6 w^2) \[Kappa]o+(-5+8 w^2) \[Mu]o) ArcCos[w])/(2 (1-w^2)^(5/2) \[Mu]o (3 \[Kappa]o+4 \[Mu]o)),(w (w Sqrt[1-w^2] (3 (-5+2 w^2) \[Kappa]o+(-17+14 w^2) \[Mu]o)+3 (3 \[Kappa]o+(5-4 w^2) \[Mu]o) ArcCos[w]))/(8 (1-w^2)^(5/2) \[Mu]o (3 \[Kappa]o+4 \[Mu]o)),(Sqrt[1-w^2] (3 (2+3 w^2+w^4) \[Kappa]o+2 (4-3 w^2+2 w^4) \[Mu]o)-3 w (3 (1+w^2) \[Kappa]o+2 \[Mu]o) ArcCos[w])/(4 (1-w^2)^(5/2) \[Mu]o (3 \[Kappa]o+4 \[Mu]o))}]];
PHill]
(* --------------------------------------------------------------------------------------------*)

(* 
Functions to compute Localization tensors in biphasic composites (A) with isotropic phases.
All computations are performed by means of the Hill's basis representation
*)

(* Fourth order identity tensor Subscript[I, 4] in Hill's basis  *)
I4 = FromIsoCoef2HillCoef[1,1];

(* 
3. Localization tensor for non-interacting inclusions with isotropic phases in Hill's basis
   for a given reference infinite medium Subscript[L, 0] and an inclusionar phase (L^(r)).     
*)
Aor[\[Kappa]o_,\[Mu]o_,\[Kappa]r_,\[Mu]r_,form_]:= Block[{Po,Ar},Po = PHillCoef[\[Kappa]o,\[Mu]o,form];
Ar = Inv4thHillCoef[I4 + Contr4thTenHillBasis[Po,(FromIsoCoef2HillCoef[3\[Kappa]r,2\[Mu]r]-FromIsoCoef2HillCoef[3\[Kappa]o,2\[Mu]o])]];
Ar]

(* 
4. Hashin-Shtrikman Scheme Localization tensors for biphasic medium with isotropic constituents (Bornert 2001) 
*)
ArHS[\[Kappa]o_,\[Mu]o_,\[Kappa]m_,\[Mu]m_,\[Kappa]f_,\[Mu]f_,form_,c_]:=Block[{Aom,Aof,Amean,AHS}, 
Aom = Aor[\[Kappa]o,\[Mu]o,\[Kappa]m,\[Mu]m,form]; Aof= Aor[\[Kappa]o,\[Mu]o,\[Kappa]f,\[Mu]f,form];
Amean = (1-c)Aom + c Aof;
AHS = Contr4thTenHillBasis[Aof,Inv4thHillCoef[Amean]];
AHS]

(*
5. Mori-Tanaka Scheme Localization tensors for biphasic medium with isotropic constituents (Bornert 2001) 
*)
ArMT[\[Kappa]m_,\[Mu]m_,\[Kappa]f_,\[Mu]f_,form_,c_]:=Block[{Ao,AMT},
Ao = Aor[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,form];
AMT = Contr4thTenHillBasis[Ao,Inv4thHillCoef[(1-c)I4 + c Ao]];
;AMT]

(*
6. Lielens Scheme Localization tensors for biphasic medium with isotropic constituents (Bornert 2001) 
*)
ArLiel[\[Kappa]m_,\[Mu]m_,\[Kappa]f_,\[Mu]f_,form_,c_]:=Block[{Ao,Aoinv,LielPar,Abar,ALiel},
LielPar =  0.5*(c+c^2);
Ao = Aor[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,form];
Aoinv = Aor[\[Kappa]f,\[Mu]f,\[Kappa]m,\[Mu]m,form];
Abar=Inv4thHillCoef[(1-LielPar)Inv4thHillCoef[Ao]+LielPar Inv4thHillCoef[Aoinv]];
ALiel = Contr4thTenHillBasis[Abar,Inv4thHillCoef[(1-c)I4 + c Abar]];
;ALiel]

(*
7. Modified Lielens Scheme Localization tensors for biphasic medium with isotropic constituents (Bornert 2001) 
*)
ArLielMod[\[Kappa]m_,\[Mu]m_,\[Kappa]f_,\[Mu]f_,form_,c_]:=Block[{ALHS,AUHS,LielPar,ALielMod},
LielPar =  0.5*(c+c^2);
ALHS = ArHS[\[Kappa]m,\[Mu]m,\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,form,c];
AUHS = ArHS[\[Kappa]f,\[Mu]f,\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,form,c];
ALielMod = Inv4thHillCoef[(1-LielPar)Inv4thHillCoef[ALHS]+LielPar Inv4thHillCoef[AUHS]];
;ALielMod]

(*
8. IDD Scheme Localization tensors for biphasic medium with isotropic constituents (Bornert 2001) 
*)
ArIDD[\[Kappa]m_,\[Mu]m_,\[Kappa]f_,\[Mu]f_,form_,c_,formD_]:=Block[{Ao,PoD,Lf,Lm,AIDD},
Ao = Aor[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,form];
PoD = PHillCoef[\[Kappa]m,\[Mu]m,formD];
Lf =  FromIsoCoef2HillCoef[3\[Kappa]f,2\[Mu]f];
Lm = FromIsoCoef2HillCoef[3\[Kappa]m,2\[Mu]m];
AIDD = Contr4thTenHillBasis[Ao,Inv4thHillCoef[I4 - c Contr4thTenHillBasis[PoD,Contr4thTenHillBasis[(Lf-Lm),Ao]]]];
;AIDD]

(* --------------------------------------------------------------------------------------------*)

(* Functions to compute effective behavior (Leff) in Hill's basis representation

IdScheme values for identification of the scheme;
1 = Diluted;
2 = Mori-Tanaka;
3 = Upper-bound Hashin-Shtrikman;
4 = Lielens Approximation;
5 = IDD Zheng & Dung;
6 = Modified Lielens;

Inclusion form;
-3 \[Equal] Spherical inclusion;
-2 \[Equal] Cylindrical inclusion;
-1 \[Equal] Penny-Shaped inclusion;
w>0 \[Equal] Ellipsoidal inclusion with aspect ratio w;

Optional parameter for characterizing spatial distribution using a Hill-Tensor Subscript[P, 0];
formD \[Equal]  Ellipsoidal distribution of aspect ratio Subscript[w, 0] 
*)

(* 
9. Function to compute effective elastic or viscoelastic behavior (Leff)

Format for entries;
MatrixProps = {Subscript[\[Kappa], m], Subscript[\[Mu], m]};
InclusionProps = {Subscript[\[Kappa], i], Subscript[\[Mu], i]};
*)
MeanFieldSchemes[IdScheme_, MatrixProps_, InclusionProps_,form_,c_,formD_:0]:=
Block[{\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,Lf,Lm,Af,Leff},
{\[Kappa]m,\[Mu]m}=MatrixProps;
{\[Kappa]f,\[Mu]f}=InclusionProps;
Lf = FromIsoCoef2HillCoef[3\[Kappa]f,2\[Mu]f];
Lm = FromIsoCoef2HillCoef[3\[Kappa]m,2\[Mu]m];
If[IdScheme==1,Af=Aor[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,form]];
If[IdScheme==2,Af=ArMT[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,form,c]];
If[IdScheme==3,Af=ArHS[\[Kappa]f,\[Mu]f,\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,form,c]];
If[IdScheme==4,Af=ArLiel[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,form,c]];
If[IdScheme==6,Af=ArLielMod[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,form,c]];
If[IdScheme==5,Af=ArIDD[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,form,c,formD]];
Leff =Lm + c Contr4thTenHillBasis[(Lf-Lm),Af ];
Leff]

(* 
10. Function to compute effective thermoelastic or thermo-viscoelastic behavior (Leff,\[Alpha]eff)

Format for entries:
MatrixProps = {Subscript[\[Kappa], m], Subscript[\[Mu], m], Subscript[\[Alpha], m]};
InclusionProps = {Subscript[\[Kappa], i], Subscript[\[Mu], i], Subscript[\[Alpha], i]};
*)
MeanFieldSchemes\[Alpha][IdScheme_, MatrixProps_, InclusionProps_,form_,c_,formD_:0]:=Block[{\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,Lf,Lm,Sf,Sm,\[Alpha]f,\[Alpha]m
,\[Alpha]mean,Smean,Leff,PRossen,\[Alpha]eff,\[Alpha]effvec},
{\[Kappa]m,\[Mu]m,\[Alpha]m}=MatrixProps;
{\[Kappa]f,\[Mu]f,\[Alpha]f}=InclusionProps;
Lf =  FromIsoCoef2HillCoef[3\[Kappa]f,2\[Mu]f];
Lm = FromIsoCoef2HillCoef[3\[Kappa]m,2\[Mu]m];
Sf =  FromIsoCoef2HillCoef[(3\[Kappa]f)^-1,(2\[Mu]f)^-1];
Sm = FromIsoCoef2HillCoef[(3\[Kappa]m)^-1,(2\[Mu]m)^-1];
\[Alpha]f = InclusionProps[[3]]*IdentityMatrix[3];
\[Alpha]m = MatrixProps[[3]]*IdentityMatrix[3];
\[Alpha]mean =  (1-c)\[Alpha]m + c \[Alpha]f;
PRossen= FromIsoCoef2HillCoef[3(1/\[Kappa]m-1/\[Kappa]f )^-1 , 2(1/\[Mu]m-1/\[Mu]f )^-1];
Smean =(1-c)Sm  + c Sf;
Leff = MeanFieldSchemes[IdScheme, {\[Kappa]m,\[Mu]m}, {\[Kappa]f,\[Mu]f},form,c,formD];
\[Alpha]eff = \[Alpha]mean + Contr4thAND2ndTen[TIso4thTensor[Contr4thTenHillBasis[(Inv4thHillCoef[Leff]-Smean),PRossen],3] ,(\[Alpha]m-\[Alpha]f)];
\[Alpha]effvec = {\[Alpha]eff[[1,1]],\[Alpha]eff[[3,3]]};
{Leff,\[Alpha]effvec}]

(* --------------------------------------------------------------------------------------------*)
(* 
Functions to deal with distributions of orientations and lengths 
*)

(* 
11. Function for Projection over the filament axis of an Elasticity tensor
*)
ProjecTIsoe3[Tens_]:= Block[{Exx, Ezz, Exz, Proxx, Prozz, Proxz,
K2,m2,l2,l1,n,g2},
Exx = {{1,0,0},{0,0,0},{0,0,0}};
Ezz = {{0,0,0},{0,0,0},{0,0,1}};
Exz = {{0,0,1},{0,0,0},{1,0,0}};
Proxx = Contr4thAND2ndTen[Tens,Exx];
Prozz = Contr4thAND2ndTen[Tens,Ezz];
Proxz = Contr4thAND2ndTen[Tens,Exz];
K2 = Proxx[[1,1]]+Proxx[[2,2]];  
m2 = Proxx[[1,1]]-Proxx[[2,2]];
l2 = Proxx[[3,3]];
l1 = Prozz[[1,1]];
n = Prozz[[3,3]];
g2 = Proxz[[1,3]];
{K2,l2,l1,n,m2,g2}]

(*
Mean over the distribution of lengths of the localization tensors Hill's coefficients
*)

(*
12. Mean "diluted" localization tensor
*)
AorW[\[Kappa]o_,\[Mu]o_,\[Kappa]r_,\[Mu]r_,wo_,cw_,wmin_,wmax_,Nw_]:= Block[{wms,fws,wtot,weigs,AoW},
{wms,fws,wtot,weigs} = WeibullArray[wo,cw,wmin,wmax,Nw];AoW = Table[0.,{i,6}];
Do[AoW+=weigs[[j]]*Aor[\[Kappa]o,\[Mu]o,\[Kappa]r,\[Mu]r,wms[[j]]],{j,Length[weigs]}]
;AoW]

(*
13. Mean Mori-Tanaka localization tensor
*)
ArMTW[\[Kappa]m_,\[Mu]m_,\[Kappa]f_,\[Mu]f_,c_,wo_,cw_,wmin_,wmax_,Nw_]:= Block[{AoW,AMTW},
AoW = AorW[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,wo,cw,wmin,wmax,Nw];
AMTW = Contr4thTenHillBasis[AoW,Inv4thHillCoef[(1-c)I4 + c AoW]];
;AMTW]

(*
14. Mean Lielens localization tensor
*)
ArLielW[\[Kappa]m_,\[Mu]m_,\[Kappa]f_,\[Mu]f_,c_,wo_,cw_,wmin_,wmax_,Nw_]:= Block[{wms,fws,wtot,weigs,ALielW},
{wms,fws,wtot,weigs} = WeibullArray[wo,cw,wmin,wmax,Nw];ALielW = Table[0.,{i,6}];
Do[ALielW+=weigs[[j]]*ArLiel[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,wms[[j]],c],{j,Length[weigs]}]
;ALielW]

(*
15. Mean modified Lielens localization tensor
*)
ArLielModW[\[Kappa]m_,\[Mu]m_,\[Kappa]f_,\[Mu]f_,c_,wo_,cw_,wmin_,wmax_,Nw_]:= Block[
{Fo,Foi,wms,fws,wtot,weigs,\[Eta],AojM1,AoijM1,ALielModW},
Fo = (1-c)I4 + c AorW[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,wo,cw,wmin,wmax,Nw];
Foi = (1-c)I4 + c AorW[\[Kappa]f,\[Mu]f,\[Kappa]m,\[Mu]m,wo,cw,wmin,wmax,Nw];
{wms,fws,wtot,weigs} = WeibullArray[wo,cw,wmin,wmax,Nw];
\[Eta] =  0.5*(c+c^2); ALielModW = Table[0.,{i,6}];
Do[
AojM1 = Inv4thHillCoef[Aor[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,wms[[j]]]]; 
AoijM1 = Inv4thHillCoef[Aor[\[Kappa]f,\[Mu]f,\[Kappa]m,\[Mu]m,wms[[j]]]];
ALielModW+=weigs[[j]]*Inv4thHillCoef[(1-\[Eta])Contr4thTenHillBasis[Fo,AojM1] + \[Eta] Contr4thTenHillBasis[Foi,AoijM1]]
,{j,Length[weigs]}];
;ALielModW]

(*
16. Mean PCW localization tensor
*)
ArPCWW[\[Kappa]m_,\[Mu]m_,\[Kappa]f_,\[Mu]f_,c_,formD_,wo_,cw_,wmin_,wmax_,Nw_]:= Block[{AoW,PoD,Lf,Lm,APCWW},
AoW = AorW[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,wo,cw,wmin,wmax,Nw];
PoD = PHillCoef[\[Kappa]m,\[Mu]m,formD*wo];
Lf =  Rationalize[FromIsoCoef2HillCoef[3\[Kappa]f,2\[Mu]f],2.*10^-50];
Lm = Rationalize[FromIsoCoef2HillCoef[3\[Kappa]m,2\[Mu]m],2.*10^-50];
APCWW = Contr4thTenHillBasis[AoW,Inv4thHillCoef[I4 - c Contr4thTenHillBasis[Contr4thTenHillBasis[PoD,(Lf-Lm)],AoW]]];
;APCWW]

(*
17. Mean IDD localization tensor
*)
ArIDDW[\[Kappa]m_,\[Mu]m_,\[Kappa]f_,\[Mu]f_,c_,formD_,wo_,cw_,wmin_,wmax_,Nw_]:= Block[
{AoW,Lf,Lm,wms,fws,wtot,weigs,Poj,Aoj,MeanIn,AIDDW},
AoW = AorW[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,wo,cw,wmin,wmax,Nw];
Lf =  FromIsoCoef2HillCoef[3\[Kappa]f,2\[Mu]f];
Lm = FromIsoCoef2HillCoef[3\[Kappa]m,2\[Mu]m];
{wms,fws,wtot,weigs} = WeibullArray[wo,cw,wmin,wmax,Nw];MeanIn=Table[0.,{i,6}];
Do[
Poj = PHillCoef[\[Kappa]m,\[Mu]m,formD*wms[[j]]];
Aoj = Aor[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,wms[[j]]];
MeanIn+=weigs[[j]]*Contr4thTenHillBasis[Contr4thTenHillBasis[Poj,(Lf-Lm)],Aoj];
,{j,Length[weigs]}];
AIDDW = Contr4thTenHillBasis[AoW,Inv4thHillCoef[I4 - c MeanIn]];
;AIDDW]

(* --------------------------------------------------------------------------------------------*)

(* Functions to compute effective behavior (Leff) in Hill's basis representation
   for a distribution of lengths (f(w)) and orientations (fn(n))

IdScheme values for identification of the scheme;
1 = Diluted;
2 = Mori-Tanaka;
4 = Lielens Approximation;
5 = IDD Zheng & Dung;
6 = Modified Lielens;
7 = PCW

form:
-3 \[Equal] Spherical inclusion;
-2 \[Equal] Cylindrical inclusion;
-1 \[Equal] Penny-Shaped inclusion;
w>0 \[Equal] Ellipsoidal inclusion with aspect ratio w;

Optional parameter for characterizing spatial distribution using a Hill-Tensor Subscript[P, 0];
formD \[Equal]  scaling factor of the spatial distribution form 

Mode:
Theoretical Advani & Tucker Law = 1
Fitted Advani & Tucker Law = 2
Deterministic approach = 3 
*)

(* 
18. Function to compute effective elastic or viscoelastic behavior (Leff)

Format for entries;
MatrixProps = {Subscript[\[Kappa], m], Subscript[\[Mu], m]};
InclusionProps = {Subscript[\[Kappa], i], Subscript[\[Mu], i]};
*)
MeanFieldSchemesATW[IdScheme_, MatrixProps_, InclusionProps_,formD_:0, MicroPars_,Nw_, Mode_]:=
Block[{cfTh, cfcomp, Ninc, ws, ns, wmArr, cWTh, woTh, mATTh, wmFit, cFit, woFit,
mATFit,cw,wo,mAT,a2,a4,Hillmean,wmin,wmax,\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,Lf,Lm,Af,Leff,LeffTen,LeffProj},
{\[Kappa]m,\[Mu]m}=MatrixProps;
{\[Kappa]f,\[Mu]f}=InclusionProps;
{cfTh, cfcomp, Ninc, ws, ns, wmArr, cWTh, woTh, mATTh, wmFit, cFit, woFit, mATFit}=MicroPars;
(* Mode setting *)
If[Mode==1, c=cfTh;cw=cWTh;wo=woTh;mAT=mATTh;{a2,a4}=OriTenATAxSym[mAT];Hillmean=MeanHillBasis[a2,a4]];
If[Mode==2, c=cfcomp;cw=cFit;wo=woFit;mAT=mATFit;{a2,a4}=OriTenATAxSym[mAT];Hillmean=MeanHillBasis[a2,a4]];
If[Mode==3, c=cfcomp;cw=cFit;wo=woFit;mAT=mATFit;{a2,a4}=OriTenFrArr[ns];Hillmean=MeanHillBasis[a2,a4]];
(* computation of mean of A coefs over lengths *)
{wmin,wmax}={Min[ws],Max[ws]};
Lf =  FromIsoCoef2HillCoef[3\[Kappa]f,2\[Mu]f];
Lm = FromIsoCoef2HillCoef[3\[Kappa]m,2\[Mu]m];
If[IdScheme==1,Af=AorW[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,wo,cw,wmin,wmax,Nw]];
If[IdScheme==2,Af=ArMTW[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,c,wo,cw,wmin,wmax,Nw]];
If[IdScheme==4,Af=ArLielW[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,c,wo,cw,wmin,wmax,Nw]];
If[IdScheme==6,Af=ArLielModW[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,c,wo,cw,wmin,wmax,Nw]];
If[IdScheme==5,Af=ArIDDW[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,c,formD,wo,cw,wmin,wmax,Nw]];
If[IdScheme==7,Af=ArPCWW[\[Kappa]m,\[Mu]m,\[Kappa]f,\[Mu]f,c,formD,wo,cw,wmin,wmax,Nw]];
Leff =Lm + c Contr4thTenHillBasis[(Lf-Lm),Af];
LeffTen = Total[Leff*Hillmean];
LeffProj = ProjecTIsoe3[LeffTen];
{LeffProj,LeffTen}]

(*
19. Effective behavior thermoelastic or thermo-viscoelastic behavior (Leff,\[Alpha]eff) and a distribution of lengths and orientations
Format for entries:
MatrixProps = {Subscript[\[Kappa], m], Subscript[\[Mu], m], Subscript[\[Alpha], m]};
InclusionProps = {Subscript[\[Kappa], i], Subscript[\[Mu], i], Subscript[\[Alpha], i]};
*)
MeanFieldSchemes\[Alpha]ATW[IdScheme_, MatrixProps_, InclusionProps_,formD_:0, MicroPars_,Nw_, Mode_]:=
Block[{cfTh, cfcomp, Ninc, ws, ns, wmArr, cWTh, woTh, mATTh, wmFit, cFit, woFit,
mATFit,cw,wo,mAT,a2,a4,Hillmean,wmin,wmax,\[Kappa]m,\[Mu]m,\[Alpha]m,\[Kappa]f,\[Mu]f,\[Alpha]f,Lf,Sf,Lm,Sm,\[Alpha]mean,Smean,PRossen,Af,Leff,LeffTen,\[Alpha]eff,\[Alpha]effvec,OP},
{\[Kappa]m,\[Mu]m,\[Alpha]m}=MatrixProps;
{\[Kappa]f,\[Mu]f,\[Alpha]f}=InclusionProps;
{cfTh, cfcomp, Ninc, ws, ns, wmArr, cWTh, woTh, mATTh, wmFit, cFit, woFit, mATFit}=MicroPars;
(* Mode setting *)
If[Mode==1, c=cfTh;cw=cWTh;wo=woTh;mAT=mATTh;{a2,a4}=OriTenATAxSym[mAT];Hillmean=MeanHillBasis[a2,a4]];
If[Mode==2, c=cfcomp;cw=cFit;wo=woFit;mAT=mATFit;{a2,a4}=OriTenATAxSym[mAT];Hillmean=MeanHillBasis[a2,a4]];
If[Mode==3, c=cfcomp;cw=cFit;wo=woFit;mAT=mATFit;{a2,a4}=OriTenFrArr[ns];Hillmean=MeanHillBasis[a2,a4]];
(* computation of mean of A coefs over lengths *)
{wmin,wmax}={Min[ws],Max[ws]};
Lf =  FromIsoCoef2HillCoef[3\[Kappa]f,2\[Mu]f];
Lm = FromIsoCoef2HillCoef[3\[Kappa]m,2\[Mu]m];
Sf =  FromIsoCoef2HillCoef[(3\[Kappa]f)^-1,(2\[Mu]f)^-1];
Sm = FromIsoCoef2HillCoef[(3\[Kappa]m)^-1,(2\[Mu]m)^-1];
\[Alpha]f = InclusionProps[[3]]*IdentityMatrix[3];
\[Alpha]m = MatrixProps[[3]]*IdentityMatrix[3];
\[Alpha]mean =  (1-c)\[Alpha]m + c \[Alpha]f;
PRossen= FromIsoCoef2HillCoef[3(1/\[Kappa]m-1/\[Kappa]f )^-1 , 2(1/\[Mu]m-1/\[Mu]f )^-1];
Smean =(1-c)Sm  + c Sf;
{Leff,LeffTen} = MeanFieldSchemesATW[IdScheme, {\[Kappa]m,\[Mu]m}, {\[Kappa]f,\[Mu]f},formD, MicroPars,Nw, Mode];
\[Alpha]eff = \[Alpha]mean + Contr4thAND2ndTen[TIso4thTensor[Contr4thTenHillBasis[(Inv4thHillCoef[Leff]-Smean),PRossen],3] ,(\[Alpha]m-\[Alpha]f)];
\[Alpha]effvec = {\[Alpha]eff[[1,1]],\[Alpha]eff[[3,3]]};
OP = {Leff,\[Alpha]effvec};
OP]

(* --------------------------------------------------------------------------------------------*)
(*
Function to compute effective thermoelastic and thermo-viscoeastic behavior by means of the correspondance principle using frequency data
*)

(*
20. Effective behavior in the frequency domain 
*)
MeanFieldSchemes\[Alpha]Freqs[IdScheme_, MatrixProps_, InclusionProps_,form_,c_,formD_:0,var_Symbol,
expminMec_,expmaxMec_,expminTh_,expmaxTh_, di_, nPartsMec_, nPartsTh_, Tolfit_:0.5]:= 
Block[{EffHill\[Infinity], EffHillF, \[Alpha]TildeF, Freqs, posMec, posTh, MatrixPropsF,
dataMec, FitMec, dataTh, FitTh, FitProps, OP,Mg,Mi,\[Tau]Mi,ErrsFit}, 
(* Computation of the long-term effective response *)
EffHill\[Infinity] = Flatten@MeanFieldSchemes\[Alpha][IdScheme, Limit[MatrixProps,var-> 0], InclusionProps,form,c,formD];
(* Initialization and setting up for frequency data treatment *)
Freqs = 10^Range[Min[expminMec,expminTh],Max[expmaxMec,expmaxTh],di];
posMec = Flatten[Position[Range[Min[expminMec,expminTh],Max[expmaxMec,expmaxTh],di],_?(#==expminMec||#==expmaxMec&)]];
posTh = Flatten[Position[Range[Min[expminMec,expminTh],Max[expmaxMec,expmaxTh],di],_?(#==expminTh||#==expmaxTh&)]];
(* Generation of frequency data *)
MatrixPropsF[s_] := MatrixProps/.{var -> s}; 
EffHillF = Table[0,{i,Length[Freqs]}];
Do[EffHillF[[i]] = Flatten@MeanFieldSchemes\[Alpha][IdScheme, MatrixPropsF[2 \[Pi] I Freqs[[i]]], InclusionProps,form,c,formD],
{i,Length[Freqs]}];
Do[EffHillF[[i]] = EffHillF[[i]]- EffHill\[Infinity],{i,1,Length[EffHillF]}];
dataMec = Table[0.,{i,6}];
Do[dataMec[[i]] =  Join[Transpose[{2 \[Pi] Freqs[[posMec[[1]];;posMec[[2]]]]}],Transpose[{Re[EffHillF[[posMec[[1]];;posMec[[2]],i]]]}],
Transpose[{Im[EffHillF[[posMec[[1]];;posMec[[2]],i]]]}],2],
{i,6}];
FitMec = Table[0.,{i,6}];
Do[FitMec[[i]] = GMMcouples[dataMec[[i]], nPartsMec, Tolfit],{i,6}];
dataTh = Table[0,{i,2}];
Do[dataTh[[i]] =  Join[Transpose[{2 \[Pi] Freqs[[posTh[[1]];;posTh[[2]]]]}],Transpose[{Re[EffHillF[[posTh[[1]];;posTh[[2]],6+i]]]}],
Transpose[{Im[EffHillF[[posTh[[1]];;posTh[[2]],6+i]]]}],2],
{i,2}];
FitTh = Table[0,{i,2}];
Do[FitTh[[i]] = DVEcouples[dataTh[[i]], nPartsTh, Tolfit],{i,2}];
FitProps = Join[FitMec,FitTh];
Mg = Table[0,{i,8}]; Mi = Mg; \[Tau]Mi = Mg;ErrsFit = Mg;
Do[
Mg[[i]] = EffHill\[Infinity][[i]] + Total[FitProps[[i,1]]];
Mi[[i]] = FitProps[[i,1]];
\[Tau]Mi[[i]]= FitProps[[i,2]];
ErrsFit[[i]]= FitProps[[i,3]];
,{i,8}];
OP = {Re@Mg,Mi,\[Tau]Mi,ErrsFit};
OP]

(*
21. Effective behavior in the frequency domain for a given distribution of lengths and orientations
*)
MeanFieldSchemes\[Alpha]FreqsATW[IdScheme_, MatrixProps_, InclusionProps_,formD_:0, MicroPars_,Nw_,Mode_,var_Symbol,
expminMec_,expmaxMec_,expminTh_,expmaxTh_, di_, nPartsMec_, nPartsTh_, Tolfit_:0.5]:= 
Block[{EffHill\[Infinity], EffHillF, \[Alpha]TildeF, Freqs, posMec, posTh, MatrixPropsF,
dataMec, FitMec, dataTh, FitTh, FitProps, OP,Mg,Mi,\[Tau]Mi,ErrsFit}, 
(* Computation of the long-term effective response *)
EffHill\[Infinity] = Flatten@MeanFieldSchemes\[Alpha]ATW[IdScheme,Limit[MatrixProps,var-> 0],InclusionProps,formD,MicroPars,Nw,Mode];
(* Initialization and setting up for frequency data treatment *)
Freqs = 10^Range[Min[expminMec,expminTh],Max[expmaxMec,expmaxTh],di];
posMec = Flatten[Position[Range[Min[expminMec,expminTh],Max[expmaxMec,expmaxTh],di],_?(#==expminMec||#==expmaxMec&)]];
posTh = Flatten[Position[Range[Min[expminMec,expminTh],Max[expmaxMec,expmaxTh],di],_?(#==expminTh||#==expmaxTh&)]];
(* Generation of frequency data *)
MatrixPropsF[s_] := MatrixProps/.{var -> s}; 
EffHillF = Table[0,{i,Length[Freqs]}];
Do[EffHillF[[i]] = Flatten@MeanFieldSchemes\[Alpha]ATW[IdScheme, MatrixPropsF[2 \[Pi] I Freqs[[i]]],InclusionProps,formD,MicroPars,Nw,Mode],
{i,Length[Freqs]}];
Do[EffHillF[[i]] = EffHillF[[i]]- EffHill\[Infinity],{i,1,Length[EffHillF]}];
dataMec = Table[0.,{i,6}];
Do[dataMec[[i]] =  Join[Transpose[{2 \[Pi] Freqs[[posMec[[1]];;posMec[[2]]]]}],Transpose[{Re[EffHillF[[posMec[[1]];;posMec[[2]],i]]]}],
Transpose[{Im[EffHillF[[posMec[[1]];;posMec[[2]],i]]]}],2],
{i,6}];
FitMec = Table[0.,{i,6}];
Do[FitMec[[i]] = GMMcouples[dataMec[[i]], nPartsMec, Tolfit],{i,6}];
dataTh = Table[0,{i,2}];
Do[dataTh[[i]] =  Join[Transpose[{2 \[Pi] Freqs[[posTh[[1]];;posTh[[2]]]]}],Transpose[{Re[EffHillF[[posTh[[1]];;posTh[[2]],6+i]]]}],
Transpose[{Im[EffHillF[[posTh[[1]];;posTh[[2]],6+i]]]}],2],
{i,2}];
FitTh = Table[0,{i,2}];
Do[FitTh[[i]] = DVEcouples[dataTh[[i]], nPartsTh, Tolfit],{i,2}];
FitProps = Join[FitMec,FitTh];
Mg = Table[0,{i,8}]; Mi = Mg; \[Tau]Mi = Mg;ErrsFit = Mg;
Do[
Mg[[i]] = EffHill\[Infinity][[i]] + Total[FitProps[[i,1]]];
Mi[[i]] = FitProps[[i,1]];
\[Tau]Mi[[i]]= FitProps[[i,2]];
ErrsFit[[i]]= FitProps[[i,3]];
,{i,8}];
OP = {Re@Mg,Mi,\[Tau]Mi,ErrsFit};
OP]


(* Functions to generate homogenized Prony series *)
HillfunMods[ListOfFunc_,var_,RelaxQ_:True]:= Block[{PronyMods,m\[Infinity],mg,mi, \[Tau]i,errors},
PronyMods = {};
If[RelaxQ==True,
Do[PronyMods= Append[PronyMods,FitPpolyGMM[ListOfFunc[[i]],var,{-17,17},7,0.05]],{i, Length[ListOfFunc]}]];
If[RelaxQ==False,
Do[PronyMods= Append[PronyMods,FitPpolyKVM[ListOfFunc[[i]],var,{-17,17},7,0.05]],{i, Length[ListOfFunc]}]];
m\[Infinity] = {};Do[m\[Infinity] = Append[m\[Infinity],PronyMods[[i]][[1]]],{i,Length[PronyMods]}];
mi = {};Do[mi = Append[mi,PronyMods[[i]][[2]]],{i,Length[PronyMods]}];
\[Tau]i = {};Do[\[Tau]i = Append[\[Tau]i,PronyMods[[i]][[3]]],{i,Length[PronyMods]}];
errors = {};Do[errors = Append[errors,PronyMods[[i]][[4]]],{i,Length[PronyMods]}];
mg={};
Do[mg = Append[mg,m\[Infinity][[i]]+Total[mi [[i]]]],{i,Length[PronyMods]}];
{mg,mi,\[Tau]i,errors}]



End[]
EndPackage[]







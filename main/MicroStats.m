(* ::Package:: *)

(* ::Section:: *)
(*Package for Microstructures Analysis *)
(*(Specific for data in Craft format )*)


(* ::Text:: *)
(*This package gives the statistical and deterministic analysis of a microstructure to be compared with Craft full-field homogenization approach *)


(* ::Input:: *)
(*(* Name: MicroStats *)*)
(*(* Author: Camilo Andr\[EAcute]s SUAREZ AFANADOR *)*)
(*(* Date : Aug 23,2021 *)*)
(* *)


BeginPackage["MicroStats`"]

(* Declaring output functions for context *)
OriTenATAxSym::usage = 
"OriTenATAxSym[n] gives the mean orientation tensors {a2,a4} for a given Advani & Tucker's Parameter n of an 
axisymetric orientation distribution around the \!\(\*OverscriptBox[\(z\), \(\[RightVector]\)]\) axis "
OriTenFrArr::usage = 
"OriTenFrArr[nList] gives the mean orientation tensors {a2,a4} for a given list of orientations"
MeanHillBasis::usage = 
"MeanHillBasis[a2,a4] gives the 6 Hill's basis tensors for a given mean concentration tensors a2, a4 "
ReadMicroATW::usage = 
"ReadMicroATW[micropath] Reads a microstructure from Craft generator in extension file.txt"
ReadMicroATWV2::usage = 
"ReadMicroATW[micropath] Reads a microstructure from Craft generator in extension file.txt"
WeibullFit::usage = 
"WeibullFit[ws,cinit=1.0,tol=1.0*^-12] gives the two parameters of the weibull law (\"{wo,c}\") from a list of aspect ratios"
AngleFra2::usage = 
" AngleFra2[a2] gives the three angles of a vector for a given second order 
orientation tensor a2 "
WeibullArray::usage = 
"WeibullArray[wo,c,wmin,wmax,Nw] gives the outputs {wms,fws,wtot,weigs} for a given weibull parameters and discretization choice 
Nw "
WeibullFun::usage = 
"WeibullFun[w,c,wo] gives the value of the PDF at a given aspect w and weibull parameters wo and c"

Begin["`Private`"]

(* Importing subpackage *)
SetDirectory[NotebookDirectory[]];
Needs["TensorCalc`"]

(* 1. Leght distribution functions *)
WeibullFun[w_,c_,wo_]:= PDF[WeibullDistribution[c,wo],w];

WeibullMeanFrwondc[wo_,c_]:= wo Gamma[1+1/c];
WeibullwoFrwmndc[wm_,c_]:=wm Gamma[1+1/c]^-1;

fc[c_,lnw_,ws_]:= 1/c + Mean[lnw] - Mean[lnw*ws^c]/Mean[ws^c];
dfc[c_,lnw_,ws_]:= -(1/c^2) - ((Mean[lnw^2*ws^c]*Mean[ws^c])-Mean[lnw*ws^c]^2)/Mean[ws^c]^2; 
WeibullFit[ws_,cinit_:1.0,tol_:1.0*^-12]:= Block[{c,lnw,wo}, c=cinit; lnw = Log[ws]; 
While[Abs[fc[c,lnw,ws]] > tol, c = c - fc[c,lnw,ws]/dfc[c,lnw,ws] ];
wo = Mean[ws*c]^(1/c);
{wo,c}]

WeibullArray[wo_,cw_,wmin_,wmax_,Nw_]:= Block[{dw, ws, wms,wtot,fws,weigs},
dw = (wmax-wmin)/Nw;ws= Range[wmin,wmax,dw];wms = (ws[[;;-2]]+ws[[2;;]])/2;  
fws = WeibullFun[wms,cw,wo]; wtot = Total[wms*fws]*dw ;
weigs = (wms/wtot)*fws*dw;
{wms,fws,wtot,weigs}]


(* 2. Orientation Tensors *)
(* Probabilistic Orientation 2nd and 4th order tensors from Advani & Tucker's Law in Voigt representation 
for an axisymmetric distribution of fibers around the Overscript[z, \[RightVector]] axis *)
OriTenATAxSym[n_]:= Block[{a2, a4, a2Rot, a4Rot}, a2 = (n+3)^-1 DiagonalMatrix[{n+1,1,1}]; 
a4 = ((n+3)(n+5))^-1 {{(n+1)(n+3),(n+1),(n+1),0,0,0},{(n+1),3,1,0,0,0},{(n+1),1,3,0,0,0},{0,0,0,1,0,0},
{0,0,0,0,(n+1),0},{0,0,0,0,0,(n+1)}};
a4Rot = Rot4thTenMatNota[a4,0,\[Pi]/2,0]; a2Rot = Rot2ndTenAxis[a2,\[Pi]/2,2];
{a2Rot, a4Rot}]

(* Orientation tensors from an array  *)
OriTenFrArr[nList_]:= Block[{N,a2,a4},N = Length[nList]; 
a2 = (1/N)Sum[Table[nList[[r,i]]nList[[r,j]],{i,3},{j,3}],{r,N}];
a4 = (1/N)Sum[Table[nList[[r,i]]nList[[r,j]]nList[[r,k]]nList[[r,l]],{i,3},{j,3},{k,3},{l,3}],{r,N}];
{a2,a4}]

(* Angles \[Theta]x,\[Theta]y,\[Theta]z from a2 *)
AngleFra2[a2_] := Block[{nvec,e1,e2,e3,\[Theta]x,\[Theta]y,\[Theta]z},nvec = Sqrt[{a2[[1,1]],a2[[2,2]],a2[[3,3]]}];
e1 = {1,0,0}; e2 = {0,1,0}; e3 = {0,0,1};
\[Theta]x = ArcCos[nvec . e1]; \[Theta]y = ArcCos[nvec . e2]; \[Theta]z = ArcCos[nvec . e3];
{\[Theta]x,\[Theta]y,\[Theta]z}]


(* Mean Hill Basis  given a2 a4 Tensors *)
Kd[i_,j_]:= KroneckerDelta[i,j];
I4Ten = (1/2)Table[Kd[i,k]Kd[j,l]+Kd[i,l]Kd[j,k],{i,3},{j,3},{k,3},{l,3}];
MeanHillBasis[a2_,a4_]:= Block[{TArr,T,MatFrTen,List1,List2,Rl,a4Ten,Id,\[Alpha],H1m,H2m,H3m,H4m,H5m,H6m},
If[Dimensions[a4]=={6,6},
TArr = Array[T,{3,3,3,3}];
MatFrTen = {{TArr[[1,1,1,1]],TArr[[1,1,2,2]],TArr[[1,1,3,3]],TArr[[1,1,2,3]],
TArr[[1,1,3,1]],TArr[[1,1,1,2]]},{TArr[[2,2,1,1]],TArr[[2,2,2,2]],TArr[[2,2,3,3]],TArr[[2,2,2,3]],TArr[[2,2,1,3]],TArr[[2,2,1,2]]},
{TArr[[3,3,1,1]],TArr[[3,3,2,2]],TArr[[3,3,3,3]],TArr[[3,3,2,3]],TArr[[3,3,1,3]],TArr[[3,3,1,2]]},{TArr[[2,3,1,1]],TArr[[2,3,2,2]],
TArr[[2,3,3,3]],TArr[[2,3,2,3]],TArr[[2,3,1,3]],TArr[[2,3,1,2]]},{TArr[[1,3,1,1]],TArr[[1,3,2,2]],TArr[[1,3,3,3]],TArr[[1,3,2,3]],
TArr[[1,3,1,3]],TArr[[1,3,1,2]]},{TArr[[1,2,1,1]],TArr[[1,2,2,2]],TArr[[1,2,3,3]],TArr[[1,2,2,3]],TArr[[1,2,1,3]],TArr[[1,2,1,2]]}};
List1 = Flatten[MatFrTen];List2 = Flatten[a4]; Rl = 0*List1; Do[Rl[[i]] = List1[[i]]-> List2[[i]],{i,Length[List1]}];
TArr = TArr/.Rl; Do[TArr[[i,j,l,k]] = TArr[[j,i,k,l]] =TArr[[i,j,k,l]],{i,3},{j,3},{k,3},{l,3}];T[1,1,1,3] = 0;
a4Ten=TArr,a4Ten=a4];
Id = IdentityMatrix[3];
\[Alpha] = Table[Kd[i,k]a2[[l,j]]+ Kd[i,l]a2[[k,j]] + a2[[i,k]]Kd[l,j]+ a2[[i,l]]Kd[k,j],{i,3},{j,3},{k,3},{l,3}];
H1m = (1/2)(Table[Id[[i,j]]Id[[k,l]] - Id[[i,j]]a2[[k,l]] - a2[[i,j]]Id[[k,l]],{i,3},{j,3},{k,3},{l,3}]+a4Ten);
H2m = Table[Id[[i,j]]a2[[k,l]],{i,3},{j,3},{k,3},{l,3}]-a4Ten;
H3m = Table[a2[[i,j]]Id[[k,l]],{i,3},{j,3},{k,3},{l,3}]-a4Ten;
H4m = a4Ten;
H5m = (1/2)(2*I4Ten - \[Alpha] + a4Ten - Table[Id[[i,j]]Id[[k,l]] - Id[[i,j]]a2[[k,l]] - a2[[i,j]]Id[[k,l]],{i,3},{j,3},{k,3},{l,3}]);
H6m = (1/2)\[Alpha] - 2*a4Ten;
{H1m,H2m,H3m,H4m,H5m,H6m}]

(* Reading Microstructure files *)
ReadMicroATW[micropath_]:=Block[{microinfo, microraw, Ninc, As, Bs, Xs, Rs, Ls,ns,cfcomp,Vol, Vf, ws, woWTh,cWTh,wm,mATTh,
 woFit,cFit,wmFit,mmin,mmax,ms,a2m,a4m,a2Arr,a4Arr,a2Norm,a4Norm,erra2, erra4, mATFit,cfTh},
microinfo = Import[micropath,"Table"]; microraw = microinfo[[14;;]]; Ninc = Dimensions[microraw][[1]]; 
As = microraw[[All,2;;4]]; Bs = microraw[[All,5;;7]]; Xs = 0.5(As+Bs); Rs = microraw[[All,8]];
Ls = ns = Table[0,{i,Length[As]}]; Do[Ls[[i]] = Norm[Bs[[i]]-As[[i]]]; ns[[i]] = (Bs[[i]]-As[[i]])/Ls[[i]],{i,Length[As]}]; 
Vf = \[Pi] Total[Ls (Rs)^2 + 4/3 Rs^3]; ws = Ls/(2 Rs);wm = Mean[ws];
Vol = microinfo[[8,4]]microinfo[[8,5]]*Read[StringToStream[StringSplit[StringSplit[microinfo[[8,3]],"="][[2]]][[1]]],Number];
cfcomp = Vf/Vol; woWTh = Read[StringToStream[StringSplit[microinfo[[4,4]],"="][[2]]]];
cWTh = Read[StringToStream[StringSplit[microinfo[[11,9]],"="][[2]]]];
mATTh = Read[StringToStream[StringSplit[microinfo[[10,4]],_?LetterQ][[18]]]];
cfTh = Read[StringToStream[Flatten[StringSplit[StringSplit[microinfo[[5,3]],_?LetterQ],"="]][[1]]]];
{woFit,cFit} = WeibullFit[ws]; wmFit = WeibullMeanFrwondc[woFit,cFit];
{mmin,mmax} = {mATTh-10, mATTh+20};ms = Table[i,{i,mmin,mmax+2,(mmax+2-mmin)/500}];
{a2Arr,a4Arr} = OriTenFrArr[ns];a4Arr = FrTenToMat[a4Arr];
{a2Norm, a4Norm} = {Norm[a2Arr], Norm[a4Arr]}; 
erra2 = erra4 = 0*ms; Do[{a2m,a4m} = OriTenATAxSym[ms[[i]]];erra2[[i]] = Norm[a2m-a2Arr]/a2Norm;
erra4[[i]] = Norm[a4m-a4Arr]/a4Norm,{i,Length[erra2]}];
mATFit = ms[[Flatten[Position[erra4,_?(#==Min[erra4]&)]][[1]]]];
{Vol, cfTh, cfcomp, Ninc, ws, Xs, Rs, Ls, ns, wm, cWTh, woWTh, mATTh, wmFit, cFit, woFit,mATFit}]

ReadMicroATWV2[micropath_,mATTh_,cWTh_,woWTh_,cfTh_]:=Block[{microinfo, microraw, Ninc, As, Bs, Xs, Rs, Ls,ns,cfcomp,Vol, Vf, ws,wm,
 woFit,cFit,wmFit,mmin,mmax,ms,a2m,a4m,a2Arr,a4Arr,a2Norm,a4Norm,erra2, erra4, mATFit},
microinfo = Import[micropath,"Table"]; microraw = microinfo[[14;;]]; Ninc = Dimensions[microraw][[1]]; 
As = microraw[[All,2;;4]]; Bs = microraw[[All,5;;7]]; Xs = 0.5(As+Bs); Rs = microraw[[All,8]];
Ls = ns = Table[0,{i,Length[As]}]; Do[Ls[[i]] = Norm[Bs[[i]]-As[[i]]]; ns[[i]] = (Bs[[i]]-As[[i]])/Ls[[i]],{i,Length[As]}]; 
Vf = \[Pi] Total[Ls (Rs)^2 + 4/3 Rs^3]; ws = Ls/(2 Rs);wm = Mean[ws];
Vol = microinfo[[8,4]]microinfo[[8,5]]*Read[StringToStream[StringSplit[StringSplit[microinfo[[8,3]],"="][[2]]][[1]]],Number];
cfcomp = Vf/Vol; {woFit,cFit} = WeibullFit[ws]; wmFit = WeibullMeanFrwondc[woFit,cFit];
{mmin,mmax} = {mATTh-10, mATTh+20};ms = Table[i,{i,mmin,mmax+2,(mmax+2-mmin)/500}];
{a2Arr,a4Arr} = OriTenFrArr[ns];a4Arr = FrTenToMat[a4Arr];
{a2Norm, a4Norm} = {Norm[a2Arr], Norm[a4Arr]}; 
erra2 = erra4 = 0*ms; Do[{a2m,a4m} = OriTenATAxSym[ms[[i]]];erra2[[i]] = Norm[a2m-a2Arr]/a2Norm;
erra4[[i]] = Norm[a4m-a4Arr]/a4Norm,{i,Length[erra2]}];
mATFit = ms[[Flatten[Position[erra4,_?(#==Min[erra4]&)]][[1]]]];
{cfTh, cfcomp, Ninc, ws, ns, wm, cWTh, woWTh, mATTh, wmFit, cFit, woFit,mATFit//N}]


End[]
EndPackage[]

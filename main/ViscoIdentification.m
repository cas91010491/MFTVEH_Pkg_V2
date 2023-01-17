(* ::Package:: *)

(* ::Section:: *)
(*Package for Identification of discrete viscoelastic spectra*)


(* ::Text:: *)
(*This first package carries with tensor calculus operations in  both harmonic and Hill's decomposition*)


(* ::Input:: *)
(*(* Name: ViscoIdentification *)*)
(*(* Author: Camilo Andr\[EAcute]s SUAREZ AFANADOR *)*)
(*(* Date : May 5,2021 *)*)
(* *)


BeginPackage["ViscoIdentification`"]

GMMcouples::usage = 
"GMMcouples[\"data\", \"n\", \"eps\"]
Gives a vector {\"Mods\",\"times\",\"Error\"} containing the discrete distribution of relaxation moduli and times, and the \!\(\*SubscriptBox[\(L\), \(2\)]\) error of identification, 
for a given DMA data matrix ({\"f\",\"RealMod\",\"ImagMod\"}), the number of data-partitions \"n\", and the tolerance \"eps\" for the KN-HW methodology,
this approach considers a preprocessed data obtained by substracting the lowest frequency \"RealMod\" from The data values "
GKVMcouples::usage = 
"GKVMcouples[\"data\", \"n\", \"eps\"]
Gives a vector {\"Mods\",\"times\",\"Error\"} containing the discrete distribution of Creep moduli and times, and the \!\(\*SubscriptBox[\(L\), \(2\)]\) error of identification, \
for a given DMA data matrix ({\"f\",\"RealMod\",\"ImagMod\"}), the number of data-partitions \"n\", and the tolerance \"eps\" for the KN-HW methodology, \
this approach considers a preprocessed data obtained by substracting the lowest frequency \"RealMod\" from The data values"
DVEcouples::usage = 
"DVEcouples[\"data\", \"n\", \"eps\"]
Gives a vector {\"Mods\",\"times\",\"Error\"} containing the discrete distribution of viscoelastic moduli and times with no constrains in the sign of
 the discrete coefficients, and the \!\(\*SubscriptBox[\(L\), \(2\)]\) error of identification, for a given DMA data matrix ({\"f\",\"RealMod\",\"ImagMod\"}),
 the number of data-partitions \"n\", and the tolerance \"eps\" for the KN-HW methodology,
this approach considers a preprocessed data obtained by substracting the lowest frequency \"RealMod\" from The data values"
FitPpolyGMM::usage = 
"FitPpolyGMM[\"Ppoly\", \"var\", \"explims\", \"n\", \"eps\"] 
Gives a vector {\"\!\(\*SubscriptBox[\(Mod\), \(\[Infinity]\)]\)\",\"Mods\",\"times\",\" Error \"} containing the long-term modulus, the discrete distribution of viscoelastic moduli and times,
 and the \!\(\*SubscriptBox[\(L\), \(2\)]\) error of identification, given: Rational polynomial \"Ppoly\" with argument \"var\", the logarithmic frequency interval {\"expmin\",\"expmax\"},
the number of data-partitions \"n\", and the tolerance \"eps\" for the KN-HW methodology in a Generalized Maxwell Model"
FitPpolyKVM::usage = 
"FitPpolyKVM[\"Ppoly\", \"var\", \"explims\", \"n\", \"eps\"] 
Gives a vector {\"\!\(\*SubscriptBox[\(Mod\), \(\[Infinity]\)]\)\",\"Mods\",\"times\",\"Error\"} containing the long-term modulus, the discrete distribution of viscoelastic moduli and times,
 and the \!\(\*SubscriptBox[\(L\), \(2\)]\) error of identification, given: Rational polynomial \"Ppoly\" with argument \"var\", the logarithmic frequency interval {\"expmin\",\"expmax\"},
the number of data-partitions \"n\", and the tolerance \"eps\" for the KN-HW methodology in a Generalized Kelvin-Voigt Model"
FitPpolyDVE::usage = 
"FitPpolyDVE[\"Ppoly\", \"var\", \"explims\", \"n\", \"eps\"] 
Gives a vector {\"\!\(\*SubscriptBox[\(Mod\), \(\[Infinity]\)]\)\",\"Mods\",\"times\",\"Error\"} containing the long-term modulus, the discrete distribution of viscoelastic moduli and times,
 and the \!\(\*SubscriptBox[\(L\), \(2\)]\) error of identification, given: Rational polynomial \"Ppoly\" with argument \"var\", the logarithmic frequency interval {\"expmin\",\"expmax\"},
the number of data-partitions \"n\", and the tolerance \"eps\" for the KN-HW methodology in a Discrete viscoelastic model with no constrains in the sign
of the coefficients distribution"
HWfun::usage = 
"HWfun[\"ws\", \"ts\", \"RealMod\", \"ImagMod\"]
Gives a list of Prony weights given : \
\"ws\", the list of sample frequencies. \
\"ts\", the list of discrete time-distribution. \
\"RealMod\", the list of storage modulus. \
\"ImagMod\", the list of Loss modulus."
DVEcouplesColl::usage = 
"DVEcouplesColl[\"data\", \"n\", \"RelaxQ\"]
Gives the Porny parameters and the approximation error ({Mods,\[Tau]vals,L2Err}) given: \
\"data\", the data in format {\[Omega]_i,Storage_i, Loss_i}. \
\"n\", the number of time relaxations for the Prony series. \
\"RelaxQ\", the Boolean diferentiating creep or relaxation modulus"


Begin["`Private`"]

(* M1 & M2 matrices *)
M1M2MatsAndNull[\[Omega]vals_,ComplexM_]:= Block[{ConjM,M1,M2,NullM1,NullM2},
ConjM= Conjugate[ComplexM];M1 = Table[0,{i,Length[ComplexM]},{j,Length[ConjM]}];M2=M1; 
Do[M1[[i,j]]=(I(ConjM[[i]]-ComplexM[[j]]))/(\[Omega]vals[[i]]+\[Omega]vals[[j]]),{i,Length[ComplexM]},{j,Length[ConjM]}];
Do[M2[[i,j]]=(ConjM[[i]]/\[Omega]vals[[i]]+ComplexM[[j]]/\[Omega]vals[[j]])/(\[Omega]vals[[i]]+\[Omega]vals[[j]]),{i,Length[ComplexM]},{j,Length[ConjM]}];
  {NullM1,NullM2} = {NullSpace[M1],NullSpace[M2]};{NullM1,NullM2}]
(* s functions *)
Ffun[NullM1_,NullM2_,s_,x_]:= Block[{f1Mod,f2Mod},f1Mod = Norm[NullM1 . (s+I x)^-1];
f2Mod = Norm[NullM2 . (s+I x)^-1];{f1Mod,f2Mod}]

(* Inflection points search functions *)
SingleInflPoints[x_,f_]:=Block[{infl},infl = {};
Do[If[Sign[f[[i]]-f[[i-1]]]==-1 && Sign[f[[i+1]]-f[[i]]]==1,infl=Append[infl,x[[i]]]],{i,2,Length[x]-1}];
infl]
CommonInflPoints[x_,f1_,f2_,eps_:0.5]:=Block[{infl1,infl2,infl},infl={};
infl1=SingleInflPoints[x,f1];infl2=SingleInflPoints[x,f2]; 
Do[If[Abs[infl1[[i]]-infl2[[j]]]<eps,infl = Append[infl,(infl1[[i]]+infl2[[j]])*0.5]],{i,Length[infl1]},{j,Length[infl2]}];
infl]

(* HW method function *)
HWfun::emptyarray = "The discret times list is empty";
HWfun[ws_,ts_List /; Length[ts]>0 ,RealMod_,ImagMod_]:=Block[{aji,bji,A,B,ModTot,Atot,Mods},
aji= Table[0,{i,Length[ws]},{j,Length[ts]}];
bji= Table[0,{i,Length[ws]},{j,Length[ts]}];
Do[aji[[j,i]]=(ws[[j]]ts[[i]])^2/(1+(ws[[j]]ts[[i]])^2),{j,Length[ws]},{i,Length[ts]}];
Do[bji[[j,i]]=(ws[[j]]ts[[i]])/(1+(ws[[j]]ts[[i]])^2),{j,Length[ws]},{i,Length[ts]}];
A = Transpose[aji];B=Transpose[bji]; ModTot= Join[RealMod,ImagMod,1];Atot = Join[A,B,2];
Mods = PseudoInverse[Atot . Transpose[Atot],Tolerance->10^(-15)] . (Atot . ModTot);
Mods]
HWfun[ws_,ts_List ,RealMod_,ImagMod_] := (Message[HWfun::emptyarray];$Failed)

(* Complex modulus computations function *)
ComplexMod[\[Omega]s_,Mi_,\[Tau]i_]:=Block[{fun},fun = Sum[Mi[[i]]((I \[Omega]s)/(I \[Omega]s  +(\[Tau]i[[i]]^-1) )),{i,Length[Mi]}];fun]

(* Function to visualize inflection points *)
PlotSfun[sl_,f1l_,f2l_]:=Block[{plot},
plot= ListLinePlot[{Thread[{sl^-1,f1l}],Thread[{sl^-1,f2l}]},
PlotStyle->{Black,Red},GridLines->Automatic,Frame->True,Axes-> False, PlotLegends->{"\!\(\*SubscriptBox[\(f\), \(1\)]\)","\!\(\*SubscriptBox[\(f\), \(2\)]\)"},
ScalingFunctions->{"Log10","Log10"},PlotLabel->"f functions vs \!\(\*SuperscriptBox[\(s\), \(-1\)]\)"]
;plot]

(* Approximation Error function *)
L2ComplexErrNorm[Testdata_,Refdata_] :=Block[{Error},Error = Sqrt[(Refdata- Testdata) . Conjugate[Refdata- Testdata]/Refdata . Conjugate[Refdata]];Error] 

(* Function to find time spectrum and associated modulus *)
(* n\[Equal]0 *)
doIfnEq0[\[Omega]vals_,ComplM_,expmin_,expmax_,eps_:0.5]:=Block[{svals,f1Modvals,f2Modvals,NullM1,NullM2,log10f1,log10f2,loginvs,plotfs,\[Tau]s},
svals= Table[10^i,{i,expmin,expmax,0.01}];f1Modvals= 0.0*svals; f2Modvals =f1Modvals;{NullM1,NullM2}=M1M2MatsAndNull[\[Omega]vals,ComplM];
Do[{f1Modvals[[i]],f2Modvals[[i]]}=Ffun[NullM1,NullM2,svals[[i]],\[Omega]vals],{i,Length[svals]}];
{log10f1,log10f2,loginvs} ={Log10[f1Modvals],Log10[f2Modvals],Log10[svals^-1]};
plotfs = PlotSfun[svals,f1Modvals,f2Modvals];\[Tau]s = 10^CommonInflPoints[loginvs,log10f1,log10f2,eps];
{plotfs,\[Tau]s}]
(* n\[NotEqual]0 *)
doIfnNeq0[\[Omega]vals_,ComplM_,expl_,expmin_,expmax_,n_,eps_:0.5]:=Block[
{\[CapitalDelta]exp,InterLen,\[Tau]s,ncount,expminInt,expmaxInt,idxs,\[Omega]in,Complin,NullM1in,NullM2in,exps,svalsin,f1Modvalsin,f2Modvalsin,log10f1in,log10f2in,loginvsin,plotfsin},
\[CapitalDelta]exp=expmax-expmin;InterLen = \[CapitalDelta]exp/n; \[Tau]s = {};plotfsin={};ncount=1;
expminInt = expmin;
While[ncount<= n,expmaxInt = expminInt +InterLen; idxs = Flatten[Position[expl,_?(expminInt<#<expmaxInt &)]];
If[ncount==1,idxs= Append[idxs,Last[idxs]+1]];If[ncount==n,idxs= Prepend[idxs,First[idxs]-1]];
If[ncount!= n && ncount!= 1,idxs= Prepend[idxs,First[idxs]-1];idxs= Append[idxs,Last[idxs]+1]];
{\[Omega]in,Complin}={idxs*0.0,idxs*0.0};
Do[{\[Omega]in[[i]],Complin[[i]]}={\[Omega]vals[[idxs[[i]]]],ComplM[[idxs[[i]]]]},{i,Length[idxs]}];
{NullM1in,NullM2in}=M1M2MatsAndNull[\[Omega]in,Complin];
If[ncount==1,exps=Table[i,{i,expminInt-5,expmaxInt+0.5,0.001}]];If[ncount==n,exps=Table[i,{i,expminInt-0.5,expmaxInt+5,0.001}]];
If[ncount!= n && ncount!= 1,exps=Table[i,{i,expminInt-0.5,expmaxInt+0.5,0.001}]];svalsin = 10^exps;
f1Modvalsin= 0.0*svalsin; f2Modvalsin =f1Modvalsin;
Do[{f1Modvalsin[[i]],f2Modvalsin[[i]]}=Ffun[NullM1in,NullM2in,svalsin[[i]],\[Omega]in],{i,Length[svalsin]}];
{log10f1in,log10f2in,loginvsin} ={Log10[f1Modvalsin],Log10[f2Modvalsin],Log10[svalsin^-1]};
plotfsin =Append[plotfsin, PlotSfun[svalsin,f1Modvalsin,f2Modvalsin]];
\[Tau]s = Append[\[Tau]s,10^CommonInflPoints[loginvsin,log10f1in,log10f2in,eps]];expminInt=expmaxInt;
ncount++];\[Tau]s=DeleteDuplicates[Sort[Flatten[\[Tau]s]],Abs[Log10[#1]-Log10[#2]]<10^-7&]
;{plotfsin,\[Tau]s}]
(* body *)
GMMcouples[data_,n_,eps_:0.5]:=Block[
{\[Omega]vals,RealM,ImagM,expl, expmin, expmax,ComplexM,plotf,\[Tau]vals,count,Mods,Negidx,testdata,L2Err},
\[Omega]vals= data[[All,1]]; RealM= data[[All,2]]; ImagM= data[[All,3]];expl= Log10[\[Omega]vals]; expmin= Min[expl];
expmax= Max[expl]; ComplexM = RealM + I ImagM;If[n==0,{plotf, \[Tau]vals}=doIfnEq0[\[Omega]vals,ComplexM,expmin,expmax,eps]];
If[n!= 0,{plotf, \[Tau]vals}=doIfnNeq0[\[Omega]vals,ComplexM,expl,expmin,expmax,n,eps]];
count=1;
While[count==1, 
Mods=HWfun[\[Omega]vals,\[Tau]vals,RealM,ImagM];
If[AnyTrue[Mods,Negative]==True,Negidx = Position[Mods,_?(#<0&)];\[Tau]vals=Delete[\[Tau]vals,Negidx];count=1];
If[AnyTrue[Mods,Negative]==False,count=-1]];
testdata=ComplexMod[\[Omega]vals,Mods,\[Tau]vals];
L2Err =Re[ L2ComplexErrNorm[testdata,ComplexM]];
{Mods,\[Tau]vals,L2Err}]
(* body for compliance mods *)
GKVMcouples::emptyarr = "No discret times left from optimization handler";
GKVMcouples[data_,n_,eps_:0.5]:=Block[
{\[Omega]vals,RealM,ImagM,expl, expmin, expmax,ComplexM,plotf,\[Tau]vals,\[Tau]valso,count,Mods,Modso,Posidx,testdata,L2Err},
\[Omega]vals= data[[All,1]]; RealM= data[[All,2]]; ImagM= data[[All,3]];
expl= Log10[\[Omega]vals]; expmin= Min[expl];expmax= Max[expl];
ComplexM = RealM + I ImagM;
If[n==0,{plotf, \[Tau]vals}=doIfnEq0[\[Omega]vals,ComplexM,expmin,expmax,eps]];
If[n!= 0,{plotf, \[Tau]vals}=doIfnNeq0[\[Omega]vals,ComplexM,expl,expmin,expmax,n,eps]];
(*count=1;*)
Mods=HWfun[\[Omega]vals,\[Tau]vals,RealM,ImagM];
Modso = Mods;
\[Tau]valso = \[Tau]vals;
While[AnyTrue[Mods,Positive]==True,
Posidx = Position[Mods,_?(#>0&)];
\[Tau]vals=Delete[\[Tau]vals,Posidx];
If[Length[\[Tau]vals]<1,Print[Column[{Style["First computed discret values and error",FontColor->Red],Modso,\[Tau]valso,Re[L2ComplexErrNorm[ComplexMod[\[Omega]vals,Modso,\[Tau]valso],ComplexM]]}]];
Print[Style[GKVMcouples::emptyarr,FontColor->Red]];Abort[]];
Mods=HWfun[\[Omega]vals,\[Tau]vals,RealM,ImagM];
];
testdata=ComplexMod[\[Omega]vals,Mods,\[Tau]vals];
L2Err =Re[ L2ComplexErrNorm[testdata,ComplexM]];
{Mods,\[Tau]vals,L2Err}]
(* body for compliance mods no restriction *)
DVEcouples[data_,n_,eps_:0.5]:=Block[
{\[Omega]vals,RealM,ImagM,expl, expmin, expmax,ComplexM,plotf,\[Tau]vals,count,Mods,Posidx,testdata,L2Err},
\[Omega]vals= data[[All,1]]; RealM= data[[All,2]]; ImagM= data[[All,3]];expl= Log10[\[Omega]vals]; expmin= Min[expl];expmax= Max[expl];
ComplexM = RealM + I ImagM;If[n==0,{plotf, \[Tau]vals}=doIfnEq0[\[Omega]vals,ComplexM,expmin,expmax,eps]];
If[n!= 0,{plotf, \[Tau]vals}=doIfnNeq0[\[Omega]vals,ComplexM,expl,expmin,expmax,n,eps]];
Mods=HWfun[\[Omega]vals,\[Tau]vals,RealM,ImagM];
testdata=ComplexMod[\[Omega]vals,Mods,\[Tau]vals];
L2Err =Re[ L2ComplexErrNorm[testdata,ComplexM]];
{Mods,\[Tau]vals,L2Err}]
(* body for relaxation modulus approximation from collocation method *)
DVEcouplesColl::notimesleft = "No discret times left from optimization handler";
DVEcouplesColl[data_,n_,RelaxQ_]:=Block[
{\[Omega]vals,RealM,ImagM,expl, expmin, expmax,ComplexM,plotf,\[Tau]vals,count,Mods,idx,testdata,L2Err,Bool1,Modso,\[Tau]valso},
\[Omega]vals= data[[All,1]]; RealM= data[[All,2]]; ImagM= data[[All,3]];expl= Log10[\[Omega]vals]; expmin= Min[expl];
expmax= Max[expl]; ComplexM = RealM + I ImagM;
\[Tau]vals=Table[10^i,{i,expmin,expmax,(expmax-expmin)/n}];
Mods=HWfun[\[Omega]vals,\[Tau]vals,RealM,ImagM];
Modso=Mods;\[Tau]valso=\[Tau]vals;
Bool1 = If[RelaxQ==True,Negative,Positive];
While[AnyTrue[Mods,Bool1]==True, 
If[RelaxQ==True, idx=Position[Mods,_?(#<0&)]];
If[RelaxQ==False, idx=Position[Mods,_?(#>0&)]];
\[Tau]vals=Delete[\[Tau]vals,idx];
If[Length[\[Tau]vals]<1,Print[Column[{Style["First computed discret values and error",FontColor->Red],Modso,\[Tau]valso,Re[L2ComplexErrNorm[ComplexMod[\[Omega]vals,Modso,\[Tau]valso],ComplexM]]}]];
Print[Style[DVEcouplesColl::notimesleft,FontColor->Red]];Abort[]];
Mods=HWfun[\[Omega]vals,\[Tau]vals,RealM,ImagM]
];
testdata=ComplexMod[\[Omega]vals,Mods,\[Tau]vals];
L2Err =Re[L2ComplexErrNorm[testdata,ComplexM]];
{Mods,\[Tau]vals,L2Err,Bool1}]


(* Function for a given p-polynomyal 
explims = {expmin,expmax}
*)
(* Function to Fit Prony-Series for a relaxation function (positivity constrain for discrete moduli distribution) *)
FitPpolyGMM[Ppoly_,var_Symbol,explims_,n_,eps_:0.5]:=Block[{m\[Infinity],data,\[Omega]sample,mi,\[Tau]i,errorL2},
m\[Infinity] = Limit[Ppoly,var-> 0]//N;\[Omega]sample=Table[10^i,{i,explims[[1]],explims[[2]],0.01}];
data=  Join[Transpose[{\[Omega]sample}],Transpose[{Re[ReplaceAll[var-> I*\[Omega]sample][Ppoly]]-m\[Infinity]}],Transpose[{Im[ReplaceAll[var-> I*\[Omega]sample][Ppoly]]}],2];
{mi,\[Tau]i,errorL2}=GMMcouples[data,n,eps];
{m\[Infinity],mi,\[Tau]i,errorL2}]
(* Function to Fit Prony-Series for a Creep function (negativity constrain for discrete moduli distribution) *)
FitPpolyKVM[Ppoly_,var_Symbol,explims_,n_,eps_:0.5]:=Block[{m\[Infinity],data,\[Omega]sample,mi,\[Tau]i,errorL2},
m\[Infinity] = Limit[Ppoly,var-> 0]//N;\[Omega]sample=Table[10^i,{i,explims[[1]],explims[[2]],0.01}];
data=  Join[Transpose[{\[Omega]sample}],Transpose[{Re[Ppoly/.{var-> I \[Omega]sample}]-m\[Infinity]}],Transpose[{Im[Ppoly/.{var-> I \[Omega]sample}]}],2];
{mi,\[Tau]i,errorL2}=GKVMcouples[data,n,eps];{m\[Infinity],mi,\[Tau]i,errorL2}]
(* Function to Fit Prony-Series for a viscoelastic function (No constrain in discrete moduli distribution) *)
FitPpolyDVE[Ppoly_,var_Symbol,explims_,n_,eps_:0.5]:=Block[{m\[Infinity],data,\[Omega]sample,mi,\[Tau]i,errorL2},
m\[Infinity] = Limit[Ppoly,var-> 0]//N;\[Omega]sample=Table[10^i,{i,explims[[1]],explims[[2]],0.01}];
data=  Join[Transpose[{\[Omega]sample}],Transpose[{Re[Ppoly/.{var-> I \[Omega]sample}]-m\[Infinity]}],Transpose[{Im[Ppoly/.{var-> I \[Omega]sample}]}],2];
{mi,\[Tau]i,errorL2}=DVEcouples[data,n,eps];{m\[Infinity],mi,\[Tau]i,errorL2}]

(*FunTest[coef_,var_]:=Block[{op},op = coef * ReplaceAll[var\[Rule] {2,3,4}][var];op]*)

End[]
EndPackage[]










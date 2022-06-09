(* ::Package:: *)

(* ::Section:: *)
(*Package for tensor calculus*)


(* ::Text:: *)
(*This first package carries with tensor calculus operations in  both harmonic and Hill's decomposition*)


(* ::Input:: *)
(*(* Name: TensorCalc *)*)
(*(* Author: Camilo Andr\[EAcute]s SUAREZ AFANADOR *)*)
(*(* Date : May 4,2021 *)*)
(* *)


BeginPackage["TensorCalc`"]

Iso4thTensor::usage = 
"Gives a fourth order Isotropic tensor for the given isotropic invariants \"\[Alpha]\" and \"\[Beta]\" "
FromIsoCoef2HillCoef::usage = 
"Gives a six length vector with the equivalent  Hill's invariants for the given isotropic invariants \"\[Alpha]\" and \"\[Beta]\" "
TIso4thTensor::usage = 
"Gives a fourth order tensor for a given vector \"aHill\" of lenthg 6 containing the Hill's basis invariants and the integer \"dir\" of the index of the principal direction (1, 2 or 3)"
Contr4thTenHillBasis::usage = 
"Gives the invariants of a 4th order tensor by contraction of two 4th order tensors given their invariants vectors \"Hill1\" and \"Hill2\" "
Contr4thAND2ndTen::usage = 
"Gives a second order tensor resulting from contraction of a 4th and 2nd order tensors \"Ten4th\" and \"Ten2nd\" "
Contr4thTen::usage = 
"Gives a fourth order tensor resulting from contraction of two 4th order tensors \"Ten1\" and \"Ten2\" "
From4thTenToHillCoef::usage = 
"Gives the vector of Hill's invariants of a 4th order transverse-isotropic tensor \"Ten4th\" with principal direction \"dir\" (1, 2 or 3)"
Inv4thHillCoef::usage = 
"Gives the vector of Hill's invariants of the inversion of a 4th order transverse-isotropic tensor \"Ten4th\" "
IsoDecomp::usage = 
"Function given the basis Tensors of the harmonic decomposition"
FrTenToMat::usage = 
"FrTenToMat[Ten] gives the matrix representation of a 4th order tensor in voigt notation as Parnell 2016 describe"
Rot4thTenMatNota::usage = 
"Rot4thTenMatNota[Mat,\[Theta]x,\[Theta]y,\[Theta]z] computes the rotation of a 4th order tensor in 6x6 notation for angles \[Theta]x,\[Theta]y and \[Theta]z in radians"
Rot2ndTenAxis::usage = 
"Rot2ndTenAxis[Mat,ANG,dir] computes the rotation of the second order tensor Mat around the direction 1,2 or 3  by  an angle ANG in radians"
FrMatToTen::usage = 
"FrMatToTen[Mat] give the fourth order Tensor from a 6x6 matrix representation of the tesor holding minor symmetries"
Inv4thMatNot::usage = 
"Inv4thMatNot[Mat] Gives the inverse of a 6x6 matrix representation (Mat) of a 4th order tensor
holding minor symmetries"
ProdQuadContrTens::usage = 
"ProdQuadContrTens[Ten1,Ten2] gives the quadruple contracted product of two fourth order tensors \"Ten1\" and \"Ten2\" "
Norm4thTen::usage = 
"Norm4thTen[Ten] gives the norm of a fourth order tensor \"Ten\" "
FromHillisoCoef2IsoCoef::usage = 
"FromHillisoCoef2IsoCoef[Hillc] gives  the two isotropic coefficients from Hill's basis coefficients representation of an isotropic 4th order tensor"


Begin["`Private`"]

(* 1.Harmonic decomposition for Tensorial quantities *)
(* Kronecker delta function *)
\[Delta]ij[i_,j_] := If[i==j,1,0]
(* 2nd order Identity tensor *)
Iij = Table[\[Delta]ij[i,j],{i,1,3},{j,1,3}];
(* 4th order Identity tensor *)
Iijkl =Table[1/2 (\[Delta]ij[i,k]\[Delta]ij[j,l]+\[Delta]ij[i,l]\[Delta]ij[j,k]),{i,1,3},{j,1,3},{k,1,3},{l,1,3}];
(* Harmonic decomposition *)
(* 4th Order Hydrostatic tensor *)
Jijkl=Table[ 1/3 (\[Delta]ij[i,j]\[Delta]ij[k,l]),{i,1,3},{j,1,3},{k,1,3},{l,1,3}];
(* 4th Order Deviatoric tensor *)
Kijkl = Iijkl - Jijkl;
(* 4th Order Isotropic Tensor in Cartesian basis *)
Iso4thTensor[\[Alpha]_,\[Beta]_]:= \[Alpha] Jijkl + \[Beta] Kijkl
IsoDecomp[]:={Jijkl,Kijkl}
(* 2.Hill's decomposition in transverse-isotropy for 4th Oder Tensorial quantities as a function of the principal direction (dir = 1|2|3 ) *)
HillBasis[dir_?IntegerQ]:=Block[{\[CapitalTheta]ij,H1,H2,H3,H4,H5,H6},
\[CapitalTheta]ij=Table[(\[Delta]ij[i,j]-\[Delta]ij[i,dir] \[Delta]ij[j,dir]),{i,3},{j,3}];
H1=Table[1/2 (\[CapitalTheta]ij[[i,j]] \[CapitalTheta]ij[[k,l]]),{i,3},{j,3},{k,3},{l,3}];
H2=Table[(\[CapitalTheta]ij[[i,j]] \[Delta]ij[k,dir] \[Delta]ij[l,dir]),{i,3},{j,3},{k,3},{l,3}];
H3=Table[(\[CapitalTheta]ij[[k,l]] \[Delta]ij[i,dir] \[Delta]ij[j,dir]),{i,3},{j,3},{k,3},{l,3}];
H4=Table[(\[Delta]ij[i,dir] \[Delta]ij[j,dir] \[Delta]ij[k,dir] \[Delta]ij[l,dir]),{i,3},{j,3},{k,3},{l,3}];
H5=Table[1/2 (\[CapitalTheta]ij[[i,k]] \[CapitalTheta]ij[[l,j]]+\[CapitalTheta]ij[[i,l]] \[CapitalTheta]ij[[k,j]]-\[CapitalTheta]ij[[i,j]] \[CapitalTheta]ij[[k,l]]),{i,3},{j,1,3},{k,3},{l,3}];
H6=Table[1/2 (\[CapitalTheta]ij[[i,k]] \[Delta]ij[l,dir] \[Delta]ij[j,dir]+\[CapitalTheta]ij[[i,l]] \[Delta]ij[k,dir] \[Delta]ij[j,dir]+\[CapitalTheta]ij[[j,k]] \[Delta]ij[l,dir] \[Delta]ij[i,dir]+\[CapitalTheta]ij[[j,l]] \[Delta]ij[k,dir] \[Delta]ij[i,dir]),{i,3},{j,3},{k,3},{l,3}];
{H1,H2,H3,H4,H5,H6}]
(* Fourth order tensor from Hill's decomposition coefficients aHill = {a1,a2,a3,a4,a5,a6} *)
TIso4thTensor[aHill_,dir_]:=Block[{H1,H2,H3,H4,H5,H6,A},
{H1,H2,H3,H4,H5,H6} =HillBasis[dir] ;
A = aHill[[1]] H1+aHill[[2]] H2+aHill[[3]] H3+aHill[[4]] H4+aHill[[5]] H5+aHill[[6]] H6;
A]

(* 3.Functions to pass from isotropic coefficients {\[Alpha],\[Beta]} to equivalent Hill's coefficients with normal dir = 3 *)
FromIsoCoef2HillCoef[\[Alpha]_,\[Beta]_]:={(2\[Alpha]+\[Beta])/3,(\[Alpha]-\[Beta])/3,(\[Alpha]-\[Beta])/3,(\[Alpha]+2\[Beta])/3,\[Beta],\[Beta]}  
FromHillisoCoef2IsoCoef[Hillc_]:={(3Hillc[[1]]-Hillc[[5]])/2,Hillc[[5]]}  

(* 4.Function to obtain the Invariants of 4th order Transverse isotropic tensor to Hill basis coefficients *)
(* 6x6 Matrix Writing of a 4th Order tensor holding minor symmetries in cartesian coordinates *)
FrTenToMat[Ten_]:=Block[{TArr,T,MatFrTen,FlatTen,FlatTen2,n,RL},TArr = Array[T,{3,3,3,3}];MatFrTen = {{TArr[[1,1,1,1]],TArr[[1,1,2,2]],TArr[[1,1,3,3]],TArr[[1,1,2,3]],
TArr[[1,1,3,1]],TArr[[1,1,1,2]]},{TArr[[2,2,1,1]],TArr[[2,2,2,2]],TArr[[2,2,3,3]],TArr[[2,2,2,3]],TArr[[2,2,1,3]],TArr[[2,2,1,2]]},{TArr[[3,3,1,1]],TArr[[3,3,2,2]],TArr[[3,3,3,3]],
TArr[[3,3,2,3]],TArr[[3,3,1,3]],TArr[[3,3,1,2]]},{TArr[[2,3,1,1]],TArr[[2,3,2,2]],TArr[[2,3,3,3]],TArr[[2,3,2,3]],TArr[[2,3,1,3]],TArr[[2,3,1,2]]},{TArr[[1,3,1,1]],TArr[[1,3,2,2]],
TArr[[1,3,3,3]],TArr[[1,3,2,3]],TArr[[1,3,1,3]],TArr[[1,3,1,2]]},{TArr[[1,2,1,1]],TArr[[1,2,2,2]],TArr[[1,2,3,3]],TArr[[1,2,2,3]],TArr[[1,2,1,3]],TArr[[1,2,1,2]]}};
FlatTen =Flatten[Ten];FlatTen2 =Flatten[TArr]; n= Length[FlatTen];RL=Table[0,{i,n}];Do[RL[[i]]=FlatTen2[[i]]-> FlatTen[[i]],{i,n}];MatFrTen/.RL]
(**)
FrMatToTen[Mat_]:=Block[{TArr,T,MatFrTen,FlatTen,FlatTen2,n,Rl,List1,List2},
TArr = Array[T,{3,3,3,3}];MatFrTen = {{TArr[[1,1,1,1]],TArr[[1,1,2,2]],TArr[[1,1,3,3]],TArr[[1,1,2,3]],
TArr[[1,1,3,1]],TArr[[1,1,1,2]]},{TArr[[2,2,1,1]],TArr[[2,2,2,2]],TArr[[2,2,3,3]],TArr[[2,2,2,3]],TArr[[2,2,1,3]],
TArr[[2,2,1,2]]},{TArr[[3,3,1,1]],TArr[[3,3,2,2]],TArr[[3,3,3,3]],TArr[[3,3,2,3]],TArr[[3,3,1,3]],TArr[[3,3,1,2]]},
{TArr[[2,3,1,1]],TArr[[2,3,2,2]],TArr[[2,3,3,3]],TArr[[2,3,2,3]],TArr[[2,3,1,3]],TArr[[2,3,1,2]]},{TArr[[1,3,1,1]],
TArr[[1,3,2,2]],TArr[[1,3,3,3]],TArr[[1,3,2,3]],TArr[[1,3,1,3]],TArr[[1,3,1,2]]},{TArr[[1,2,1,1]],TArr[[1,2,2,2]],
TArr[[1,2,3,3]],TArr[[1,2,2,3]],TArr[[1,2,1,3]],TArr[[1,2,1,2]]}};
List1 = Flatten[MatFrTen];List2 = Flatten[Mat]; Rl = 0*List1; Do[Rl[[i]] = List1[[i]]-> List2[[i]],{i,Length[List1]}];
TArr = TArr/.Rl; Do[TArr[[i,j,l,k]] = TArr[[j,i,k,l]] =TArr[[i,j,k,l]],{i,3},{j,3},{k,3},{l,3}];T[1,1,1,3] = 0;
TArr]

(* From 4th to Hill Coefs with principal direction dir = 1,2 or 3 *)
From4thTenToHillCoef[Ten4th_,dir_]:=Block[{Mat,Hcoef,l2,l3,mu},Mat=FrTenToMat[Ten4th];
If[dir==3, Hcoef={2(Mat[[1,1]]-Mat[[6,6]]), Mat[[1,3]],Mat[[3,1]],Mat[[3,3]],2Mat[[6,6]],2Mat[[4,4]]}];
If[dir==2, Hcoef={2(Mat[[1,1]]-Mat[[5,5]]), Mat[[1,2]],Mat[[2,1]],Mat[[2,2]],2Mat[[5,5]],2Mat[[4,4]]}];
If[dir==1, Hcoef={2(Mat[[2,2]]-Mat[[4,4]]), Mat[[1,2]],Mat[[2,1]],Mat[[1,1]],2Mat[[4,4]],2Mat[[6,6]]}];
Hcoef]

(* 5.Contractions of tensors *)
(* Contraction between 4th order tensors *)
Contr4thTen[Ten1_,Ten2_]:=Block[{Ten},Ten = Table[0,{i,3},{j,3},{k,3},{l,3}];
Do[Ten[[i,j,k,l]]=Sum[Sum[Ten1[[i,j,m,n]]Ten2[[n,m,k,l]],{m,3}],{n,3}],{i,3},{j,3},{k,3},{l,3}];
Ten]
(* Contraction between 4th order tensors in Hill Basis *)
Contr4thTenHillBasis[Hill1_,Hill2_]:={(Hill1[[1]]*Hill2[[1]]+2Hill1[[2]]*Hill2[[3]]),(Hill1[[1]]*Hill2[[2]]+Hill1[[2]]*Hill2[[4]]),
(Hill1[[3]]*Hill2[[1]]+Hill1[[4]]*Hill2[[3]]),(2Hill1[[3]]*Hill2[[2]]+Hill1[[4]]*Hill2[[4]]),(Hill1[[5]]*Hill2[[5]]),(Hill1[[6]]*Hill2[[6]])}
(* Contraction between 4th and 2nd order tensors  *)
Contr4thAND2ndTen[Ten4th_,Ten2nd_]:= Block[{OUT},OUT=Table[0,{i,3},{j,3}];
Do[OUT [[i,j]]=Sum[Sum[Ten4th[[i,j,k,l]]Ten2nd[[k,l]],{k,3}],{l,3}],{i,3},{j,3}];
OUT]
(* The norm of a 4th order tensor *)
ProdQuadContrTens[Ten1_,Ten2_]:=Block[{Norm2}, Norm2 = Sum[Sum[Sum[Sum[Ten1[[i,j,k,l]]*Ten2[[i,j,k,l]],{i,3}],{j,3}],{k,3}],{l,3}];
Norm2]
Norm4thTen[Ten_]:= Block[{norm},norm = Sqrt[ProdQuadContrTens[Ten,Ten]];
norm]



(* 6.Inversion of 4th order tensors in Hill's basis *)
Inv4thHillCoef[aHill_] := Block[{\[CapitalDelta]Hill, coef}, \[CapitalDelta]Hill = (aHill[[1]] aHill[[4]])/2 - (aHill[[2]] aHill[[3]]);
coef = {aHill[[4]]/(2 \[CapitalDelta]Hill), -aHill[[2]]/(2 \[CapitalDelta]Hill), -aHill[[3]]/(2 \[CapitalDelta]Hill), 
aHill[[1]]/(2 \[CapitalDelta]Hill), 1/aHill[[5]], 1/aHill[[6]]}; coef]

(* Inversion of 4th order tensors in Matrix notation *)
Inv4thMatNot[Mat_]:= Block[{W,M,InvMat}, W = {{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},
{0,0,0,2,0,0},{0,0,0,0,2,0},{0,0,0,0,0,2}}; M = Inverse[W];
InvMat = M . Inverse[Mat] . M;
InvMat]

(* Rotations of 4th and 2nd order tensors in voigt representation *)

(* 4th order tensors *)
GenRot[ANG_] := {{Cos[ANG]^2 ,Sin[ANG]^2,0,0,0,Sin[2ANG]},{Sin[ANG]^2 ,Cos[ANG]^2,0,0,0,-Sin[2ANG]},
{0,0,1,0,0,0},{0,0,0,Cos[ANG],-Sin[ANG],0},{0,0,0,Sin[ANG],Cos[ANG],0},{-1/2 Sin[2ANG],1/2 Sin[2ANG],0,0,0,Cos[2ANG]}}
GenRotP[ANG_] := {{Cos[ANG]^2 ,0,Sin[ANG]^2,0,-Sin[2ANG],0},{0,1,0,0,0,0},{Sin[ANG]^2 ,0,Cos[ANG]^2,0,Sin[2ANG],0},
{0,0,0,Cos[ANG],0,Sin[ANG]},{1/2 Sin[2ANG],0,-1/2 Sin[2ANG],0,Cos[2ANG],0},{0,0,0,-Sin[ANG],0,Cos[ANG]}}
RotOP4thTenMat6x6[\[Theta]x_,\[Theta]y_,\[Theta]z_]:= GenRot[\[Theta]x] . GenRotP[\[Theta]y] . GenRot[\[Theta]z]

(* 2nd order tensor around Subscript[e, 1],Subscript[e, 2] or Subscript[e, 3] *)
RotMat3Dy[ANG_,dir_] := Block[{vecdir, RotMat},vecdir = {0,0,0}; vecdir[[dir]]=1;RotMat=RotationMatrix[ANG,vecdir];
RotMat]

(* Rotation operations in 4th and 2nd order tensor in matrix representations  *)

Rot4thTenMatNota[Mat_,\[Theta]x_,\[Theta]y_,\[Theta]z_]:=RotOP4thTenMat6x6[\[Theta]x,\[Theta]y,\[Theta]z] . Mat . Transpose[RotOP4thTenMat6x6[\[Theta]x,\[Theta]y,\[Theta]z]]
Rot2ndTenAxis[Mat_,ANG_,dir_]:=RotMat3Dy[ANG,dir] . Mat . Transpose[RotMat3Dy[ANG,dir]]


End[]
EndPackage[]

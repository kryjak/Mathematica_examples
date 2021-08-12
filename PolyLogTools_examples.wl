(* ::Package:: *)

$PolyLogPath = SetDirectory["/home/jkrys/packages/PolyLogTools"]
<<PolyLogTools`;


HPLToG[HPL[{2},x]]


ShuffleG[G[a,b,z]*G[c,d,z]]


HToG[H[0,-1,1,2,z]]
GToLi[G[0,-1,-1,1,z]]


?*Fibration*


LiToG[Li[{1},{z}]]


G[0,z]


-PolyLog[1,1-z]
LiToG[-Li[{1},{1-z}]]
ToFibrationBasis[G[1,1-z], {z}] 


(* input the expression to be tested and convert to Li notation *)
(* PolyLog[m, z] = Li[{m},{z}] with Log[z] = -PolyLog[1, 1-z] *)
PolyLog[2,z]+PolyLog[2,1-z]+Log[1-z]*Log[z]-Pi^2/6 /. {PolyLog[2,z_]:>Li[{2},{z}], Log[z_]:>Li[{1},{1-z}]}
(* now convert to G language *)
LiToG[Li[{2},{z}]+Li[{2},{1-z}]+Li[{1},{z}]*Li[{1},{1-z}]-Pi^2/6]
(* fibration basis *)
G[1,1-z]*G[1,z]-G[0,1,1-z]-G[0,1,z]-Pi^2/6 // ToFibrationBasis[#,{z}]&
(* evaluate the product *)
ShuffleG[%]


coeff = Coefficient[diff //.ids //. polylogids, eps, 0] /. {spzA[1]->287/2000,spzA[2]->3287/10000,spzA[3]->3247/5000,spzB[1]->200/287} ;
coeff = coeff /. PolyLog[2,x_]:>Li[{2},{x}]/. {PolyLog[2,z_]:>Li[{2},{z}], Log[z_]:>Li[{1},{1-z}]} ;


LiToG[coeff] /.(f:(s|z))[x__]:>ToExpression@StringJoin[ToString[f], StringJoin@@IntegerString/@{x}];
ToFibrationBasis[%,{z1,z2,z3,s12,s23,s45,s5}, FitValue->{z1->1/10,z2->3/4,z3->3/20,s12->1/100,s23->3/111,s45->15,s5->10}];
Collect[ShuffleG[%/.{AAA:>-1, nc:>1}], {G[__],z1,z2,z3}, Simplify]


DecomposeToLyndonWords[G[1,0,0,z]]
ExtractZeroes[G[1,0,0,z]]
ShuffleG[DecomposeToLyndonWords[G[b,a,a,z], Alphabet->{a,b}]]


ShuffleG[G[0,z1] G[1,z1]]

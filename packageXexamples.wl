(* ::Package:: *)

(* ::Input:: *)
(*(*<<FeynCalc`*)*)
(*<< X`*)


(* ::Input:: *)
(*vacuumPol =-e^2*I^2/(16Pi^2)LoopIntegrate[Spur[\[Gamma] . k +m \[DoubleStruckOne], Subscript[\[Gamma], \[Nu]], \[Gamma] . (k+p)+m \[DoubleStruckOne], Subscript[\[Gamma], \[Mu]]], k, {k,m},{k+p,m}];*)
(*LoopRefine[vacuumPol] /. e^2->\[Alpha]*4Pi ;*)
(*\[CapitalPi]=Coefficient[%,1/\[Epsilon]]// Simplify (* + finite terms *)*)
(**)


(* ::Input:: *)
(*selfEnergy = -I^2e^2/(16Pi^2) LoopIntegrate[DiracMatrix[Subscript[\[Gamma], \[Mu]],\[Gamma] . (k+p)+m \[DoubleStruckOne],Subscript[\[Gamma], \[Mu]]],k,{k,0},{k+p,m}];*)
(*LoopRefine[selfEnergy] /. e^2->\[Alpha]*4Pi;*)
(*\[CapitalSigma] = Coefficient[%,1/\[Epsilon]] // Simplify (* + finite terms *)*)
(**)


(* ::Input:: *)
(*vertex = -I^2/(16Pi^2e)*e^3 LoopIntegrate[DiracMatrix[Subscript[\[Gamma], \[Nu]],\[Gamma] . (-k-p1)+m \[DoubleStruckOne], Subscript[\[Gamma], \[Mu]], \[Gamma] . (-k+p2 )+m \[DoubleStruckOne], Subscript[\[Gamma], \[Nu]]],k, {k,0},{k+p1,m},{k-p2,m}] ;*)
(*LoopRefine[vertex] /. e^2 -> \[Alpha]*4Pi ;*)
(*\[CapitalLambda] = Coefficient[%,1/\[Epsilon]] // Simplify (* + finite terms *)*)


(* ::Input:: *)
(*?LoopIntegrate*)
(*?DiracMatrix*)


(* ::Input:: *)
(*LoopIntegrate[k . p,k,{k,0},{k-p,0}]//LoopRefine*)
(*LoopIntegrate[1,k,{k,0},{k-p,0}]//LoopRefine*)


(* ::Input:: *)
(*?PVB*)


(* ::Input:: *)
(*LoopIntegrate[1,k,{k,M}] - LoopIntegrate[1,k,{k,m}]-(p^2-m^2+M^2)LoopIntegrate[1,k,{k,M},{k+p,m}] *)
(*% // LoopRefine  *)


(* ::Input:: *)
(*LoopIntegrate[1,k,{k,0}] // LoopRefine*)


(* ::Input:: *)
(*ScalarC0[0,2p1 . p2,0,m,m1,m2] // ExpandScalarC0*)


(* ::Input:: *)
(*LoopRefine[ScalarC0[m,2p1 . p2,1,0,0,0], ExplicitC0->All]*)


(* ::Input:: *)
(*PVB[0,0,s,0,0] // LoopRefine*)


(* ::Input:: *)
(*y =1/\[CapitalLambda]*Integrate[1/Sqrt[a^2-1/\[CapitalLambda]],a]*)


(* ::Input:: *)
(*x := -1/Sqrt[\[CapitalLambda]]Tanh[Sqrt[\[CapitalLambda] t + c]]*)
(*(D[x,t]/x)^2 // FullSimplify*)

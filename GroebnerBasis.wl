(* ::Package:: *)

(* ::Subsection:: *)
(*Basics of GroebnerBasis*)


?GroebnerBasis
?PolynomialReduce


(* calculate the Grobner basis of some polynomials *)
polys = {2x^2+3y^2,xy-4x}
gb = GroebnerBasis[polys,{x,y}]
 (* need to specify which variable is 'greater' for the lexicographic ordering *)
(* we can also use e.g. MonomialOrder->DegreeReverseLexicographic *)


(* define some polynomial to be reduced *)
pol = x^3*y^2 + x*y^7 +y^3


(* if we reduce wrt to polys, the reduction result depends on the path - i.e. whether we reduce wrt p1, then p2, or p2, then p1 *)
PolynomialReduce[pol,polys[[1]],{x,y}][[2]];
PolynomialReduce[%,polys[[2]],{x,y}][[2]]//Expand
PolynomialReduce[pol,polys[[2]],{x,y}][[2]];
PolynomialReduce[%,polys[[1]],{x,y}][[2]]//Expand


(* reduction wrt to both polynomials in polys will select one of these paths *)
PolynomialReduce[pol,polys,{x,y}]
(* check *)
(%[[1]] . polys+%[[2]]//Expand)==pol


(* if we reduce wrt to Groebner basis, the result is the same no matter which path we take *)
PolynomialReduce[pol,gb[[1]],{x,y}][[2]];
PolynomialReduce[%,gb[[2]],{x,y}][[2]]//Expand
PolynomialReduce[pol,gb[[2]],{x,y}][[2]];
PolynomialReduce[%,gb[[1]],{x,y}][[2]]//Expand


(* and the reduction wrt to both polynomials in gb is consistent *)
PolynomialReduce[pol,gb,{x,y}]//Expand
(* check *)
(%[[1]] . gb+%[[2]]//Expand)==pol


(* ::Subsection:: *)
(*MultivariateApart*)


Needs["MultivariateApart`"]


?MultivariateApart`*


(* expression to be partial-fractioned *)
rr = (2y-x)/(y(x+y)(y-x))
(* extract denominator factors *)
denfacs = Denominator[rr]/.x_*y__:>{x,y}
(* define inverse denominators - I think the ordering of q1, q2, q3,... doesn't matter as long as we are consistent throughout later steps *)
invdens = ToExpression/@Table["q"<>ToString[i],{i,Length[denfacs]}]
(* define rules for later to transform the reduced polynomial back into normal notation *)
rules = Thread[RuleDelayed[invdens,#]]&@(1/denfacs)


(* calculate monomial block ordering following Algorithm 3 in the paper *)
blockordering = ApartOrder[denfacs,invdens]


(* Groebner basis of the denominator factors using the block-ordering *)
gb = ApartBasis[denfacs,invdens,blockordering]
(* reduce wrt to the Groebner basis and block-ordering *)
(* polynomial is supplied in the 'inverse' form, i.e. Num*q1*q2*q3... *)
ApartReduce[Numerator[rr]*Times@@invdens,gb,blockordering]
ans = %/.rules
(* the MultivariateApart command achieves the same result in one line *)
ans - MultivariateApart[rr]//Simplify


(* we can also tryo to perform the reduction using the built-in Mathematica commands *)
(* first, we need to define the ideal using denfacs and invdens *)
ideal = Table[invdens[[i]]*denfacs[[i]]-1,{i,Length[denfacs]}]
(* we also need to convert the block-ordering into a weight matrix - HOW IS THIS DONE? *)
weights = DegRevLexWeights[blockordering]


gb2 = GroebnerBasis[ideal,{q1,q2,q3,x,y},MonomialOrder->weights]
PolynomialReduce[Numerator[rr]*Times@@invdens,gb2,{q1,q2,q3,x,y},MonomialOrder->weights][[2]]/.rules
(* this is the same as the result from the MultivariateApart package, but looks different *)
(* am I doing something wrong or is it some subtlety of MMA's Groebner Basis function? *)
%-ans//Simplify


(* if we use the Groebner Basis from MMA back in MultivariateApart, it gives the answer in the same form: *)
ApartReduce[Numerator[rr]*Times@@invdens,gb2,blockordering]/.rules
(* so it seems there is some subtlety of MMA's PolynomialReduce *)

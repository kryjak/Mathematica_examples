(* ::Package:: *)

(* ::Input:: *)
(* *)


(* ::Section::Closed:: *)
(*Slot #*)


Array[Plus[##] &, {2, 2}] 


Array[Plus[#1,#2] &, {2, 2}] 


Array[Plus[##,#2] &, {2, 2}] 


Array[Plus[## #2] &, {2,2}]


Plus[{2,2} 2]


(* ::Section::Closed:: *)
(*Plotting*)


RegionPlot[y<=1-x, {x,-.01,1},{y,-.01,1}, Axes:>True, AxesLabel->{x,y}]
RegionPlot3D[z<=1-x-y, {x,0,1},{y,0,1},{z,0,1}, AxesLabel->Automatic]


Plot[1-x, {x, 0, 1}, Filling -> Bottom, AxesLabel->{x,y}, AspectRatio->Automatic] 
Plot3D[1-x-y, {x, 0, 1},{y, 0, 1}, 
RegionFunction->Function[{x,y,z}, 0<=z<=1], 
PlotStyle->Directive[Opacity[0.5,Blue]],
Filling -> Bottom, 
FillingStyle->Directive[Opacity[0.1,Blue]],
(* Mesh\[RuleDelayed]None, *) (* The mesh on the surface can be turned off *)
BoxRatios:>{1,1,1}, 
AxesLabel->{Subscript[x,1],Subscript[x,2],Subscript[x,3]},
BoundaryStyle->Directive[Red, Thick]]


(* ::Section::Closed:: *)
(*Dual basis vectors*)


vec = {r*Sin[\[Theta]]*Cos[\[Phi]], r*Sin[\[Theta]]*Sin[\[Phi]], r*Cos[\[Theta]] }
rmag = Simplify[Norm[vec], Assumptions->{r>0,\[Theta]\[Element]Reals,\[Phi]\[Element]Reals}]


pars = {{r, Sqrt[x^2+y^2+z^2]},
		{\[Theta], ArcCos[z/Sqrt[x^2+y^2+z^2]]},
		{\[Phi], ArcTan[y/x]}
		}
coords = {
x:> r*Sin[\[Theta]]*Cos[\[Phi]],
y:> r*Sin[\[Theta]]*Sin[\[Phi]],
z:> r*Cos[\[Theta]]
}


e[x_] := D[vec, pars[[x,1]] ]
Do[Print[Subscript[Style[e,Bold], pars[[i,1]]], " = ", e[i]], {i,1,3}]


g = Table[e[i] . e[j], {i, 3}, {j, 3}] // Simplify ;
g // MatrixForm


f[a_] := Simplify[\!\(
\*SubscriptBox[\(\[Del]\), \({x, y, z}\)]\(pars[\([a, 2]\)]\)\) /. coords,
				  Assumptions->{r>0,0<\[Theta]<Pi,\[Phi]\[Element]Reals}
				  ] 
Do[Print[Style[e,Bold]^pars[[i,1]]," = ",f[i]], {i,1,3}]


check[i_] := Sum[ g[[i,j]]*f[j], {j,1,Length[g]} ] (* contract g_ij * e^j *)
Do[ Print[Subscript["e", pars[[i,1]]],(" = e")^pars[[i,1]]," ---> " ,check[i]===e[i]], {i,1,Length[g]} ]


(* ::Section::Closed:: *)
(*Tree form*)


TreeForm[1+a+b^2]


f[a_] := Sin[a(2*Pi-x)] 
g[a_] := Sin[-a*x]
Plot[{f[3], g[3]}, {x,-5,5}]


(* ::Section::Closed:: *)
(*Sound*)


EmitSound[Sound[SoundNote[{"Snare", "BassDrum"}]]]
Pause[0.08]
EmitSound[Sound[SoundNote["Snare"]]]
Pause[0.4]
Do[{EmitSound[Sound[SoundNote[{"RideCymbal", "Snare", "BassDrum2"}]]], 
Pause[0.04],
EmitSound[Sound[SoundNote["BassDrum"]]]
}, {20}]
Pause[1/2]
EmitSound[Sound[SoundNote[{"CrashCymbal", "BassDrum2"}]]]


(* ::Section::Closed:: *)
(*Inverting a system*)


SetDirectory[NotebookDirectory[]]
<< "result.txt";
data = %[[1,2,1;;3]]


xs = Transpose[data][[1]]
ys = Transpose[data][[2]]
data
coeffs = NSolve[
		{a/(xs[[1]]^2) + b/xs[[1]] + c == ys[[1]], 
         a/(xs[[2]]^2) + b/xs[[2]] + c == ys[[2]],
         a/(xs[[3]]^2) + b/xs[[3]] + c == ys[[3]]
        }, 
	{a,b,c}
	];
	
Plot[a/(x^2) + b/x + c /. coeffs, {x,-0.2, 0}]


(* ::Section::Closed:: *)
(*Numerical errors*)


f[x_] := 1/\[Pi] Cos[80 Sin[x] - x];

res = NIntegrate[f[x], {x, 0, 2Pi}, 
  Method -> {"Trapezoidal", "SymbolicProcessing" -> 0},
  IntegrationMonitor :> ((errors = Through[#1@"Error"]) &)];

res
Length[errors]
Total[errors]

(* 
  -0.112115          <-- integral value
  2.67841*10^-15     <-- error estimate by NIntegrate
  4.996*10^-16       <-- actual error
*)


f[x_] := 1/\[Pi] Cos[80 Sin[x] - x];

res2[a_] := NIntegrate[f[x], {x, 0, a*Pi},
 PrecisionGoal -> 7,
 IntegrationMonitor :> ((errors = Through[#1@"Error"]) &)]
res2[2]
errors;
Length[errors]
Total[errors] (* sums the errors for all subregions *)
(*
  -0.112115          <-- integral value
  175                <-- length of errors = number of subregions
  1.1067647*10^-8    <-- overall error estimate
*)



numvalues[{s12_,s23_}] := Module[{x,y,i},
x = ConstantArray[0, {Length[epslist], 2} ] ;
Do[x[[i]] = {epslist[[i]], res[s12,s23,epslist[[i]]]}, {i,1,Length[epslist]}];
{{s12, s23}, x}
]

num = Map[numvalues, points]
num[[1,2,1,2]]


f[x_] := 1/\[Pi] Cos[80 Sin[x] - x];

iRegionMethods = {"Axis", "Boundaries", "Dimension", "Error", 
  "GetRule", "Integral", "Integrand", "WorkingPrecision"}; 
 
res = Reap@NIntegrate[f[x], {x, 0, Pi}, PrecisionGoal -> 1.1, 
   Method -> "AdaptiveMonteCarlo", 
   IntegrationMonitor :> 
    Function[{iregs}, 
     Sow[Association /@ 
       Transpose[
        Map[Thread[# -> Through[iregs[#]]] &, iRegionMethods]]]]];
Dataset[Flatten[res[[2]]]]


numericalint00[\[Epsilon]_,s12_,s23_] :=
        (Gamma[1 - 2 \[Epsilon]] Gamma[
        2 + \[Epsilon]] (1/(-s23 x z +
        s12 y (-1 + x + y + z)))^(2 + \[Epsilon]))/(
        Gamma[1 - \[Epsilon]]^2 Gamma[1 + \[Epsilon]])

cut = 10^(-10)

numericalint[\[Epsilon]_,s12_,s23_] :=  NIntegrate[numericalint00[\[Epsilon],s12,s23], {x,cut,1-cut}, {y,cut,1-cut-x}, {z,cut,1-cut-x-y},
MinRecursion->10^5, MaxRecursion->10^12,
Method -> {"GlobalAdaptive", "MaxErrorIncreases" ->10^4},
PrecisionGoal->4,
WorkingPrecision->20,
IntegrationMonitor :> ((errors = Through[#1@"Error"]) &)]

res = numericalint[-1/20,-3,-1]
err = Total[errors]
Print[res]
Print[err]
SetDirectory[NotebookDirectory[]]
Save["test.txt", {res,err}]


numericalint00[\[Epsilon]_,s12_,s23_] :=
        (Gamma[1 - 2 \[Epsilon]] Gamma[
        2 + \[Epsilon]] (1/(-s23 x z +
        s12 y (-1 + x + y + z)))^(2 + \[Epsilon]))/(
        Gamma[1 - \[Epsilon]]^2 Gamma[1 + \[Epsilon]])

cut = 10^(-100)

numericalint[\[Epsilon]_,s12_,s23_] :=  NIntegrate[numericalint00[\[Epsilon],s12,s23], {x,cut,1-cut}, {y,cut,1-cut-x}, {z,cut,1-cut-x-y},
MinRecursion->10^5, MaxRecursion->10^12,
Method -> {"GlobalAdaptive", "MaxErrorIncreases" ->10^7, "SingularityHandler" -> "DuffyCoordinates"},
PrecisionGoal->4,
WorkingPrecision->200,
IntegrationMonitor :> ((errors = Through[#1@"Error"]) &)]

SetDirectory[NotebookDirectory[]]
res = numericalint[-1/10,-3/2,-3/2]
err = Total[errors]
Print[res]
Print[err]

Save["test.txt", {res,err}]


(* ::Section::Closed:: *)
(*Parallel evaluation*)


ParallelEvaluate[$ProcessID]
ParallelEvaluate[$MachineName]


Select[Range[2000], PrimeQ[2^#-1]&] // Timing


Parallelize[Select[Range[2000], PrimeQ[2^#-1]&]] // Timing


numericalint00[\[Epsilon]_,s12_,s23_] :=
        (Gamma[1 - 2 \[Epsilon]] Gamma[
        2 + \[Epsilon]] (1/(-s23 x z +
        s12 y (-1 + x + y + z)))^(2 + \[Epsilon]))/(
        Gamma[1 - \[Epsilon]]^2 Gamma[1 + \[Epsilon]])

cut = 10^(-10);


Head[cut]


DistributeDefinitions[numericalint00,cut];
ParallelEvaluate[Head[numericalint00]]


ClealAll
f[x_] := 1/\[Pi] Cos[80 Sin[x] - x]*Exp[-x^1.5]
test[a_] := NIntegrate[f[x], {x, 0, a*Pi}, PrecisionGoal -> 3]
Table[test[i], {i, 1, 10}] // Timing

Parallelize[Table[test[i], {i, 1, 10}]] // Timing


Table[numericalint00[i,-3,-1], {i,-0.9,0,0.2}] // Timing
Parallelize[Table[numericalint00[i,-3,-1], {i,-0.9,0,0.2}]] // Timing


numericalint[\[Epsilon]_,s12_,s23_] :=  NIntegrate[numericalint00[\[Epsilon],s12,s23], {x,cut,1-cut}, {y,cut,1-cut-x}, {z,cut,1-cut-x-y},
MinRecursion->10^5, MaxRecursion->10^12,
Method -> {"GlobalAdaptive", "MaxErrorIncreases" ->10^4, "SingularityHandler" -> "DuffyCoordinates"},
PrecisionGoal->4,
WorkingPrecision->20,
IntegrationMonitor :> ((errors = Through[#1@"Error"]) &)]


numericalint[-1/10,-3/2,-3/2]


ParallelEvaluate[numericalint[-1/10,-3/2,-3/2]]


f[x_] := x^2


Parallelize[f[4]]


(* ::Section::Closed:: *)
(*Similarity Transformation*)


(* ::Input:: *)
(*a = {{2,1},{-1,-1}} ;*)
(*s=Transpose[Eigenvectors[a]] ; (* Must be a transpose since the eigenvectors are given as row vectors rather than column ones *)*)
(*sinv = Inverse[s] ;*)
(*sinv . a . s // Simplify*)
(*Eigenvalues[a]*)
(*(* Eigenvalues and their order agree with the entries of the diagonalised matrix *)*)


(* ::Input:: *)
(*b ={{-2,5},{-1,3}} ;*)
(*t = Transpose[Eigenvectors[b]] ;*)
(*tinv = Inverse[t] ;*)
(*tinv . b . t // Simplify*)
(*Eigenvalues[b]*)


(* ::Input:: *)
(*p = s . tinv// Simplify*)
(*a . p - p . b (* This property is also true, see https://math.stackexchange.com/questions/625925/how-to-compute-the-similarity-transformation-matrix *)*)


(* ::Section::Closed:: *)
(*Group Theory Tut 1 Q13*)


ClearAll;
u = {{a,b},{-b\[Conjugate], a\[Conjugate]}} (*A matrix in the fundamental representation of SU(2) *)
ubar = u\[Conjugate] (* Antifundamental rep *)
tprod = ArrayFlatten[u\[TensorProduct]ubar]  (* tensor product - Subscript[D, 2\[TensorProduct]Overscript[2, _]] representation*)


id = {1,0,0,1} (* The identity in 2D leaves two vectors invariant, {1,0} and {0,1}, so in the combined representation, the vector {1,0,0,1} will also be invariant *)
tprod . id /. x_*x_\[Conjugate] :> Abs[x] /. (* Multiply id by the tensor product to check. replace aa^* = |a(|^2) *)
 Abs[a]+Abs[b] :> 1// MatrixForm (* Recall |a|^2 + |b(|^2) = 1 *)


sim = Transpose[{1/Sqrt[2]{1,0,0,1} ,1/Sqrt[2]{-1,0,0,-1},{0,1,0,0},{0,0,1,0}}];
MatrixForm[sim]


Inverse[sim] . tprod . sim


sim2=Transpose[Eigenvectors[tprod]] ; (* Must be a transpose since the eigenvectors are given as row vectors rather than column ones *)
sim2inv = Inverse[sim2] ;
diag=sim2inv . tprod . sim2  //. x_*x_\[Conjugate] :> Abs[x] /. (* Multiply id by the tensor product to check. replace aa^* = |a(|^2) *)
 Abs[a]+Abs[b] :> 1 


(* ::Text:: *)
(**)


(* ::Section::Closed:: *)
(*Avoid long expressions by using # *)


list = {{1,2,3},{3,4,5},{}};
(* declare the If statement as a function using &, then Map list onto the slots # *)
(* use {list} if a nested list is meant to be treated as one *) 
(* if list is not empty, return list *)
(* if list is empty, ##&[] replaces Null output with a 'vanishing function' *)
(* https://mathematica.stackexchange.com/questions/3700/how-to-avoid-returning-a-null-if-there-is-no-else-condition-in-an-if-construct *)
If[# != {}, #, ##&[]] & /@ list 
If[# != {}, #, ##&[]] & /@ {list}


(* ::Section::Closed:: *)
(*Save a no - argument command for later*)


dates ={};
update := AppendTo[dates,CurrentDate[]] (* needs to be := rather than = *)
update;
dates


Do[update,{i,10}];
dates


(* ::Text:: *)
(*Perhaps with an optional argument:*)


(* https://mathematica.stackexchange.com/questions/2651/how-to-pass-a-symbol-name-to-a-function-with-any-of-the-hold-attributes/2677#2677 *)
(* https://mathematica.stackexchange.com/questions/17767/how-to-modify-function-argument *)
Clear[list]
update2[list_List] := Append[list, 6]
list={};
update2[list]
update2[%]


(* ::Section::Closed:: *)
(*NumbersBases+module+error message+grid axes labels*)


(* based on https://reference.wolfram.com/language/howto/PutHeadingsInATable.html 
			https://reference.wolfram.com/language/workflow/SetUpErrorCheckingAndMessagesInAFunction.html *)

Clear[NumbersBases]
NumbersBases::usage = "This could be a documentation for how to use this function";
(* Create a message to be printed in case an invalid base is used *)
NumbersBases::maxbase = "Maximum base `1` is larger than the allowed value - capping at 36 <-> z";

(* x_:value indicates an optional argument with default 'value'. 
   Some arguments can be specified as optional without giving the default value as Mma has some default optoins built in 
   https://reference.wolfram.com/language/tutorial/Patterns.html#17673 *)
NumbersBases[nmin_:0, nmax_, maxbase_Integer] := Module[{bases,range,data,grid},
(* Range of bases used, change notation to hexadecimal after base 10 *)
If[maxbase<=10, bases=Range[2,maxbase],
	If[maxbase<=36, bases=Join[Range[2,10], Alphabet[][[;;maxbase-10]]],
		Message[NumbersBases::maxbase,maxbase];
		bases=Join[Range[2,10], Alphabet[]]
	  ]
];
range=nmax-nmin+1;
(* Populate the table with numbers nmin-nmax in bases 2-maxbases (possible bases are 2-36 *)
data = Table[BaseForm[i,j],{i,nmin,nmax},{j,2,Length[bases]+1}];
(* Prepend a first element to the list which contains all the horizontal headings (in this case, bases) *)
data = Prepend[data,bases];
(* Prepend another `list` with one element only - the name of the heading.
   Can ensure correct placement by prepending an empty list of length equal to half of the length of the horizontal headings.
   The -1 is just to shift it nicely into the middle of the table.
   Note that Table, Range etc. can be specified with a non-integer length, e.g. Range[10.9]=Range[10.1] - like an in-built Floor function. *)
data = Prepend[data,Join[Table["",Length[bases]/2-1],{"\nBase"}]];
(* Prepend the number to the each row (apart from the first two which are just the header 'Base' and a list of the bases *)
data = MapThread[Prepend,{data,Join[{"",""},Range[nmin,nmax]]}];
(* Prepend the table title and vertical header *)
data =  MapThread[Prepend,{data,Join[{"NUMBERS/ \n BASES"},Table["", (range+2)/2],{"Numbers"},Table["",Ceiling[range/2]-1]]}];
(* Display as a grid, with frames only for the two selected rows *)
(* Elements of {{a,b},{c,d}}\[Rule]True refer to: a/b = first/last rows, c/d = start/end columns *)
grid = Grid[data,{Frame->{None,None,{{{2,2},{2,-1}}->True,{{2,-1},{2,2}}->True}}}];
Return[grid];
]

?NumbersBases
NumbersBases[-2,9,39]


(* ::Input::Initialization:: *)
(* to convert from other bases to decimal, one should use base^^number *)
16^^8B54240883FA007706B800000000C383FA027706B801000000C353BB01000000B9010000008D041983FA0376078BD989


(* ::Section::Closed:: *)
(*Name+StringJoin+Capitalize*)


pos  = {{13,1,14,27,3}, {14,23,29,25}};
name = StringJoin /@ (Alphabet["Polish"][[#]] & /@ pos) // Capitalize (* related function: ToUpperCase/ToLowerCase *)


(* ::Section::Closed:: *)
(*Timing for very fast functions*)


(* If a naive attempt at timing a function returns time < $TimeUnit, then this result is not precise *)
$TimeUnit
Integrate[Sin[x],{x,-Pi,Pi}] // AbsoluteTiming
repeat = 10^3;
Do[Integrate[Sin[x],{x,-Pi,Pi}],repeat] // AbsoluteTiming 
%[[1]]/repeat (* average timing *)
(* this timing will be quicker than the original computation above, because the Kernel stores some 
   variables internally during repeatead evaluation for further use, so the subsequent computations are quicker *)


(* ::Section::Closed:: *)
(*Assign a value if symbol not already defined*)


(* ::Input::Initialization:: *)
ClearAll[condassign];
condassign[var_Symbol,value_:{1,2}]:=(var=value)
condassign[var_,value_:2020]:=var

Clear[a];                 
condassign[a]             (* Out: {1,2} because a had no value *)
condassign[a]             (* Out: {1,2} because a now has a value *)
b := Integrate[Sin[x],x]   (* Out: -Cos[x] because b had an assigned value *)
condassign[b]


(* ::Section::Closed:: *)
(*{{x1,x2}, {y1,y2}} --> {{x1,y1},{x2,y2}}*)


xvals = {a,b,c};
yvals = {d,e,f};
Riffle[xvals,yvals]
Partition[%,2]


(* ::Section::Closed:: *)
(*MaTeX*)


(* http://szhorvat.net/pelican/latex-typesetting-in-mathematica.html *)
(* https://github.com/szhorvat/MaTeX *)
(* Help\[Rule]Wolfram Documentation\[Rule]"MaTeX" *)


ResourceFunction["MaTeXInstall"][] (* install MaTeX *)


<<MaTeX`
ConfigureMaTeX[]
SetOptions[MaTeX, Magnification->1.2] (* default axis labels too small *)


"\\sum_{k=1}^\\infty \\frac{1}{k^2}" // MaTeX (* "" will take the content as literal TeX input - need \\ to escape \sin etc. *)
Sum[1/k^2,{k,1,Infinity}] // MaTeX            (* convert a Mathematica expression to TeX input *)
HoldForm[Sum[1/k^2,{k,1,Infinity}]] // MaTeX  (* need to use HoldForm to prevent Mathematica from evaluating the expression *)
MaTeX[Sin[a*x^2]]
MaTeX["\\sin(ax^2)"]                          (* the spacing is different *)               


texStyle = {FontFamily -> "Latin Modern Roman", FontSize -> 12}; (* use the LaTeX font *)

ContourPlot[x^2 + y^4 == 1, {x, -1.2, 1.2}, {y, -1.2, 1.2},
 BaseStyle -> texStyle,
 Epilog -> {
     Arrow[{{0.1, 0.3}, {0.5, 0.80}}],
     Inset[MaTeX["x^2+y^4=1", Magnification -> 2], {0.1, 0.3}, Scaled[{0.5, 1}]] (* can also specify Magnification locally *)
    }]
    
Plot[Sin[x], {x, 0, 2 Pi},
 Frame -> True, FrameStyle -> BlackFrame,
 FrameTicks -> {{Automatic, None},
                {Table[{x, MaTeX[x, "DisplayStyle" -> False]}, {x, Pi/4 Range[0, 8]}], None}},
 FrameLabel -> MaTeX /@ {"x", "\\sin x"}, (* \\ needed to escape \ *)
 BaseStyle -> texStyle]


(* MaTeX only supports the kind of TeX input that would be valid in the $...$ environment, i.e. inline. *)
MaTeX["
 \\begin{equation}
 x^2
 \\end{equation}
"]
(* However, "aligned" works: *)
MaTeX["
 \\begin{aligned}
 x+y&=z \\\\
 u+v&=w
 \\end{aligned}
"]


(* ::Input:: *)
(**)


(* ::Section::Closed:: *)
(*FiniteFlow*)


(* ::Subsection::Closed:: *)
(*FFAlgMul*)


<<FiniteFlow`


FFNewGraph[testgraph, input, {a,b,c}];
FFAlgRatFunEval[testgraph, vec1, {input}, {a,b,c}, {a,b,c}];
FFAlgRatFunEval[testgraph, vec2, {input}, {a,b,c}, {1,2,3}];

FFAlgMatMul[testgraph, matmul, {vec1, vec2}, 3,1,3];
FFGraphOutput[testgraph,matmul];

FFReconstructFunction[testgraph, {a,b,c}]


(* ::Subsection::Closed:: *)
(*FFAlgTake*)


<< FiniteFlow`
FFNewGraph[graph, in , {x, y}];
list = Table[x^i, {i, 20}]
FFAlgRatFunEval[graph, evalNode, {in}, {x, y}, list];
allCoeffs = f[#] & /@ Range[FFNParsOut[graph, evalNode]];
selectRules = {2, 6, 7, 15};
selectedCoeffs = allCoeffs[[selectRules]];
FFAlgTake[graph, indep, {evalNode}, {allCoeffs} -> selectedCoeffs];
FFGraphOutput[graph, indep];
FFReconstructFunction[graph, {x, y}]


eqs = {x + 2 y + 3 z == 0, 3 x + 4 y + 5 z == 0, 6 x + 7 y + 8 z == 0};
(*eqs = {x + 2 y + 3 z + 7 w== 4, 3 x + 4 y + 5 z +4w == 6, -19 x + 7 y + 8 z + 3w == 9, 12 x + 5 y + 8z + 11w == 0};*)
vars = Variables[eqs/.a_==b_:>a]
Solve[eqs, vars]


FFNewGraph[testgraph,input,vars];
FFAlgDenseSolver[testgraph,solve,{input},vars,eqs,vars];
FFGraphOutput[testgraph,solve];
(*FFSolverOnlyHomogeneous[testgraph,solve]*)
solverlearn = FFDenseSolverLearn[testgraph,vars]
{depvars, indepvars} = {"DepVars", "IndepVars"} /. solverlearn;
res = FFReconstructFunction[testgraph,vars]

Print["IndepEqs = ", FFSolverIndepEqs[testgraph,solve]]
Print["NIndepEqs = ", FFSolverNIndepEqs[testgraph,solve]]
Print["NParsOut = ", FFNParsOut[testgraph]]

resmat = ArrayReshape[res, {Length[depvars],Length[indepvars]+1}];
resmat // MatrixForm
Do[Print[depvars[[i]]," = ", resmat[[i,;;-2]] . indepvars + resmat[[i,-1]]],{i,Length[depvars]}];


(* ::Section::Closed:: *)
(*More complicated slots *)


Function[#, Select[Times@@Boole[FreeQ[#, #] &/@{Log[x],Log[x]^2}]==0]]
Map[Function[x, Select[x, Times@@Boole[FreeQ[x, #] &/@ {a,b,c}]==1]], {a,c,d}]
Function[x, Select[{a,c,d}, Times@@Boole[FreeQ[{a,c,d}, #] &/@ {a,b,c}]==1]]


(* ::Section::Closed:: *)
(*TensorContract*)


(* ::Input:: *)
(*(* the TensorContract[tensor, {{i,j},{k,l}...}] function contracts indices at positions {i,j}, {k,l},... of tensor(ijkl...) *)*)
(*(* e.g. To contract i and j of T(ijk), where i,j=1...n, we need to use TensorContract[T, {1,2}] *)*)
(*(* which leaves a rank-1 tensor with entries T` = {T(111)+T(221)+...+T(nn1), T(112)+T(222)+...+T(nn2) ,     ...n-times...     , T(11n)+T(22n)+...+T(nnn)}*)*)
(*(* e.g. if we have T(ijk), i,j,k=1,2, and contract i and j, we will obtain a two element rank-1 tensor (a vector) with elements: *)*)
(*(* T` = {T(111)+T(221), T(112)+T(222)} *)*)


(* ::Input:: *)
(*tensor = {{{a[1,1,1],a[1,1,2]},{a[1,2,1],a[1,2,2]}},{{a[2,1,1],a[2,1,2]},{a[2,2,1],a[2,2,2]}}};*)


(* ::Input:: *)
(*TensorContract[tensor,{1,2}]*)
(*%=={tensor[[1,1,1]]+tensor[[2,2,1]],tensor[[1,1,2]]+tensor[[2,2,2]]}*)


(* ::Input:: *)
(*(* of course the dimensionality of each tensor level needs to be the same *)*)
(*(* e.g. this will not work, because index i has depth 2, but index j has depth 3, so we cannot contract these two indices *)*)


(* ::Input:: *)
(*TensorContract[{{a,b,c},{c,d,e}},{{1,2}}]*)


(* ::Input:: *)
(*(* we can implement the matrix determinant using the Levi-Civita tensor *)*)
(*(* take a square matrix with rows W(1)...W(n) and form the outer product with \[Epsilon] in n-dim *)*)
(*(* we get a rank-2n tensor \[Epsilon](ijk...n)W(n+1 n+2 ... 2n) *)*)
(*(* contract it pairwise: i with n+1, j with n+2, n with 2n *)*)
(*mat = {{1,3},{4,5}};*)
(*Det[mat]*)
(*LeviCivitaTensor[2,List]\[TensorProduct]mat[[1]]\[TensorProduct]mat[[2]];*)
(*TensorContract[%,{{1,3},{2,4}}]*)


(* ::Input:: *)
(*mat = {{a,b,c},{d,e,f},{g,h,j}};*)
(*Det[mat]*)
(*LeviCivitaTensor[3,List]\[TensorProduct]mat[[1]]\[TensorProduct]mat[[2]]\[TensorProduct]mat[[3]];*)
(*TensorContract[%,{{1,4},{2,5},{3,6}}]*)


(* ::Subsubsection::Closed:: *)
(*tr5 is parity-odd*)


(* ::Input:: *)
(*(*tr5 = 4*i*eps(a,b,c,d)*p1(a)*p2(b)*p3(c)*p4(d)*)*)
(*moms = Table[p[i,j],{i,4},{j,4}]*)
(*tensor = LeviCivitaTensor[4,List]\[TensorProduct]moms[[1]]\[TensorProduct]moms[[2]]\[TensorProduct]moms[[3]]\[TensorProduct]moms[[4]];*)
(*tr5 = 4*I*TensorContract[tensor,{{1,5},{2,6},{3,7},{4,8}}];*)
(*tr5parity =4*I*TensorContract[tensor/.p[i_,j_]:>-p[i,j]/;j!=1,{{1,5},{2,6},{3,7},{4,8}}]//Simplify;*)
(*tr5-(-tr5parity)//Simplify*)


tr5 + (tr5 /. {p[1,x_]:>p[3,x], p[3,x_]:>p[1,x]}) // Simplify
tr5 + (tr5 /. {p[1,x_]:>p[2,x], p[2,x_]:>p[1,x]}) // Simplify


(* ::Subsubsection::Closed:: *)
(*tr5^2 = Gram determinant*)


(* ::Input:: *)
(*Vmat = moms//Transpose;*)
(*g = {{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}};*)
(*GramDet = Det[2*Transpose[Vmat] . g . Vmat]//Simplify;*)
(*tr5^2-GramDet//Simplify*)


(* ::Section::Closed:: *)
(*Colour-ordering at tree-level (Schwartz example)*)


(* ::Input:: *)
(*f[a_,b_,A,c_,d_,A]:=tr[a,b,c,d]-tr[b,a,c,d]-tr[a,b,d,c]+tr[b,a,d,c] (* Eq. 27.77 *)*)


(* ::Input:: *)
(*(* the full 4-gluon tree amplitude is the sum of 3 channels *)*)
(*(* Mstilde is the colour-stripped amplitude *)*)
(*(* t-channel related to s-channel by a 2<->4 swap, u-channel by a 2<->3 swap *)*)
(*Ms=f[1,2,A,3,4,A]*Mstilde[1,2,3,4]; (* Eq 27.78 *)*)
(*Mt=f[1,4,A,3,2,A]*Mstilde[1,4,3,2];*)
(*Mu=f[1,3,A,2,4,A]*Mstilde[1,3,2,4];*)


(* ::Input:: *)
(*(* cycle trace so that 1 is always first *)*)
(*tracerules = {tr[2,1,3,4]->tr[1,3,4,2],tr[2,1,4,3]->tr[1,4,3,2],tr[3,1,2,4]->tr[1,2,4,3],tr[3,1,4,2]->tr[1,4,2,3],tr[4,1,2,3]->tr[1,2,3,4],tr[4,1,3,2]->tr[1,3,2,4]};*)


(* ::Input:: *)
(*(* in the s-channel, we can use the antisymmetry of the swaps: 1<->2, 3<->4 *)*)
(*Collect[Ms+Mt+Mu /. tracerules , {tr[__]}] /.{-Mstilde[x___,3,y___,4,z___]-> Mstilde[x,4,y,3,z], -Mstilde[x___,4,y___,3,z___]-> Mstilde[x,3,y,4,z]}*)


(* ::Input:: *)
(*(* there are 6 traces, each come with 2 colour-stripped amplitudes in the s-channel *)*)
(*(* alternatively, we can rewrite each sum as Mstilde[1,j,k,l]+Mttilde[1,j,k,l] by using the 2<->4 swap again *) (* Eq. 27.81 *)*)
(*(* this means that the amplitude is written as: Sum[Mtilde[1,j,k,l]*tr[1,j,k,l], {permutations of {j,k,l}}], where Mtilde = Mstilde + Mttilde are the colour-ordered partial amplitudes *)*)
(*(* Hence, each trace receives a contribution from the partial amplitudes which contains a sum of PLANAR diagrams only within a given particle ordering *)*)


(* ::Section::Closed:: *)
(*Generating momentum twistors*)


(* ::Input:: *)
(*(* this follows Appendix A of [hep-ph/1310.1051] and shows how to generate 5pt massless kinematics *)*)


(* ::Input:: *)
(*(* start by defining a parametrisation of momentum twistors *)*)
(*(* there are many ways to do it and it's not clear which one will lead to most compact expressions *)*)
(*Zmat = {{1,0,1/x1,1/x1+1/x2,1/x1+1/x2+1/x3}, {0,1,1,1,1}, {0,0,0,x4,1}, {0,0,1,1,x5/x4}};*)


(* ::Input:: *)
(*Zmat//MatrixForm (* not there are 4N-N-10 = 5 variables (momenta invariant under Poincare - 10-dim symmetry, and little group scaling - U(1) for each momentum *)*)


(* ::Input:: *)
(*(* define holomorphic spinors *)*)
(*lambda[i_] := Thread@Zmat[[;;2,i]]*)
(*mu[i_] := Thread@Zmat[[3;;,i]]*)
(*mu[0] = mu[5];*)
(*mu[6] = mu[1];*)
(*lambda[6]=lambda[1];*)
(*lambda[0]=lambda[5];*)


(* ::Input:: *)
(*(* 2D Levi-Civita *)*)
(*epsmat = {{0,1},{-1,0}};*)


(* ::Input:: *)
(*(* define the angle bracket *)*)
(*angle[i_,j_] := lambda[i] . epsmat . lambda[j]*)


(* ::Input:: *)
(*(* define the anti-holomorphic spinors *)*)
(*lambdatilde[i_]:=(angle[i,i+1]*mu[i-1]+angle[i+1,i-1]*mu[i]+angle[i-1,i]*mu[i+1])/(angle[i,i+1]*angle[i-1,i])*)


(* ::Input:: *)
(*(* define the square bracket *)*)
(*square[i_,j_] := lambdatilde[i] . epsmat . lambdatilde[j]*)


(* ::Input:: *)
(*(* off by some minus signs, but that's the general idea *)*)
(*s[i_,j_]:=angle[i,j]*square[j,i]*)


(* ::Input:: *)
(*{s[1,2],s[2,3],s[3,4],s[4,5],s[1,5]} // Simplify*)


(* ::Section::Closed:: *)
(*Tensor decomposition of an amplitude *)


(* ::Subsection:: *)
(*5pt or more*)


(* the goal of the following is to derive the basis needed in tensor decomposition of a 5pt gluon amplitude, as described in [hep-ph/1906.03298] *)


(* for n>=5, the 4D space is spanned by the external momenta: *)
ps = p /@ Range[4]


(* for an n-pt amplitude, the tensor basis is composed of rank-n tensors T[mu1,mu2,...,mun] *)
(* so for n=5, we have 5 free indices *)
(* we can construct 3 possible types of monomials with exactly 5 indices: *)
(* 1. p1[mu1]p2[mu2]p3[mu3]p4[mu4]p5[mu5] - only external momenta *)
(* 2. p1[mu1]p2[mu2]p3[mu3]g[mu4,mu5] - one g[mu,nu] and 3 external momenta *)
(* 3. p1[mu1]g[mu2,mu3]g[mu4,mu5] - two g[mu,nu] and one external momentum *)
(* NOTE: the momenta here are drawn from the {p1,...p4} set, that is excluding p5 ! *)
(* we will now construct all the possible the monomials of each type *)


(* p1[mu1]p2[mu2]p3[mu3]p4[mu4]p5[mu5] *)
fiveps = Times@@Table[mu[i,#[[i]]],{i,Length[#]}]&/@Tuples[ps,5];


(* p1[mu1]p2[mu2]p3[mu3]g[mu4,mu5] *)
mus = Table[mu[i],{i,Range[5]}];
musforg = Subsets[mus, {2}];
gs = g@@@musforg;
musforp = Complement[mus,#]&/@musforg;
threeps1g = gs*Table[Times@@@(Times[ii,#]&/@Tuples[ps,3]/.mu[i_]*p[j_]->mu[i,p[j]]), {ii,musforp}]//Flatten;


(* p1[mu1]g[mu2,mu3]g[mu4,mu5] *)
mus = Table[mu[i],{i,Range[5]}];
musforg1 = Subsets[mus, {2}];
gs1 = g@@@musforg1;
musforg2 = Complement[mus,#]&/@musforg1;
musforg2 = Subsets[#,{2}]&/@musforg2;
gs2 = (g@@@#) &/@ musforg2;
gs1gs2 = Table[gs1[[ii]]*#&/@gs2[[ii]], {ii,Length[gs1]}]//Flatten;
muforp = Flatten[Complement[mus, #] &/@ Flatten[Table[Join[musforg1[[ii]],#]&/@musforg2[[ii]],{ii,Length@musforg1}], 1]/.mu[i_]:>(mu[i,#]&/@ps), 1];
onep2gs = gs1gs2*muforp//Flatten//DeleteDuplicates;
(* delete duplicates becase g[mu[1],mu[2]]*g[mu[3],mu[4]] = g[mu[3],mu[4]]*g[mu[1],mu[2]], so we're overcounting by 2 *)


Length/@{fiveps,threeps1g,onep2gs}
{4^5, 4^3*Binomial[5,2], 4*Binomial[5,2]*Binomial[3,2]/2}


(* so overall we have 1724 tensor basis elements *)
tensorbasis = Join[fiveps,threeps1g,onep2gs];
Length[tensorbasis]


(* but the basis is contracted with external polarisation vectors, see e.g. (Eq. 2.7), which will impose further constraints *)
epsilons = Times@@(mu[#,eps[p[#]]] &/@ Range[5])


(* multiply the basis by the epsilons, contract the indices and impose constraints *)
SetAttributes[dot,Orderless]
fixgauge = {
	mu[i_,eps[p[j_]]]*mu[i_,p[k_]]:>dot[eps[p[j]],p[k]], (* contract indices *)
	dot[eps[p[i_]],p[i_]]:>0, (* transversality condition eps(i).p(i)=0 *)
	x___*dot[eps[p[5]],p[4]]*y___->x*y*Table[dot[eps[p[5]],p[i]],{i,3}], (* needed to impose eps(5).p(5)=0 as p(5) does not appear explicitly *)
	Table[dot[eps[p[i]],p[i+1]]:>0, {i,Range[4]}], (* fix the gauge with the cyclic choice *) (* I don't understand why this is a gauge choice *)
	dot[eps[p[5]],p[1]]:>0,
	dot[eps[p[4]],p[3]]:>0 (* this is required to bring the number of basis elements from 185 to 142 - but where does it come from ??? *)
	} // Flatten;


(* after imposing the constraints, the basis has 142 elements *)
tensorbasis = epsilons*tensorbasis //. fixgauge // Flatten;
tensorbasis = DeleteCases[tensorbasis,0] // DeleteDuplicates;
tensorbasis//Length


(* but because these are physical projectors, i.e. for external states in 4D, we reject all terms depend on g[mu,nu] *)
(* now only 32 elements remain, in agreement with (Eq. 4.3) *)
basisfreeofg = Select[tensorbasis, FreeQ[Head/@Variables[#], g] &]//DeleteDuplicates;
basisfreeofg = basisfreeofg /. dot[eps[p[i_]],p[j_]]:>mu[i,p[j]];
basisfreeofg // Length


(* we can check to make sure these are the same as given in the ancilliary files of the paper *)
Tensors5g={
p1[mu2]*p1[mu3]*p1[mu4]*p2[mu5]*p3[mu1],p1[mu2]*p1[mu3]*p2[mu4]*p2[mu5]*p3[mu1],p1[mu2]*p1[mu3]*p1[mu4]*p3[mu1]*p3[mu5],p1[mu2]*p1[mu3]*p2[mu4]*p3[mu1]*p3[mu5],
p1[mu2]*p1[mu4]*p2[mu3]*p2[mu5]*p3[mu1],p1[mu2]*p2[mu3]*p2[mu4]*p2[mu5]*p3[mu1],p1[mu2]*p1[mu4]*p2[mu3]*p3[mu1]*p3[mu5],p1[mu2]*p2[mu3]*p2[mu4]*p3[mu1]*p3[mu5],
p1[mu3]*p1[mu4]*p2[mu5]*p3[mu1]*p4[mu2],p1[mu3]*p2[mu4]*p2[mu5]*p3[mu1]*p4[mu2],p1[mu3]*p1[mu4]*p3[mu1]*p3[mu5]*p4[mu2],p1[mu3]*p2[mu4]*p3[mu1]*p3[mu5]*p4[mu2],
p1[mu4]*p2[mu3]*p2[mu5]*p3[mu1]*p4[mu2],p2[mu3]*p2[mu4]*p2[mu5]*p3[mu1]*p4[mu2],p1[mu4]*p2[mu3]*p3[mu1]*p3[mu5]*p4[mu2],p2[mu3]*p2[mu4]*p3[mu1]*p3[mu5]*p4[mu2],
p1[mu2]*p1[mu3]*p1[mu4]*p2[mu5]*p4[mu1],p1[mu2]*p1[mu3]*p2[mu4]*p2[mu5]*p4[mu1],p1[mu2]*p1[mu3]*p1[mu4]*p3[mu5]*p4[mu1],p1[mu2]*p1[mu3]*p2[mu4]*p3[mu5]*p4[mu1],
p1[mu2]*p1[mu4]*p2[mu3]*p2[mu5]*p4[mu1],p1[mu2]*p2[mu3]*p2[mu4]*p2[mu5]*p4[mu1],p1[mu2]*p1[mu4]*p2[mu3]*p3[mu5]*p4[mu1],p1[mu2]*p2[mu3]*p2[mu4]*p3[mu5]*p4[mu1],
p1[mu3]*p1[mu4]*p2[mu5]*p4[mu1]*p4[mu2],p1[mu3]*p2[mu4]*p2[mu5]*p4[mu1]*p4[mu2],p1[mu3]*p1[mu4]*p3[mu5]*p4[mu1]*p4[mu2],p1[mu3]*p2[mu4]*p3[mu5]*p4[mu1]*p4[mu2],
p1[mu4]*p2[mu3]*p2[mu5]*p4[mu1]*p4[mu2],p2[mu3]*p2[mu4]*p2[mu5]*p4[mu1]*p4[mu2],p1[mu4]*p2[mu3]*p3[mu5]*p4[mu1]*p4[mu2],p2[mu3]*p2[mu4]*p3[mu5]*p4[mu1]*p4[mu2]};

SubsetQ[basisfreeofg /. mu[i_,p[j_]]:> ToExpression["p"<>ToString[j]][ToExpression["mu"<>ToString[i]]], Tensors5g]


(* ::Subsection:: *)
(*4pt example*)


(* for n=4, the 4D space is spanned by 3 of the external momenta: *)
ps = p /@ Range[3]


(* for an n-pt amplitude, the tensor basis is composed of rank-n tensors T[mu1,mu2,...,mun] *)
(* so for n=4, we have 4 free indices *)
(* we can construct 3 possible types of monomials with exactly 4 indices: (see Eq.3.3 of [hep-ph/0202266]) *)
(* 1. p1[mu1]p2[mu2]p3[mu3]p4[mu4] - only external momenta *)
(* 2. p1[mu1]p2[mu2]g[mu3,mu4] - one g[mu,nu] and 2 external momenta *)
(* 3. g[mu1,mu2]g[mu3,mu4] - two g[mu,nu] *)
(* NOTE: the momenta here are drawn from the {p1,...p3} set, that is excluding p4 ! *)
(* we will now construct all the possible the monomials of each type *)


(* p1[mu1]p2[mu2]p3[mu3]p4[mu4] *)
fourps = Times@@Table[mu[i,#[[i]]],{i,Length[#]}]&/@Tuples[ps,4];


(* p1[mu1]p2[mu2]g[mu3,mu4] *)
mus = Table[mu[i],{i,Range[4]}];
musforg = Subsets[mus, {2}];
gs = g@@@musforg;
musforp = Complement[mus,#]&/@musforg;
twops1g = gs*Table[Times@@@(Times[ii,#]&/@Tuples[ps,2]/.mu[i_]*p[j_]->mu[i,p[j]]), {ii,musforp}]//Flatten;


(* g[mu1,mu2]g[mu3,mu4] *)
mus = Table[mu[i],{i,Range[4]}];
musforg1 = Subsets[mus, {2}];
gs1 = g@@@musforg1;
musforg2 = Complement[mus,#]&/@musforg1;
musforg2 = Subsets[#,{2}]&/@musforg2;
gs2 = (g@@@#) &/@ musforg2;
gs1gs2 = Table[gs1[[ii]]*#&/@gs2[[ii]], {ii,Length[gs1]}]//Flatten // DeleteDuplicates
(*muforp = Flatten[Complement[mus, #] &/@ Flatten[Table[Join[musforg1[[ii]],#]&/@musforg2[[ii]],{ii,Length@musforg1}], 1]/.mu[i_]:>(mu[i,#]&/@ps), 1];*)
(*onep2gs = gs1gs2*muforp//Flatten//DeleteDuplicates;*)
(* delete duplicates becase g[mu[1],mu[2]]*g[mu[3],mu[4]] = g[mu[3],mu[4]]*g[mu[1],mu[2]], so we're overcounting by 2 *)


Length/@{fourps,twops1g,gs1gs2}
{3^4, 3^2*Binomial[4,2], Binomial[4,2]/2}


(* so overall we have 1724 tensor basis elements *)
tensorbasis = Join[fourps,twops1g,gs1gs2];
Length[tensorbasis]


(* but the basis is contracted with external polarisation vectors, which will impose further constraints *)
epsilons = Times@@(mu[#,eps[p[#]]] &/@ Range[4])


(* multiply the basis by the epsilons, contract the indices and impose constraints *)
SetAttributes[dot,Orderless]
fixgauge = {
	mu[i_,eps[p[j_]]]*mu[i_,p[k_]]:>dot[eps[p[j]],p[k]], (* contract indices *)
	dot[eps[p[i_]],p[i_]]:>0, (* transversality condition eps(i).p(i)=0 *)
	x___*dot[eps[p[4]],p[3]]*y___->x*y*Table[dot[eps[p[4]],p[i]],{i,2}], (* needed to impose eps(5).p(5)=0 as p(5) does not appear explicitly *)
	Table[dot[eps[p[i]],p[i+1]]:>0, {i,Range[4]}], (* fix the gauge with the cyclic choice *) (* I don't understand why this is a gauge choice *)
	dot[eps[p[4]],p[1]]:>0,
	dot[eps[p[3]],p[2]]:>0 (* this is required to bring the number of basis elements from 185 to 142 - but where does it come from ??? *)
	} // Flatten;


(* after imposing the constraints, the basis has 142 elements *)
tensorbasis = epsilons*tensorbasis //. fixgauge // Flatten;
tensorbasis = DeleteCases[tensorbasis,0] // DeleteDuplicates;
tensorbasis//Length


(* ::Section::Closed:: *)
(*Changing variables within FF*)


<<FiniteFlow`


(* ::Subsection::Closed:: *)
(*A general example*)


FFNewGraph[graphXY,in,{x,y}];
FFAlgRatFunEval[graphXY,evalXY,{in},{x,y},{x^4+y^8}];
FFGraphOutput[graphXY,evalXY];

(* {x,y}->{a^2,b^2} *)
FFNewGraph[graphAB,in,{a,b}];
FFAlgRatFunEval[graphAB,evalAB,{in},{a,b},{a^2,b^2}];
FFAlgSimpleSubgraph[graphAB,subnode,{evalAB},graphXY];

FFGraphOutput[graphAB,subnode];
FFReconstructFunction[graphAB,{a,b}]


(* ::Subsection::Closed:: *)
(*sij -> MTs example*)


process = "H3g";
psmode = "PSanalytic";


<< "InitTwoLoopToolsFF.m"
<<"~/gitrepos/myfiniteflowexamples/setupfiles/Setup_H3g.m"
Get["~/gitrepos/myfiniteflowexamples/amplitudes/GlobalMapProcessSetup.m"];


(* this is how we can perform a variable change within the FF setup *)
(* first, create a graph in terms of the original variables: *)
vvars = {s12,s23,s4};
FFNewGraph[sijgraph,in,vvars];
FFAlgRatFunEval[sijgraph,evalnode,{in},vvars,{(s12^4+s4^2)/s23}];
FFGraphOutput[sijgraph,evalnode];

(* then get the list of old variables in terms of the news ones *)
(* in this case, we're converting sijs to MTs *)
sijstomts = GetMomentumTwistorExpression[{s[1,2],s[2,3],s[4]},PSanalytic]
vvarsmt = Variables[sijstomts]

(* initiate a graph in the new variables *)
FFNewGraph[mtgraph,in2,vvarsmt];
(* evaluate the substitution as a node *)
FFAlgRatFunEval[mtgraph, newvars, {in2}, vvarsmt, sijstomts];
(* use that node as input for the subgraph where we evaluate the original graph in terms of sijs *)
FFAlgSimpleSubgraph[mtgraph,subnode,{newvars},sijgraph];
FFGraphOutput[mtgraph,subnode];

(* check *)
res = FFReconstructFunction[mtgraph,vvarsmt];
exprmt = GetMomentumTwistorExpression[{(s[1,2]^4+s[4]^2)/s[2,3]},PSanalytic];
res/exprmt // Simplify


(* ::Section::Closed:: *)
(*GetSlice*)


<<FiniteFlow`
Get["~/gitrepos/myfiniteflowexamples/TwoLoopTools/BCFWreconstructionTools/BCFWReconstructFunction.wl"];
Get["~/gitrepos/myfiniteflowexamples/TwoLoopTools/BCFWreconstructionTools/ReconstructFunctionApart.wl"];


expr = 1/s12;
var = {s12};
FFNewGraph[graph,in,var];
FFAlgRatFunEval[graph,evalnode,{in},var,{expr}];
FFGraphOutput[graph,evalnode];


(* get slice evaluates the expression on a univariate slice in xx *)
(* each graph variable is replaced with A + B*tt, where A,B - random primes *)
(* this sliced graph is then reconstructed - modulo FFPrimeNo[0] *)
{slicerules, coeffs} = GetSlice[graph,var]
(* we can check that the numbers in coeffs correspond to the slice variables by doing rational reconstruction: *)
s12 /. slicerules
coeffs /. x_Integer :> FFRatRec[x,FFPrimeNo[0]] // Simplify


process = "H3g";
psmode  = "PSanalytic";
<<"InitTwoLoopToolsFF.m"
Get["amplitudes/GlobalMapProcessSetup.m"];


(* ::Text:: *)
(*Now let's look at something more complicated. When using RecontructFunctionFactors, we need to caluclate the slice of the graph first.*)
(*Below we show how to understand the output of that slice.*)


FFDeleteGraph[graph]
expr = GetMomentumTwistorExpression[s[1,3]^3/(eps*s[4]^2), PSanalytic];
vars = {eps,ex[1],ex[2],ex[5]};
FFNewGraph[graph,in,vars];
FFAlgRatFunEval[graph,evalnode,{in},vars,{expr}];
FFGraphOutput[graph,evalnode];


coeffansatz = {s[1,2],s[2,3],s[1,3],s[4]};
ReconstructFunctionFactors[graph,Join[{eps},PSanalytic],PSanalytic,PSvalues,"CoefficientAnsatz"->coeffansatz];


(* the starting point is the FF graph *)
(* in this case, we know the function analytically: *)
expr = GetMomentumTwistorExpression[s[1,3]^3/(eps*s[4]^2),PSanalytic]
(* first, the univariate slice is calculated: *)
slice = {eps->33845041+64920803 xx,ex[1]->93589099+32934527 xx,ex[2]->23383127+29051563 xx,ex[5]->76172237+74713663 xx};
expr = expr /. slice;
(* here the analytic form is known, so it comes out in a factorised form, that's why I'm expanding it *)
expr = Expand[Numerator[expr]]/Expand[Denominator[expr]]
(* the slice is normalised such that the coefficient of the xx^0 term in the denominator is always 1 *)
(* remember that everything has to be evaluated mod p *)
norm = Coefficient[Denominator[expr],xx,0]
expr = Expand[Numerator[expr]/norm]/Expand[Denominator[expr]/norm] /. x_Rational:>FFRatMod[x,FFPrimeNo[0]]


(* ::Section::Closed:: *)
(*Factor*)


(* factor a polynomial *)
Factor[x^10-1]
(* its coefficients are taken modulo a number: *)
Factor[x^10-1,Modulus->7]
(* what's going on here? *)
Factor[11+3x,Modulus->7]//Expand


(* ::Section::Closed:: *)
(*Understanding sparse multiplication and FF degrees*)


list1 = {(x^4+x^2)/(x*y^3),y/x,0,0,0,0};
list2 = {(x+y^2)/(x^3+y^4),0,0,0,0,z^17};
MatrixForm/@{ArrayReshape[list1,{2,3}],ArrayReshape[list2,{3,2}]}
vars = {x,y,z};
FFNewGraph[graph,in,vars];
FFAlgRatFunEval[graph,evalnode1,{in},vars,list1];
FFAlgRatFunEval[graph,evalnode2,{in},vars,list2];
FFAlgMatMul[graph,matnodeA,{evalnode1, evalnode2},2,3,2];
FFGraphOutput[graph,matnodeA];
FFReconstructFunction[graph,vars]

(* or use sparse multiplication: *)
list3 = DeleteCases[list1, 0]
list4 = DeleteCases[list2, 0]
FFAlgRatFunEval[graph,evalnode3,{in},vars,list3];
FFAlgRatFunEval[graph,evalnode4,{in},vars,list4];
FFAlgSparseMatMul[graph,matnodeB,{evalnode3, evalnode4},2,3,2,{{1,2},{}},{{1},{},{2}}];
FFGraphOutput[graph,matnodeB];
FFReconstructFunction[graph,vars]

degs = FFAllDegrees[graph];
FFDumpDegrees[graph,"degfile.m"];


(* FFDumpDegrees should be imported with Import["filename", "Integer64"] *)
(* length of the list is interpreted as: *)
1+1+(1+1+(1+1+1+1)*FFNParsOut[graph,in])*FFNParsOut[graph]
(* which stands for: *)
(* #vars in + output node length + (maxdegnum + maxdegdenom + (maxdegnum in var i + mindegnum in var i + maxdegdenom in var i + maxdegdenom in var i)*#vars in)*output node length *)


degs = Import["degfile.m","Integer64"];
degs
nvars = degs[[1]];
parsout = degs[[2]];
stepsize = 2+4*nvars;
totaldegsnum = Table[degs[[3+stepsize*i]],{i,0,parsout-1}];
totaldegsdenom = Table[degs[[4+stepsize*i]],{i,0,parsout-1}];
maxdegsnum = Association@Table[Rule[vars[[i]], Table[degs[[5+4*(i-1)+stepsize*j]],{j,0,parsout-1}]],{i,nvars}];
maxdegsdenom = Association@Table[Rule[vars[[i]], Table[degs[[7+4*(i-1)+stepsize*j]],{j,0,parsout-1}]],{i,nvars}];


Print["totaldegsnum = ", totaldegsnum]
Print["totaldegsdenom = ", totaldegsdenom]
Print["maxdegsnum = ", maxdegsnum]
Print["maxdegsdenom = ", maxdegsdenom]
Print["Max num total degree = ", Max@totaldegsnum]
Print["Max denom total degree = ", Max@totaldegsdenom]
Print["Max num degree in ", #,": ", Max@maxdegsnum[#]] &/@ vars;
Print["Max denom degree in ", #, ": ", Max@maxdegsdenom[#]] &/@ vars;


(* ::Section::Closed:: *)
(*Generating trace structures*)


(* generate 4! = 24 permutations with a1 fixed at the first position *)
allperms = Prepend[#, 1] &/@ Permutations[{2,3,4,5}];
(* remove 12 of them by applying the reflection identity *)
indepperms = {};
Do[tmp = RotateRight@Reverse@i; If[!MemberQ[indepperms,tmp], AppendTo[indepperms,i]], {i,allperms}]
(* these are also used in the colour run file *)
Print["independent permutations: ", indepperms]


(* ::Section::Closed:: *)
(*Sorting lists by weight*)


ints = {j[t322ZZZMp2314,0,2,0,1,1,1,1,0,0],j[t331ZZMZp1243,2,1,0,0,1,1,1,0,0],j[t331ZZMZp1342,2,1,0,0,1,1,1,0,0]}
fams = ints /. j[t_,x__]:>t
Do[weight[fams[[-ii]]]=ii, {ii,Length@fams}]
SortBy[{j[t322ZZZMp2314,0,2,0,1,1,1,1,0,0],j[t331ZZMZp1243,2,1,0,0,1,1,1,0,0],j[t331ZZMZp1342,2,1,0,0,1,1,1,0,0]}, weight[#/.j[t_,x__]:>t] &]


(* ::Section::Closed:: *)
(*Testing Frules*)


<<"InitTwoLoopToolsFF.m"


(* Frules serve to translate the topo[...] notation onto LiteRed's t***p*** notation *)
(* moreover, we implement by hand the mapping from bubble insertions onto penta/hexa-triangles *)


(* we can check this mappings by making sure that the joint set of {props, ISPs} is consistent *)
(* first, import the ISPs *)


(* ::Subsection::Closed:: *)
(*H3g*)


UserDefinedISPs[topo[{{p1_},{p2_}},{{},{},{}},{{p3_},{p4_}},{{},{},{}},{}]] := Join[{prop[(k[1]+p[p4])^2],prop[(k[2]+p[p1])^2]},{dot[k[1],w[1,p[1],p[2],p[3]]],dot[k[2],w[1,p[1],p[2],p[3]]]}];
UserDefinedISPs[topo[{{p1_},{p2_},{p3_}},{{},{},{}},{{p4_}},{{},{},{}},{}]] := Join[{prop[(k[2]+p[p1])^2],prop[(k[2]+p[p1]+p[p2])^2]},{dot[k[1],w[1,p[1],p[2],p[3]]],dot[k[2],w[1,p[1],p[2],p[3]]]}];
UserDefinedISPs[topo[{{p1_},{p2_},{p3_},{p4_}},{{},{},{}},{},{{},{},{}},{}]] := Join[{prop[(k[2]+p[p1])^2],prop[(k[2]+p[p1]+p[p2])^2],prop[(k[2]+p[p1]+p[p2]+p[p3])^2]},{dot[k[1],w[1,p[1],p[2],p[3]]],dot[k[2],w[1,p[1],p[2],p[3]]]}];
UserDefinedISPs[topo[{{p1_},{p2_}},{{},{},{}},{{p3_}},{{},{},{}},{{p4_}}]] := Join[{prop[(k[1]+p[p4])^2],prop[(k[2]+p[p1])^2]},{dot[k[1],w[1,p[1],p[2],p[3]]],dot[k[2],w[1,p[1],p[2],p[3]]]}];


(* now check the loop momenta in the two families *)
(* remember to impose momentum conservation and a possible FruleShift *)


bubble = topo[{{1}, {4}, {2}, {3}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}];
propsbub = Select[Join[MakePropagators[#],UserDefinedISPs[#]]&@bubble, FreeQ[#,w] &] /. p[4]->-p[1]-p[2]-p[3] // DeleteDuplicates
pentatri = topo[{{1}, {4}, {2}}, {{}, {}, {}}, {{3}}, {{}, {}, {}}, {}];
propspt = Select[Join[MakePropagators[#],UserDefinedISPs[#]]&@pentatri, FreeQ[#,w] &] /. p[4]->-p[1]-p[2]-p[3]  (*/. {k[1] -> k[1] - p[1], k[2] -> k[2] + p[1]} *)
propsbub[[{1, 2, 3, 4, 5, 9, 6, 7, 8}]]===propspt


bubbles = {
		topo[{{4},{1},{2},{3}},{{},{},{}},{},{{},{},{}},{}], (* t511pijkl->t421MZZZpijkl *)
        topo[{{4},{1},{3},{2}},{{},{},{}},{},{{},{},{}},{}], (* {1, 2, 3, 4, 5, 9, 6, 7, 8} *)
        topo[{{4},{2},{1},{3}},{{},{},{}},{},{{},{},{}},{}],
        topo[{{4},{2},{3},{1}},{{},{},{}},{},{{},{},{}},{}],
        topo[{{4},{3},{1},{2}},{{},{},{}},{},{{},{},{}},{}],
        topo[{{4},{3},{2},{1}},{{},{},{}},{},{{},{},{}},{}],
       
        topo[{{1},{4},{2},{3}},{{},{},{}},{},{{},{},{}},{}], (* t511pijkl->t421ZMZZpijkl *)
        topo[{{3},{4},{1},{2}},{{},{},{}},{},{{},{},{}},{}], (* {1, 2, 3, 4, 5, 9, 6, 7, 8} *)
        topo[{{3},{4},{2},{1}},{{},{},{}},{},{{},{},{}},{}],
       
        topo[{{1},{4},{3},{2}},{{},{},{}},{},{{},{},{}},{}], (* t511pijkl->t421MZZZpjkli *)
        topo[{{2},{4},{1},{3}},{{},{},{}},{},{{},{},{}},{}], (* {2, 3, 4, 1, 7, 5, 6, 8, 9} *)
        topo[{{2},{4},{3},{1}},{{},{},{}},{},{{},{},{}},{}]
       };


Do[Print[
"F["<>ToString[tt]<>", p__] :> F["<>(tt/. topo[x__]:>"t421MZZZp"<>StringJoin[ToString/@Flatten@DeleteCases[x,{}]])<>", Sequence @@ {p}[[{1, 2, 3, 4, 5, 9, 6, 7, 8}]]],"
],{tt,bubbles[[;;6]]}]

Do[Print[
"F["<>ToString[tt]<>", p__] :> F["<>(tt/. topo[x__]:>"t421ZMZZp"<>StringJoin[ToString/@Flatten@DeleteCases[x,{}]])<>", Sequence @@ {p}[[{1, 2, 3, 4, 5, 9, 6, 7, 8}]]],"
],{tt,bubbles[[7;;9]]}]

Do[Print[ (* note the RotateLeft here *)
"F["<>ToString[tt]<>", p__] :> F["<>(tt/. topo[x__]:>"t421MZZZp"<>StringJoin[ToString/@RotateLeft[Flatten@DeleteCases[x,{}]]])<>", Sequence @@ {p}[[{2, 3, 4, 1, 7, 5, 6, 8, 9}]]],"
],{tt,bubbles[[10;;]]}]

(*  FruleShift[topo[{{1}, {4}, {3}, {2}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}]] = {k[1] -> k[1] - p[1], k[2] -> k[2] + p[1]};
    FruleShift[topo[{{2}, {4}, {1}, {3}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}]] = {k[1] -> k[1] - p[2], k[2] -> k[2] + p[2]};
    FruleShift[topo[{{2}, {4}, {3}, {1}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}]] = {k[1] -> k[1] - p[2], k[2] -> k[2] + p[2]}; *)


(* ::Subsection::Closed:: *)
(*H4g*)


UserDefinedISPs[topo[{{p1_},{p2_},{p3_}},{{},{},{}},{{p4_},{p5_}},{{},{},{}},{}]] := {prop[(k[1]+p[p5])^2],prop[(k[2]+p[p1])^2],prop[(k[2]+p[p1]+p[p2])^2]};
UserDefinedISPs[topo[{{p1_},{p2_},{p3_},{p4_}},{{},{},{}},{{p5_}},{{},{},{}},{}]] := {prop[(k[2]+p[p1])^2],prop[(k[2]+p[p1]+p[p2])^2],prop[(k[2]+p[p1]+p[p2]+p[p3])^2]};
UserDefinedISPs[topo[{{p1_},{p2_},{p3_},{p4_},{p5_}},{{},{},{}},{},{{},{},{}},{}]] := {prop[(k[2]+p[p1])^2],prop[(k[2]+p[p1]+p[p2])^2],prop[(k[2]+p[p1]+p[p2]+p[p3])^2],prop[(k[2]-p[p5])^2]};

UserDefinedISPs[topo[{{p1_},{p2_},{p3_}},{{},{},{}},{{p4_}},{{},{},{}},{{p5_}}]] := {prop[(k[1]+k[2]-p[p1])^2],prop[(k[1]+k[2]-p[p1]-p[p2])^2],prop[(k[1]+k[2]-p[p1]-p[p2]-p[p3])^2]};
UserDefinedISPs[topo[{{p1_},{p2_}},{{},{},{}},{{p4_},{p5_}},{{},{},{}},{{p3_}}]] := {prop[(k[1]+p[p5])^2],prop[(k[2]+p[p1])^2],prop[(k[2]+p[p1]+p[p2])^2]};


UserDefinedISPs[topo[{{p1_},{p2_},{p3_}},{{},{},{},{}},{{p4_},{p5_}}]] := {prop[(k[2]+p[p1])^2],prop[(k[2]+p[p1]+p[p2])^2],prop[(k[1]+p[p5])^2],prop[(k[1]+k[2])^2]};
UserDefinedISPs[topo[{{p1_},{p2_},{p3_},{p4_}},{{},{},{},{}},{{p5_}}]] := {prop[(k[2]+p[p1])^2],prop[(k[2]+p[p1]+p[p2])^2],prop[(k[2]+p[p1]+p[p2]+p[p3])^2],prop[(k[1]+k[2])^2]};


bubbles = {topo[{{5}, {1}, {2}, {3}, {4}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}], (* t611pijklm->t521MZZZZpijklm *)
        topo[{{5}, {1}, {4}, {3}, {2}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],    (* {1, 2, 3, 4, 5, 6, 11, 7 ,8 ,9, 10} *)
        topo[{{5}, {2}, {1}, {4}, {3}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],
        topo[{{5}, {2}, {3}, {4}, {1}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],
        topo[{{5}, {3}, {2}, {1}, {4}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],
        topo[{{5}, {3}, {4}, {1}, {2}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],
        topo[{{5}, {4}, {1}, {2}, {3}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],
        topo[{{5}, {4}, {3}, {2}, {1}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],

        topo[{{1}, {5}, {2}, {3}, {4}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],   (* t611pijklm->t521MZZZZpijklm *)
        topo[{{1}, {5}, {4}, {3}, {2}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],   (* {1, 2, 3, 4, 5, 6, 11, 7 ,8 ,9, 10} *)
        topo[{{2}, {5}, {1}, {4}, {3}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],
        topo[{{2}, {5}, {3}, {4}, {1}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],
        topo[{{3}, {5}, {2}, {1}, {4}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],
        topo[{{3}, {5}, {4}, {1}, {2}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],
        topo[{{4}, {5}, {1}, {2}, {3}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],
        topo[{{4}, {5}, {3}, {2}, {1}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],

        topo[{{1}, {2}, {5}, {3}, {4}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],  (* t611pijklm->t521ZMZZZpjklmi *)
        topo[{{2}, {3}, {5}, {4}, {1}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],  (* {2, 3, 4, 5, 1, 8, 6, 7, 9, 10, 11} *)
        topo[{{3}, {2}, {5}, {1}, {4}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}],
        topo[{{3}, {4}, {5}, {1}, {2}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}]};


bubble = bubbles[[20]]
propsbub = Select[Join[MakePropagators[#],UserDefinedISPs[#]]&@bubble, FreeQ[#,w] &] /. p[5]->-p[1]-p[2]-p[3]-p[4] // DeleteDuplicates
hexatri = topo[{{4}, {5}, {1}, {2}}, {{}, {}, {}}, {{3}}, {{}, {}, {}}, {}];
propsht = Select[Join[MakePropagators[#],UserDefinedISPs[#]]&@hexatri, FreeQ[#,w] &] /. p[5]->-p[1]-p[2]-p[3]-p[4]  /. {k[1] -> k[1] - p[3], k[2] -> k[2] + p[3]} 
propsbub[[{2, 3, 4, 5, 1, 8, 6, 7, 9, 10, 11}]]===propsht


Do[Print[
"F["<>ToString[tt]<>", p__] :> F["<>(tt/. topo[x__]:>"t521MZZZZp"<>StringJoin[ToString/@Flatten@DeleteCases[x,{}]])<>", Sequence @@ {p}[[{1, 2, 3, 4, 5, 6, 11, 7 ,8 ,9, 10}]]],"
],{tt,bubbles[[;;8]]}]

Do[Print[
"F["<>ToString[tt]<>", p__] :> F["<>(tt/. topo[x__]:>"t521ZMZZZp"<>StringJoin[ToString/@Flatten@DeleteCases[x,{}]])<>", Sequence @@ {p}[[{1, 2, 3, 4, 5, 6, 11, 7 ,8 ,9, 10}]]],"
],{tt,bubbles[[9;;16]]}]

Do[Print[
"F["<>ToString[tt]<>", p__] :> F["<>(tt/. topo[x__]:>"t521ZMZZZp"<>StringJoin[ToString/@RotateLeft[Flatten@DeleteCases[x,{}]]])<>", Sequence @@ {p}[[{2, 3, 4, 5, 1, 8, 6, 7, 9, 10, 11}]]],"
],{tt,bubbles[[17;;20]]}]

    FruleShift[topo[{{1}, {2}, {5}, {3}, {4}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}]] = {k[1] -> k[1] - p[1], k[2] -> k[2] + p[1]};
    FruleShift[topo[{{2}, {3}, {5}, {4}, {1}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}]] = {k[1] -> k[1] - p[2], k[2] -> k[2] + p[2]};
    FruleShift[topo[{{3}, {2}, {5}, {1}, {4}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}]] = {k[1] -> k[1] - p[3], k[2] -> k[2] + p[3]};
    FruleShift[topo[{{3}, {4}, {5}, {1}, {2}}, {{}, {}, {}}, {}, {{}, {}, {}}, {}]] = {k[1] -> k[1] - p[3], k[2] -> k[2] + p[3]};


(* ::Section::Closed:: *)
(*Drawing ints in j[...] notation*)


MomMode = PSanalytic;
<<"InitTwoLoopToolsFF.m";
Get["setupfiles/Setup_phi3g.m"];


UserDefinedISPs[topo[{{p1_},{p2_}},{{},{},{}},{{p3_},{p4_}},{{},{},{}},{}]] := Join[{prop[(k[1]+p[p4])^2],prop[(k[2]+p[p1])^2]},{dot[k[1],w[1,p[1],p[2],p[3]]],dot[k[2],w[1,p[1],p[2],p[3]]]}];


topo2name = {
 topo[{{1}, {2}}, {{}, {}, {}}, {{3}, {4}}, {{}, {}, {}}, {}]->t331ZZZMp1234
};
name2topo = topo2name /. Rule[a_,b_]:>Rule[b,a];


jToINT[int_j]:=Module[{pows,tp,tp2,props,pows1,pows2,num,res},
tp=int[[1]];
tp2=tp /. name2topo;
props=DeleteDuplicates[MakePropagators[tp2]];
pows=int[[2;;]] /. j->List;
pows1=pows[[1;;Length[props]]];
pows2=pows[[Length[props]+1;;]];
num=Times@@(Select[UserDefinedISPs[tp2],FreeQ[#,w]&]^(-pows2));
res=INT[ToDotProducts[num],pows1,tp2];
Return[res];
];


j[t331ZZZMp1234,1,1,1,1,1,1,0,0,0]/. int_j:>jToINT[int] // SquashTopology // SymbolForm


(* ::Section::Closed:: *)
(*Spurious zeroes in when processing DiagramNumerators*)


Quit[];


<<"InitTwoLoopToolsFF.m"
<<"setupfiles/Setup_H3g.m"


(* ::Text:: *)
(*Sometimes, SortLoopNumerator does not recognise that the coefficient of a paritcular loop monomial is zero.*)
(*This happens when the coefficient is a really complicated cancellation that hasn't been simplified to 0.*)


(*monomialbasis = <<"helicityamplitudes/H3g/monomialbasis.m";
mockrules = <<"helicityamplitudes/H3g/mockrules.m";*)
num = <<"helicityamplitudes/H3g/num.m";


{loopvars, monomialbasis} = GetVarsAndMonoms[num[[11]]]
mockrules =Table[Rule[monomialbasis[[i]],mon[i]], {i,Length[monomialbasis]}];
FreeQ[num[[6;;8]] /. mockrules, k]


num[[8]]+num[[11]]


(* ::Text:: *)
(*We can check if the coefficient is indeed spurious. We multiply each loop monomial by some random flag and then evaluate the coefficient numerically to see if it's really zero.*)
(*If it is, we use GetVarsAndMonoms to extract that loop monomial.*)
(*The use of the flag is necessary, because otherwise GetVarsAndMonoms will just simplify this term and say there were no loop monomials to begin with.*)
(*But we want the output to be an extra rule that we can apply to num in SLN.*)


CheckSpuriousMonoms[expr_]:= Module[{tmp,vars,rules,spurious},
tmp = Select[expr, !FreeQ[#, k] &];
vars = Select[Variables[tmp], !FreeQ[#,k] &];
rules = RuleDelayed[#,c[RandomInteger[{9,99999}]]*#]&/@vars;
Print["rules = ", rules];
tmp = tmp /. rules /. (f:(spA|spB|spAB|s|dot))[x__]:>GetMomentumTwistorExpression[f[x],PSnumeric];
Print[Simplify[tmp/.c[_]:>1]];
If[Simplify[tmp/.c[_]:>1]=!=0, Print["Could not substitute all monomials into mock variables"](*;Quit[]*)];
spurious = GetVarsAndMonoms[tmp][[2]];
Return[Thread@Rule[spurious, 0]]
];


CheckSpuriousMonoms[num[[8]]+num[[11]]]


(* ::Section::Closed:: *)
(*Simplifying with assumptions*)


Simplify[Log[-m2]+Log[-s]+Log[-t], s<0&&t<0]


(* ::Section::Closed:: *)
(*Selecting simple numerical sij's*)


<<"InitTwoLoopToolsFF.m"
<<"setupfiles/Setup_phi4g.m"	
<<"setupfiles/Setup_5pt_extra.m"


bits = 10^9;
Do[
expr = GetMomentumTwistorExpression[{s[1,2],s[2,3],s[3,4],s[4,5],s[1,5],s[5]},PSanalytic] /. ex[1]->1 /. Table[ex[ii]->RandomPrime[{9,99}]/RandomPrime[{9,99}], {ii,2,6}];
newbits = Total[BitLength[Numerator[#]]+BitLength[Denominator[#]]&/@expr];
If[newbits<bits&&Cases[expr[[2;;]], 1]=={}&&Cases[expr[[2;;]], 0]=={}&&AllTrue[expr, #>0 &], Print[expr]; bits=newbits];
,{jj,100}]


(* ::Section::Closed:: *)
(*trp, trm, tr5*)


<<"InitTwoLoopToolsFF.m"
<<"setupfiles/Setup_phi4g.m"	
<<"setupfiles/Setup_5pt_extra.m"


(* p1p2 = (|1>[1| + |1]<1|)*((|2>[2| + |2]<2|) = |1>[12]<2| + |1]<12>[2| *)
(* p1p2p3p4 = |1>[12]<23>[34]<4| + |1]<12>[23]<34>[4| *)
(* P+p1p2p3p4 = |1>[12]<23>[34]<4| *)
(* P-p1p2p3p4 = |1]<12>[23]<34>[4| *)
(* tr(P+,p1p2p3p4) = [12]<23>[34]<41> *)
(* tr(P-,p1p2p3p4) = <12>[23]<34>[41] *)
(* helicity projectors: P+ = (1+g5)/2, P- = (1-g5)/2, so that g5 = P+ - P- *)
(* tr(g5,p1p2p3p4) = tr(P+,p1p2p3p4) - tr(P-,p1p2p3p4) = [12]<23>[34]<41> - <12>[23]<34>[41] *)
GetMomentumTwistorExpression[trp[1,2,3,4] - spB[1,2]*spA[2,3]*spB[3,4]*spA[4,1], PSanalytic]
GetMomentumTwistorExpression[trm[1,2,3,4] - spA[1,2]*spB[2,3]*spA[3,4]*spB[4,1], PSanalytic]
GetMomentumTwistorExpression[tr5[1,2,3,4] - (trp[1,2,3,4]-trm[1,2,3,4]), PSanalytic]
GetMomentumTwistorExpression[trp[1,2,3,4] - (-s[1,2]*s[2,3]), PSanalytic]
GetMomentumTwistorExpression[trm[1,2,3,4] - (-s[1,2]*s[2,3]), PSanalytic]


GetMomentumTwistorExpression[spAB[2,5,2]-2*dot[p[2],p[5]],PSanalytic]
GetMomentumTwistorExpression[spAB[2,5,2]-(s[2,5]-s[5]),PSanalytic]


<<"InitTwoLoopToolsFF.m"
<<"setupfiles/Setup_4g.m"	
<<"setupfiles/Setup_4pt_extra.m"


(* https://arxiv.org/pdf/hep-th/9611127.pdf *)
GetMomentumTwistorExpression[trm[1,2,3,4]/(spA[1,2]*spA[2,3]*spA[3,4]*spA[4,1]),PSanalytic]
(* https://arxiv.org/pdf/0810.2964.pdf *)
GetMomentumTwistorExpression[spB[1,2]*spB[3,4]/(spA[1,2]*spA[3,4]), PSanalytic]
GetMomentumTwistorExpression[-s[1,2]*s[2,3]/(spA[1,2]*spA[2,3]*spA[3,4]*spA[4,1]), PSanalytic]
GetMomentumTwistorExpression[-s[1,2]*s[2,3]-trm[1,2,3,4], PSanalytic]


GetMomentumTwistorExpression[spAB[1,2,3]*spAB[3,2,1] - trm[1,2,3,2], PSanalytic]
GetMomentumTwistorExpression[trm[1,2,3,2] - (2*spAB[1,2,1]*dot[p[2],p[3]]-s[1,3]*dot[p[2],p[2]]), PSanalytic]
(* for massive particles, we cannot make the spAB <-> trm swap, but the RHS still holds *)
GetMomentumTwistorExpression[spAB[1,5,3]*spAB[3,5,1] - (2*spAB[1,5,1]*dot[p[5],p[3]]-s[1,3]*dot[p[5],p[5]]), PSanalytic]


<<"InitTwoLoopToolsFF.m"
<<"setupfiles/Setup_phi3g.m"
<<"setupfiles/Setup_4pt_extra.m"


(*vars = Variables[A2cc /. {MyLog[x_,y_]->Log[x/y], MyDiLog[s_,t_]->PolyLog[2,1-s/t]}][[;;-5]] /. {mu->1, delta->0};*)
vars = {Log[-(1/s[4])],Log[s[1,2]/s[4]],Log[s[2,3]/s[4]],Log[s[1,3]/s[4]],Log[s[2,4]/s[1,4]],Log[s[1,4]/s[3,4]],Log[s[3,4]/s[2,4]],PolyLog[2,1-s[4]/s[1,4]],PolyLog[2,1-s[4]/s[2,4]],PolyLog[2,1-s[4]/s[3,4]]}
(*vars = Insert[vars,{Log[s[1,2]/s[4]],Log[s[2,3]/s[4]],Log[s[1,3]/s[4]]},2] // Flatten*)
varslogs4 = Simplify@GetMomentumTwistorExpression[vars[[1]], PSanalytic] /. Log[x_]->x;
varslogs = Simplify@GetMomentumTwistorExpression[vars[[2;;7]], PSanalytic] /. Log[x_]->x;
varsdilogs = Simplify@GetMomentumTwistorExpression[vars[[8;;10]], PSanalytic] /. PolyLog[2,x_]->x;


Reduce[varslogs4>0&&varslogs>0&&varsdilogs<1,{ex[1],ex[2],ex[5]}]
newpoint = {ex[1]->-1, ex[2]->3, ex[5]->6}


vars
GetMomentumTwistorExpression[vars,PSanalytic] /. newpoint // N


invariants = {s12, s15, s23, s34, s45};
rules = {p1->p5, p2->p1, p3->p2, p4->p3, p5->p4} /. Rule[a_,b_]:>Rule[ToExpression[StringDrop[ToString[a],1]],ToExpression[StringDrop[ToString[b],1]]]
digits = IntegerDigits[ToExpression[StringDrop[ToString/@invariants, 1]]] /. rules /. x_Integer:>ToString[x]
ToExpression/@(StringJoin["s", #]&/@StringJoin/@Sort/@%)


(* ::Section:: *)
(*DEs, alphabet and symbol*)


(* ::Subsection::Closed:: *)
(*Finding the alphabet*)


(* ::Text:: *)
(*First, obtain the differential equation matrices for all relevant families using, e.g. packages/finiteflow-mathtools/examples/differential_equations/differential_equations.wl*)
(*For each family, there will be one matrix for each invariant (here: {s12, s23, s4}).*)
(*Then, load all the DE matrices, Flatten and Together to bring everything under a common denominator.*)
(*Select denominator factors (FactorList) - they will correspond to the letters of the alphabet.*)
(*If there are no square roots, this will suffice to get the alphabet, otherwise we have to be more careful.*)
(*letters = {s12,s23,s4,s12+s23,s12-s4,s23-s4,s12+s23-s4,-s23+s4,-s12+s4} - 9 letters selected, but we can see the last two are superfluous.*)
(*The letters need to be linearly independent:*)
(*Sum[c[i]*dLog[w[i]],  {i, Length[letters]}]= 0   <=>  c[i]=0  for all c[i].*)
(*Now, df[x,y] = D[f[x,y],x] dx + D[f[x,y],y] dy   , where d is the differential, D is the normal derivative. In this case:*)
(*dLog[w(s12,s23,s4)] = 1/w*dw(s12,s23,s4) = 1/w*(Dw/d(s12) + Dw/d(s23) + Dw/d(s4)).*)
(*This gives 3 separate equations, one for each invariant.*)
(*It is not enough to solve the system of these 3 equations just as it is, because it will consider c[i] as functions of {s12, s23, s4}, whereas c[i] can have any numerical value.*)
(*We therefore need to sample the system on some numerical points. Here, 5 samples is enough.*)


letters = {s12,s23,s4,s12+s23,s12-s4,s23-s4,s12+s23-s4,-s23+s4,-s12+s4};
(*letters = letters[[;;7]]*)
system = {
Sum[c[ii]/letters[[ii]]*D[letters[[ii]], s12], {ii,Length[letters]}]==0, 
Sum[c[ii]/letters[[ii]]*D[letters[[ii]], s23], {ii,Length[letters]}]==0,
Sum[c[ii]/letters[[ii]]*D[letters[[ii]], s4], {ii,Length[letters]}]==0
};

Solve[
Flatten@Table[system /. Thread[{s12,s23,s4}->GetMomentumTwistorExpression[{s[1,2],s[2,3],s[4]}, PSanalytic] /. Thread[{ex[1],ex[2],ex[5]}:>RandomInteger[{9,999}]/1000]], {ii,5}],
c/@Range[Length[letters]]]


(* ::Text:: *)
(*We can see that the last two letters are indeed eliminated. If we repeat the above Solve with letters = letters[[;;7]], no solutions will be found, which means that the first 7 letters are independent.*)


letters = letters[[;;7]]


(* ::Text:: *)
(*In general, however, this alphabet may not be enough, because we obtained it by looking at DEs of families in their 'main' permutations, e.g. 1234 for ZZZM, 1423 for ZMZZ, etc.*)
(*The alphabet needs to be closed under all possible permutations of external legs.*)
(*In this case, it is already closed, but let's see how we would extend if it wasn't.*)
(*We need to permute the massless {1,2,3} legs and keep the massive leg {4} fixed.*)
(*We generate all permutations of the 7 letters above and again check for their linear independence.*)


rules = Thread@Rule[{1,2,3},#]&/@Permutations[{1,2,3}];
newletters = (letters /. {s12->s[1,2],s23->s[2,3],s4->s[4]} /. #)&/@ rules /. s[x_,y_]:>s[y,x] /; x>y /. s[1,3]->s[4]-s[1,2]-s[2,3] // Flatten // DeleteDuplicates;
newletters = newletters /. {s[1,2]->s12,s[2,3]->s23,s[4]->s4}

(*newletters = newletters[[;;7]]*)

system = {
Sum[c[ii]/newletters[[ii]]*D[newletters[[ii]], s12], {ii,Length[newletters]}]==0, 
Sum[c[ii]/newletters[[ii]]*D[newletters[[ii]], s23], {ii,Length[newletters]}]==0,
Sum[c[ii]/newletters[[ii]]*D[newletters[[ii]], s4], {ii,Length[newletters]}]==0
} 

Solve[
Flatten@Table[system /. Thread[{s12,s23,s4}->GetMomentumTwistorExpression[{s[1,2],s[2,3],s[4]}, PSanalytic] /. Thread[{ex[1],ex[2],ex[5]}:>RandomInteger[{9,999}]/1000]], {ii,5}],
c/@Range[Length[newletters]]]


(* ::Text:: *)
(*The {1,2,3} permutations generate 13 letters, but again, only 7 of them are independent and they are the same as the original ones:*)


newletters = newletters[[;;7]];
ContainsExactly[newletters, letters]


(* ::Text:: *)
(*However, in a general case, this wouldn't be true and the permutation procedure would increase the number of letters.*)


newletters


(* ::Subsection:: *)
(*Finding Atilde*)


(* ::Text:: *)
(*Say we have only three letters, dependent on 2 invariants:*)


letters = {s12,s23,s12-s23};


(* ::Text:: *)
(*We're trying to find Atilde, but we don't know what a[i] are:*)


(*Atilde = Sum[a[i]*Log[letters[[i]]], {i, Length[letters]}]*)
(*As12 = Sum[a[i]*D[Log[letters[[i]]], s12], {i, Length[letters]}]*)
(*As23 = Sum[a[i]*D[Log[letters[[i]]], s23], {i, Length[letters]}]*)


(* ::Text:: *)
(*These matrices wrt to each invariant can be obtained by first differentiating and then IBP-reducing the resulting expressions back onto the same set of MIs.*)
(*In practice, this is done in differential_equations.wl with a combination of LiteRed and FiniteFlow.*)


As12 = {{1/s12-1/(s12-s23),3/(s12-s23)},{-(2/s12),3/s12-7/(s12-s23)}};
As23 = {{1/(s12-s23)+4/s23,-(3/(s12-s23))+5/s23},{3/s23,7/(s12-s23)-3/s23}};


(* ::Text:: *)
(*Now we're going to perform the linear fit by numerically sampling equations for each entry of a[k] separately.*)


sols = {};
Do[
(* system for each entry of a[k] separately *)
system = 
{
As12[[i,j]] == Sum[a[k,i,j]*D[Log[letters[[k]]], s12], {k, Length[letters]}],
As23[[i,j]] == Sum[a[k,i,j]*D[Log[letters[[k]]], s23], {k, Length[letters]}]
};

(* solve each such system by numerically sampling both sides at random points *)
(* may have to experiment with how many times we need to sample the system *)
sol = Solve[Flatten@Table[system /. Thread[{s12,s23}->GetMomentumTwistorExpression[{s[1,2],s[2,3]}, PSanalytic] /. Thread[{ex[1],ex[2],ex[5]}:>RandomInteger[{9,999}]/1000]], {ii,Length[letters]}]];
AppendTo[sols, sol[[1]]];
,{i, Length[As23]}, {j,Length[As23]}]


(* put the matrices back together *)
sols;
mat = Table[a[k,i,j], {k, Length[letters]}, {i, Length[As12]}, {j,Length[As12]}] /. Flatten[sols];
(* construct Atilde *)
Atilde = Sum[mat[[i]]*Log[letters[[i]]], {i, Length[letters]}];


(* ::Text:: *)
(*Perform sanity checks:*)


DeleteDuplicates@Flatten@Simplify[D[Atilde, s12]-As12]=={0}
DeleteDuplicates@Flatten@Simplify[D[Atilde, s23]-As23]=={0}


D[PolyLog[2,x], x] - D[Log[x],x]*PolyLog[1,x]


D[PolyLog[2,x/y], x]+D[PolyLog[2,x/y], y] - (D[Log[x],x]-D[Log[y],y])*PolyLog[1,x/y] // Simplify


D[PolyLog[2,x+y^2], x]+D[PolyLog[2,x+y^2], y]
D[Log[x+y^2],x]*PolyLog[1,x+y^2]


Integrate[1/(x+y^2), x]
Integrate[2y/(x+y^2), y]


(* ::Subsection:: *)
(*Table with a variable number of iterators*)


tab[n_] := Table[f[a,b,...], {a, 3}, {b, 3}, ...]
f[##] &@@ ii/@Range[5]
{##, 3} &/@ ii/@Range[5]
tab[n_] := Table[f[##] &@@ ii/@Range[5], {##, 3} &/@ ii/@Range[5]]
tab[1]





(* ::Section:: *)
(*Selecting numerical sij values*)


<<"InitTwoLoopToolsFF.m"
<<"setupfiles/Setup_phi4g.m"


RationalQ[x_]:=(Head[x]===Rational||IntegerQ[x]);


Do[
rnd=Join[RandomInteger[{9,99}, 5], {RandomPrime[{9,99}]}];
eqs = Thread[rnd- GetMomentumTwistorExpression[{s[1,2],s[2,3],s[3,4],s[4,5],s[1,5],s[5]},PSanalytic] == 0];
sol = Solve[eqs,ex/@Range[6]];
(* select sij such that the corresponding ex[i] are real and positive *)
If[AllTrue[Flatten[sol][[;;,2]], RationalQ[#]&&#>0 & ], Print[rnd]]
,{ii,1000}]




# Comments on permutations

In the computation of scattering amplitudes, we most often use master integral bases which satisfy 'uniform transcendentality'. These bases are supplied for only one permutation of a given family, so to calculate an amplitude, in general we need to obtain the permuted bases for all the required permutations. This involves permuting the sij's, including any potential square roots and tr5. The latter is a pseudo-scalar invariant:

$$tr_5 = 4i\epsilon_{\mu\nu\rho\sigma} p_1^\mu p_2^\nu p_3^\rho p_4^\sigma = [12]<23>[34]<41> - <12>[23]<34>[41]$$

whose sign might change under a given permutation. As we will see below, this can lead to some practical complications.

We have two main approaches for dealing with IBPs and the reduction of an amplitude onto master integrals:
1. The **non-alt approach** - the IBP relations are generated for all the family permutations required.
2. The **alt approach** - the IBP relations are generated only for the permutation which defines a given family, e.g. the *mzz pentagon-box*, or *t431MZZZ* in our notation. The IBPs are solved only for this permutation and to get the solution in other permutations, the sij's are permuted *numerically* within a `FiniteFlow` graph. The master integrals for these other permutations are obtained simply by relabelling the family name from the main permutation to the target one (this point will become important later).

In general, the *alt* approach is preferred as it involves smaller sizes of the IBP system, quicker solution, etc. 

## Expressing the UT basis in our notation

The UT bases provided by `2005.04195` and `2107.14180` are in kinematics where q1^2 != 0. We need to translate them onto our kinematics, where p5^2 != 0. In addition, we also want these results to be compatible with the `PentagonFunctions++` package of `2110.10111`. This package always evaluates the pentagon functions in the `45->123` channel. Combining these two requirements, we decide to use the map: `swap15 = {q1->p5, q5->p1, q2->p4, q4->p2, q3->p3}` (yes, I know, stupid name), which ensures that p5 is massive and that a process we compute in the `12->345` channel can be evaluated using `PentagonFunctions++`.

In practice, this amounts to permuting the sij's according to `swap15` and mapping the integrals as defined in the papers onto the integrals as defined in our notation. This can be done by finding `LiteRed`'s jExtRules from former to latter. At this point, we have UT master integral bases written in our kinematics and according to our definition of propagators.

Note that this step does not involve swapping the sign of tr5. Schematically, we can see this from:
```
tr5 = 4*I*eps(1,2,3,4)

tr5(q) /. swap15 = 4*I*eps(5,4,3,2) = 
-4*I*eps(1,4,3,2) = +4*I*eps(1,2,3,4) = tr5(p)
```
where we used momentum conservation in the third equality. The proper way to do this is to explicitly compute tr5(q) in terms of the 4-momentum components, apply `swap15` and verify that it equals the explicit expression for tr5(p).

## Permuting the UT basis definitions into all needed permutations

Now it's time to obtain the UT bases for other families needed in the amplitude. Why do we need that information in the first place? In the *non-alt* approach, we need to generate the IBPs in all permutations, so we need to have the corresponding integral bases. For the *alt* approach, it might seem that this is not needed as the IBPs are generated only in one permutation per family. This is not true, however, for the following reason.

In general, MIs for one permutation might overlap with MIs for other permutations (e.g, a bubble in s123=s45 is the same as bubble in s312=s54). In order to reduce this degeneracy of MIs, we need to find relations between them. Then, the set of MIs for all permutations considered together is a subset of their union. Note that this procedure of finding the **mappings** between MIs is not necessary in the *non-alt* approach, because there the IBP relations are generated for all the perms together, so any mappings between MIs are automatically found.

### Note on the mappings between MIs
Working out these mappings in the *alt* approach is not strictly necessary, depending on what we want to do. It will reduce the number of coefficients that we work with, but is not necessary for the result to be correct. Indeed, we can compute the amplitude without using the mappings. If we stop at the level of MIs, they will not be linearly independent (i.e. they will not be the *true* MIs). This is fine if we want to, e.g. evaluate such a result using `AMFlow`, but we have to remember that if we're trying to perform the 'gauge check', such a result in terms of MIs will not give us 0, because there are missing relations between the MIs that have not been implemented using the mappings. On the other hand, if we expand the amplitude onto pentagon functions, we don't need to do anything extra, because pentagon functions are closed under all permutations, so any relations that have not been imposed at the level of MIs, will now be automatically imposed at the level of pentagon functions. In general, it is a good idea to work out the mappings, as this can drastically reduce the number of coefficients processed in the `FF` graph. For the mappings, we need to generate the IBPs in all permutations, but we can usually do it for a rank lower than that needed in the amplitude. We can also do this task numerically, because mappings between UT integrals do not have any s_ij or eps dependence, only simple rational numbers (the most complicated I've ever seen was 57/2). If that's not the case, we must have messed up either the procedure of determining the mappings, or the MIs we're working with are not actually UT (do check this with differential equations).

### Figuring out the correct permutation

Once we've understood that, regardless of the *alt* or *non-alt* approach, we need to generate the basis definitions for all family permutations present in the amplitudes, we have to carefully work out the correct way to do this. We want to permute the family which was obtained in the 'Expressing the UT basis in our notation' section (henceforth referred to as **main topology** or **main family**) onto the target topology.

One easy mistake we have to avoid is that the main topology permutation for massive processes is not just p12345, as it might have been for massless processes. This is because, if we have a massive external leg (in our case, p5), we have to keep track of its position within the family. For example, we don't have just one pentagon-box family anymore - we have 3, depending on the position of p5. To make this more concrete: if we want to obtain the permuted MI basis for p53412, and the main family is p53214, we need to permute:

```p53214 -> p53412```

Note that at one point we used to define the main families according to the massive leg inserted between the massless ones in the 1234 permutation, for example MZZZZ was p51234, ZMZZZ was p15234 and ZZZZM was p12345. However, because we need to use `swap15` to uniformly map our 12 channel onto the 45 channel of `PentagonFunctions++`, our 'main families' will be more complicated than that and can be worked out by applying `swap15` to diagrams in Fig. 1 of `2005.04195` and Fig. 1 of `2107.14180` (*remember to match the position of the propagators*). In particular, note that e.g. t422MZZZZ and t431MZZZZ might not have the same main family permutation. 

Once we've figured out the correct permutation, we can permute the labels of the j[...] integrals.

***Don't forget to also permute the sij's and tr5!*** (more on this below)

## UT vs root-free basis
At 5pt, we have square roots and tr5 included in the UT basis definitions. We want to get rid of them, because `FF` cannot handle roots. We therefore prepare a root-free/semi-UT/semi-canonical basis which is obtained by stripping the full basis off all such square roots and tr5. 

For the explicit square roots, we just normalise them away and save these prefactors in a separate file. Then, when expanding the root-free basis onto special functions, we need to remember to divide the expansion by the corresponding normalisation.

**Important note:** make sure to first permute the UT basis and *then* strip it off square roots. If we do it the other way around, we'll obtain permuted MIs normalised by incorrect, unpermuted square roots. (*NB: the number of square roots will not explode after considering all permutations - the three explicit square roots map onto each other under permutations of massless legs.*)

We also need to handle tr5 and here the situation is more nuanced. We could go about it in several ways:
1. Explicitly include tr5 in the sij list: this is a bad idea - it increases the complexity of the `FF` graph and of course `FF` does not know that in fact tr5 is a function of MTs.
2. Evaluate tr5 in terms of MTs and absorb it into the coefficients - better, but we don't want to separate it from the pentagon functions, because then we need to swap the sign all parity-odd pentagon functions when evaluating them in `PentagonFunctions++`. Note that this is only possible because MTs rationalise tr5.
3. Keep tr5 and F[i,j] together, so that in the expansion to pentagon functions, mono[...] is parity-even and we don't ever have to swap the sign of F during evaluation.

In the 5pt 1-mass computations, we choose approach 3.

### GMN_alt the sign of tr5

The fact that the MIs that are used by GMN_alt are square-root free has some serious implications. GMN_alt does not know about tr5 - it work with parity-even scalar master integrals. Therefore, when it permutes the solution of the IBP system from mi[perm] to mi[sigma(perm)], **this is simply a relabelling, rather than an actual permutation**. Strictly speaking, we should add a bit of code in GMN_alt to make this relabelling mi[perm] -> -mi[sigma(perm)], where sigma - an odd permutation from main topo onto the target. An alternative approach - which is what we employ - is to deliberately skip the tr5 permutation when generating the UT bases in all needed permutations.

Note: this problem is not present in the `non-alt` approach, where the IBPs are generated for all the required permutations together and there is no further trickery going on.

***this fixes the GMN_alt reduction onto MIs, which can later be evaluated using AMFlow through their expressions in terms of j[...].***
***describe the problem with expansion to pentagon functions - we should, in fact permute the tr5 when generating the permuted UT basis, but this sign should go into the prefactors that we include together with the F[i, j] in mono[...], NOT in the definition of the MIs like mi[fam, i]->...***

### Aligning the loop-momenta - the MU[i, j] issue
As explained in the section 'Expressing the UT basis in our notation', we employ `swap15` to translate the UT basis as given in literature in terms of our kinematics and our propagator definitions. This means that, in addition to swapping the labels of the external momenta, we should also make sure the loop momenta are aligned between the two families.

For the pentagon-boxes, the relationship between the loop momenta (Abreu et al used `l`, we use `k`) is `l1->+k1, l2->-k2`, while for the non-planar hexa-boxes, it is `l1->+k1, l2->-k1-k2`. In principle, supplying these maps to LiteRed should not be necessary - it should be able to find the jExtRules between the two families. In practice, it seems that we might need to help it by supplying this input, otherwise not all jExtRules are found. Curiously, modifying these maps does not change the jExtRules.

There is, however, another subtlety which *is* affected by the propagator misalignment. Note that we define the extra-dimensional scalar products as:
`
MU[i, j] = - k_i^(d-4) . k_j^(d-4)
`
whereas in Abreu et al, they are:
`
Muij = + l_i^(d-4) . l_j^(d-4)
`. Therefore, taking into account both the misalignment of the propagators and the different sign of MU, the translation *for the penta-boxes* is:
```
Mu11 = - MU[1, 1] 
Mu22 = - MU[2, 2]	(* pentagon-box MU[i, j] translation rules *)
Mu12 = + MU[1, 2]
```
For the *non-planar hexa-boxes*, these rules are more complicated:
```
Mu11 = - MU[1, 1]	(* non-planar hexa-box MU[i, j] translation rules *)
Mu22 = - MU[1, 1] - MU[2, 2] - 2*MU[1, 2] 
Mu12 = + MU[1, 1] + MU[1, 2]
```
*These rules could be simplified if we changed the propagator definitions of the hexa-boxes on the k2 branch, but this is dangerous and would probably require adjusting other hard-coded bits of the code.*

 For the 1-loop pentagon1M, the `swap15` we use actually ensures that the `{2, 3, 4, 5}` permutation from Abreu et al is equivalent to our main perm `pentagon1Mp54321` and `l1->+k1`. So the translation rule is simply:
```
Mu11 = - MU[1, 1]	(* 1L pentagon1M MU[i, j] translation rules *)
```

Where does this become important? For example, when performing a check on a MI which has an explicit Muij insertion (see below for various checks on integrals), we need to use the above rules to express the extra-dimensional scalar products in our notation. However, this was not necessary when translating the Abreu et al basis into our kinematics, because the files supplied where given with Muij expressed in terms of sij's and propagators. Therefore, these translation rules were not needed there. **Ask Simone to expand on this**

**Note: As of 10/03/23, I have not yet worked implemented the UT basis for the non-planar double-pentagon. The non-planar hexagon-box should be ready to go, but I have never actually calculated an amplitude with it, so this needs to be checked further.**

## Penta-triangles (t521)

The penta-triangle topologies are a special case because, even though they appear when generating the amplitude, they can be mapped onto the pentagon-boxes. This is done at the stage of the IBP relations. As such, we don't need to express this family in terms of the pentagon functions, but we still need to have a basis for it so that we can generate the IBPs in the first place. 

### Defining the MainTopo

This can again be done by using LiteRed's jExtRules. In practice, I took the basis definitions for all the t431 permutations I needed for 2q2aW amplitudes, and I found the jExtRules between them and a single t521MZZZZ/t521ZMZZZ/t521ZZZM permutation that we define as the 'main topo'. This choice is arbitrary, but I chose a permutation that actually appears in the 2q2aW amplitude. The hard-coded choices for the t521 main topos are:
```
t521MZZZZp53214
t521ZMZZZp35124 (* chosen t521 main topos *)
t521ZZZZMp21345
```

Once the jExtRules have been found, we can then express the t431 MIs as t521 integrals in their corresponding main topos. **(Do not forget to include the normalisations as well!)** In general, this will lead to more integrals than the true number of masters for the t521 families. Therefore, as a last step we perform a mock IBP reduction to determine the final set of MIs. We also check that these MIs are indeed UT by deriving the DEs and checking their properties.

### 1-loop^2 t521ZZZZM integrals

For the t521ZZZZM 'subfamily', not all MIs can be obtained in the way described above. This is because we have a special subset of these integrals which are 1-loop^2, with the massive leg being on the triangle side, so that the integral becomes (bubble in s5)*(pentagon1M). Note that this is not possible for t521MZZZZ, t521ZMZZZ or the fully massless 5pt case, because it is not possible to formulate a 'bubble\*pentagon'-type integral without the mass.

To take care of this special case, we generate candidate masters by multiplying the 13 1-loop pentagon1M MIs by an s5-bubble. These 2L candidates are added onto the set of integrals translated from the t431 families, as described above. In this way, we obtain 71 + 13 = 84 candidate MIs. As in the previous section, we then perform a mock IBP reduction to determine the final set of MIs (which for t521ZZZZM is 69) and check the DEs.

## Expanding master integrals onto pentagon functions
Once we're confident that we have the UT basis for all needed permutations (including the t521 families), we need to express all the masters integrals in terms of `F[i, j]`, constants, square roots and tr5. In other words, we need to match our MIs onto the expressions provided in [PentagonFunctions-m1-datafiles](https://gitlab.com/pentagon-functions/PentagonFunctions-m1-datafiles).

This is more or less the reverse of the process which we employed to translate the Abreu et al UT bases in our notation, except that now we're working with MIs in many permutations, not just the 'main perm'. First, we need to apply the inverse of `swap15` to express our permutation, e.g. `t431MZZZZp51324` in the p1^2 kinematics. Then, we need to make sure that the propagators are aligned as well (this involves flipping the diagram upside down). We can then read off the corresponding permutation for the `mi[fam][{i, j, k, l}][i]` notation, where in perm the first leg is omitted, because it is always the massive p1. For the pentagon-boxes we read the permutation anticlockwise, starting from p1 (see Fig 1. of 2005.04195). This is achieved by applying the `translate[fam]` rul	e. *Again, note that this applies only to the penta-boxes. For the non-planar hexa-box, we still start at p1, but proceed clockwise - see Fig. 1 of 2107.14180.*

If we don't perform this translation task correctly, we will obtain incorrect expansions of the MIs int terms of pentagon functions. This was a mistake I made at one point when computing a 2L 2q2aW +-+-1 amplitude, where the 1/eps^3 contained incorrect F[i, j].

For reasons explained in the previous section, it might happen that our reduced amplitude will contain some t521ZZZZM MIs. These integrals cannot be mapped onto t431 ones, however we can obtain their expansion to pentagon functions by considering the product of the expansions of the s5-bubble and the right pentagon1M integral.

## Debugging

If things go wrong (and they almost certainly will...), there are various checks we can perform in order to help ourselves narrow down the problem.

### Comparing results obtained in multiple ways

1. **AMFlow evaluation ** -
pick a simple integral in the j[...] notation and evaluate it with AMFlow. Pick a phase-space point which is in the s12=t45 channel. 
Remember that to restore our normalisation of Feynman integrals, results from AMFlow have to be multiplied by a factor of `Exp[EulerGamma*eps]` *per loop*.

2. **Reduction onto Pentagon functions** -
reduce the integrals through 	GMN_alt/GMN onto pentagon functions. The result is written in terms of mono[...], which contain F[i, j], constants, square roots and tr5. The F[i, j] can be evaluated with `PentagonFunctions++`, but remember that the this package uses **(note the ordering)**:
$$
	\{p_1^2, t_{12}, t_{23}, t_{34}, t_{45}, t_{15}\}
$$
which in our Mandelstams is:
$$
	\{s_5, s_{45}, s_{34}, s_{23}, s_{12}, s_{15}\},
$$
so make sure to use the correct numerical values for the invariants. 
$tr_5$ can be evaluated through its expression in terms of momentum twistors, which can be obtained by solving the system $\{s_{12} = c_1, \dots \, \}$. Remember to choose the solution which insures that $\text{Im}(tr_5) > 0$, which is required for `PentagonFunctions++`.


3. **Reduction onto root-free MIs through GMN_alt/GMN** - 
the result can be expressed back in terms of `j[...]` using the basis definitions. This can then be evaluated using `AMFlow`.

4. **Reduction onto root-free MIs - no mappings used ('alt' only)** - 
same as above, but without using the mappings between MIs. This allows us to isolate the problem - does it occur already at the level of basis definitions or only when MIs from different families are mapped onto each other? This is impossible in the 'non-alt' approach as there all the families are considered together.

### Permuting the MI basis onto a different perm and checking with AMFlow
We start with the definition of a permuted master integral:
$$
	\text{mi}[\sigma(\text{fam})](X) = \text{mi}[\text{fam}](\sigma(X))\,,
$$
where $X$ is a phase-space point. Hence, we can verify whether the basis has been permuted correctly by evaluating the original basis at a permuted PS-point (RHS) and comparing it with the permuted basis evaluated at the original PS-point (LHS).

In particular, if we're trying to debug an incorrect $tr_5$ sign, we should take $\sigma$ to be an odd permutation, such that $\sigma(tr_5) = -tr_5\,$.
Done:
1. Store root vectors in Phi+ in a list called roots. Store all linear combinations of root vectors with coefficients of 1 in a list called sym_vec_combs
2. Calculate and store the Ualpha under which the possible root vector combinations are not stable.

To do:
1. Make matrices for the root vectors, torus, and unipotent root groups
8. Find out if there is a way to tell which orbits lie in a family with one another based on their root vector combination, rather than by calculating the B-orbits by brute force.
2. Calculate all possible B-Orbits
3. For each B-orbit, pull the coefficients on each root vector. Test if each coefficient is zero, and add a list of all coefficients of 0 to a list that stores the root vectors in the zero set. Do the same for coefficients that have only entries from T and none from the Ualpha
4. __Automate the discovery of algebraic relationships between the coordinates to give the polynomial defining equations in the zero set and the nonzero set, probably using an algorithm in elimination theory. "The basic idea of the implicitization problem is to convert the parametrization into definining equations for V" (Cox Little O'Shea pp. 133. N.t. they use V to refer to the polynomial zero set)__
5. __Find a way to tell when an "orbit" is infinite. It's probably when there's an algebraic dependence relation among the coordinates that involves a coefficient on one of the root vectors (in addition to or instead of the variables standing for an arbitrary scalar entry in T or Ualpha). Find a rigorous way of stating and defending this. Modality is probably the useful concept here. Since the orbits are varieties, is it possible to calculate the modality of the individual "orbits" to see if they turn out to be families of orbits? Then an orbit family would have modality 1, and a singleton orbit would have modality zero.__
6. __Automate the reverse inclusion calculations in Python. Is it possible that it actually isn't necessary to check the reverse inclusion? When would it happen that the B-orbit is in the intersection of the zero set and the nonzero set but not vice versa?__
7. __Find a way to tell when we have calculated all the B-orbits. Automate proof of nilradical exhaustion. Is there a way we can tell how many infinite "orbits" and how many singleton "orbits" there are before we calculate them? If so, we could just keep count and do a proof of exhaustion of n when we've got enough.__
9. Put in the option to display a B-orbit in equation form

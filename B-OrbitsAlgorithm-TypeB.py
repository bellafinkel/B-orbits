#!/usr/bin/env python
# coding: utf-8

# Done:
# 1. Store root vectors in Phi+ in a list called roots. Store all linear combinations of root vectors with coefficients of 1 in a list called sym_vec_combs
# 2. Calculate and store the Ualpha under which the possible root vector combinations are not stable.
# 
# To do:
# 1. Make matrices for the root vectors, torus, and unipotent root groups
#     1.i Change the data type of the root vectors from numpy arrays to matrices using np.asmatrix https://numpy.org/doc/stable/reference/generated/numpy.asmatrix.html
#     1.ii Write make_rvecs_A(l)
#     
# 8. Find out if there is a way to tell which vectors lie in the same orbit based on their root vector combination, rather than by calculating the B-orbits by brute force.
# 2. Calculate all possible B-Orbits
# 3. For each B-orbit, pull the coefficients on each root vector. Test if each coefficient is zero, and add a list of all coefficients of 0 to a list that stores the root vectors in the zero set. Do the same for coefficients that have only entries from T and none from the Ualpha
# 4. __Automate the discovery of algebraic relationships between the coordinates to give the polynomial defining equations in the zero set and the nonzero set, probably using an algorithm in elimination theory. "The basic idea of the implicitization problem is to convert the parametrization into definining equations for V" (Cox Little O'Shea pp. 133. N.t. they use V to refer to the polynomial zero set)__
# 5. __Find a way to tell when an "orbit" is infinite. It's probably when there's an algebraic dependence relation among the coordinates that involves a coefficient on one of the root vectors (in addition to or instead of the variables standing for an arbitrary scalar entry in T or Ualpha). Find a rigorous way of stating and defending this. Modality is probably the useful concept here. Since the orbits are varieties, is it possible to calculate the modality of the individual "orbits" to see if they turn out to be families of orbits? Then an orbit family would have modality 1, and a singleton orbit would have modality zero.__
# 6. __Automate the reverse inclusion calculations in Python. Is it possible that it actually isn't necessary to check the reverse inclusion? When would it happen that the B-orbit is in the intersection of the zero set and the nonzero set but not vice versa?__
# 7. __Find a way to tell when we have calculated all the B-orbits. Automate proof of nilradical exhaustion. Is there a way we can tell how many infinite "orbits" and how many singleton "orbits" there are before we calculate them? If so, we could just keep count and do a proof of exhaustion of n when we've got enough.__
# 9. Put in the option to display a B-orbit in equation form
# 10. Implement class structure
# 11. Add error handling

# In[2]:


import numpy as np
from itertools import combinations
from collections import OrderedDict


# In[3]:


"""
Make matrices for the root vectors for the positive roots in type B (i.e., the root vectors in the nilradical).
Store the root vector matrices in a dictionary called rvecs where the keys are the eigenvector names given on page 131
of Erdmann and the values are the eigenvectors as numpy arrays
"""
def make_rvecs_Bv1(l):
    rvecs = {}

    M = np.zeros((2*l+1,2*l+1),dtype = int)
    # Add root vectors b_i for roots of the form epsilon_i to rvecs
    for i in np.arange(1,l+1,1):
        # b_ij
        M[i,0] = 1
        M[0,l+i] = -1
        hilo = "b" + str(i)
        rvecs.setdefault(hilo,M)
        M = np.zeros((2*l+1,2*l+1),dtype = int)
        #print(i)
        
    for i in np.arange(1,l,1):
        for j in np.arange(1,l+1,1):
            # p_ij
            if i<j:
                M[i,l+j] = 1
                M[j,l+i] = -1
                hilo = "p" + str(i) + str(j)
                rvecs.setdefault(hilo,M)
                M = np.zeros((2*l+1,2*l+1),dtype = int)
            # m_ij
            if i != j and i<j:
                M[i,j] = 1
                M[l+j,l+i] = -1
                hilo = "m" + str(i) + str(j)
                rvecs.setdefault(hilo,M)
                M = np.zeros((2*l+1,2*l+1),dtype = int)
    #print(i,j)
    return rvecs

rvecs = make_rvecs_Bv1(3)
print(len(rvecs))


# In[4]:


# Take the type and rank of the Lie group
L = RootSystem(["B",3]).root_space()
roots=list(L.positive_roots())
print(roots)
#roots[0].to_vector()
roots = [root.to_vector() for root in roots]
roots


# In[5]:


"""
Make matrices for the root vectors for the positive roots in type B (i.e., the root vectors in the nilradical).
Store the root vector matrices in a dictionary called rvecs where the keys are the eigenvector names given on page 131
of Erdmann and the values are the eigenvectors as numpy arrays
"""
def make_rvecs_B(l):
    rvecs = {}

    M = np.zeros((2*l+1,2*l+1),dtype = int)
    # Add root vectors b_i for roots of the form epsilon_i to rvecs
    for i in np.arange(1,l+1,1):
        # b_ij
        M[i,0] = 1
        M[0,l+i] = -1
        hilo = "b" + str(i)
        coef = np.zeros(l, dtype = int) # Array of length the number of roots in the base
        coef[i-1:l] = 1
        coef = tuple(coef)
        #print(i,coef)
        # Assign the list consisting of the root vector M and the root coef to the key hilo in the dictionary rvecs
        #rvecs.setdefault(hilo,[coef,M])
        rvecs.setdefault(coef,M)
        M = np.zeros((2*l+1,2*l+1),dtype = int)
        #print(i)
        
    for i in np.arange(1,l,1):
        for j in np.arange(1,l+1,1):
            # p_ij
            if i<j:
                M[i,l+j] = 1
                M[j,l+i] = -1
                hilo = "p" + str(i) + str(j)
                coef = np.zeros(l, dtype = int) # Array of length the number of roots in the base
                coef[i-1:j] = 1
                coef[j-1:l] = 2
                coef = tuple(coef)
                #print(i,j,coef)
                #rvecs.setdefault(hilo,[coef,M])
                rvecs.setdefault(coef,M)
                M = np.zeros((2*l+1,2*l+1),dtype = int)
            # m_ij
            if i != j and i<j:
                M[i,j] = 1
                M[l+j,l+i] = -1
                hilo = "m" + str(i) + str(j)
                coef = np.zeros(l, dtype = int) # Array of length the number of roots in the base
                coef[i-1:j-1] = 1 
                coef = tuple(coef)
                #rvecs.setdefault(hilo,[coef,M])
                rvecs.setdefault(coef,M)
                M = np.zeros((2*l+1,2*l+1), dtype = int)
    #print(i,j)
    return rvecs

rvecs = make_rvecs_B(3)
print(rvecs)

# Check that make_rvecs_B(l) computes the same roots as Sage (it looks like it does)
bee = []
for value in rvecs.values():
    bee.append(value[0])
set(bee)
x=[tuple(root) for root in roots]
set(x)==set(bee)
# In[73]:


"""
Make a matrix for the torus in type B
"""
# Make a matrix for the torus given group type and rank
# Input for Gtype may be A,B,C,D
# None of these are accurate except the one for type B
def make_T(Gtype,l):
    if Gtype == 'A':
        n = l + 1
        #syms=np.array([var("q"+str(i)) for i in range(n)])
        #T = np.diag(syms)
    elif Gtype == 'B':
        n = 2*l+1 
        syms = np.array([var("q"+str(i)) for i in range(l)], dtype = object)
        symsinv = np.array([var("q"+str(i))^(-1) for i in range(l)], dtype = object)
        syms = np.insert(syms,0,1)
        syms = np.insert(syms,l+1,symsinv)
        #print(syms)
        T = np.diag(syms)
    elif Gtype == 'C':
        n = 2*l
        #syms=np.array([var("q"+str(i)) for i in range(n)])
        #T = np.diag(syms)
    elif Gtype == 'D':
        n = 2*l
        #syms=np.array([var("q"+str(i)) for i in range(n)])
        #T = np.diag(syms)
    return T


# In[8]:


"""
Make matrices for each unipotent root group of the positive roots for 
type B
"""
def make_Ua(l):
    Ua = {}
    i = 0
    for key, value in rvecs.items():
        # Make the numpy array a Sage matrix so that it can be properly 
        # exponentiated (equiv. of MatrixExp), then take the matrix exponential
        name = "t" + str(key).replace("(", "").replace(")", "").replace(" ", "").replace(",", "")
        #print(name)
        name = var(name)
        U = matrix(value*name)
        U = e^U
        Ua.setdefault(key,U)
    return Ua
U = make_Ua(3)


# In[9]:


"""
Create a list storing all combinations we want to sum of root vectors as symbols, 
excluding the empty sum but including single root vectors
"""

roots = np.array(list(rvecs.keys()))
for t in range(len(roots)):
    roots[t] = list(roots[t])
m = len(roots)

sym_vec_combs = []
for i in range(0, m):
    sym_vec_combs += list(combinations(roots,i))
print(len(sym_vec_combs))


#         # If an entry involves some expression of t, change it to 
#         # another variable. This doesn't work and we also probably
#         # don't want it
#         for i in range(2*l+1):
#             for j in range(2*l+1):
#                 if U[i,j] != 0 and U[i,j] != 1:
#                     U[i,j] = name
#                     Ua[key] = U

# In[10]:


roots = list(rvecs.keys())
for t in range(len(roots)):
    roots[t] = list(roots[t])
print(roots)
#print(type(x[0]))


# In[11]:


# Take the type and rank of the Lie group
L = RootSystem(["B",3]).root_space()
roots=list(L.positive_roots())
print(roots)
#roots[0].to_vector()
roots = [root.to_vector() for root in roots]
roots


# In[12]:


"""
Create a list of all combinations we want to sum of root vectors as symbols, 
excluding the empty sum but including single root vectors
"""
x_list=roots
sym_vec_combs = []
for i in range(1, len(x_list) + 1):
    sym_vec_combs += list(combinations(x_list,i))
print(len(sym_vec_combs))


# In[13]:


"""
Fix a root x in Phi+, then loop through each root xi in Phi+ to check
whether x+xi is in Phi+. If x+xi is in Phi+, then x is not stable 
under Ui and we place Ui in a list with the same indexing as roots.

Index zero corresponds to alpha1, etc.
"""
Ualpha_list=[ [] for _ in range(len(roots)) ]
N = range(len(roots))
for f in N:
    for t in N:
        if roots[f]+roots[t] in roots:
            Ualpha_list[f].append(roots[t])

print(Ualpha_list)


# In[14]:


tuple(Ualpha_list[0][0])


# In[15]:


"""
Umat is basically Ualpha_list but with matrices
"""
Umat = []

for i in range(len(Ualpha_list)):
    dum = []
    for j in range(len(Ualpha_list[i])):
        #print(Ualpha_list[i][j])
        #print(type(U[tuple(Ualpha_list[0][0])]))
        dum.append(np.array(U[tuple(Ualpha_list[i][j])]))
    Umat.append(dum)


# In[16]:


Umat[0]


# In[17]:


print(type(np.array(sym_vec_combs[10])))


# In[18]:


Ualpha_list


# In[19]:


"""
Create a list xUlist to store the Ualpha under which the possible 
root vector combinations (stored in sym_vec_combs) are not stable.
Shares index with sym_vec_combs.
"""
xUlist = [ [] for _ in range(len(sym_vec_combs)) ]
xUlist_mat = [ [] for _ in range(len(sym_vec_combs)) ]
for k in range(len(sym_vec_combs)):
    for i in range(len(Ualpha_list)):
        if roots[i] in sym_vec_combs[k]:
            xUlist[k].extend(Ualpha_list[i])
            # xUlist_mat to store the Ualpha as matrices under which the possible 
            # root vector combinations (stored in sym_vec_combs) are not stable
            xUlist_mat[k].extend(Umat[i])


# In[20]:


xUlist_mat[0][0]


# In[56]:


"""
THIS WOULD BE BETTER AS A DICTIONARY

Save the products of U_alpha as matrices in a list (This is important)
THESE DIFFER FROM THE ONES CALCULATED IN B_OrbitAutomation_v3.1, MAYBE BECAUSE THE MATRICES ARE BEING MULTIPLIED IN
DIFFERENT ORDERS. WHY DOESN'T ORDER MATTER FOR THE PURPOSE OF CALCULATING B-ORBITS?
"""

Uprod = []
l = 7

for i in range(len(xUlist_mat)):
    dum = np.eye(l, dtype = int)
    for j in range(len(xUlist_mat[i])):
        new_mat = xUlist_mat[i][j]
        dum = np.matmul(dum,new_mat)
    Uprod.append(dum)
Uprod[1]


# In[22]:


Uprod[2]


# In[47]:


type(sym_vec_combs[0][0])


# In[36]:


rvecs[tuple(sym_vec_combs[0][0])]


# In[113]:


"""
THIS WOULD BE BETTER AS A DICTIONARY

Save the sums of the root vector combinations

using the matrices stored as values in the rvecs dictionary
"""
l = 7
num = len(sym_vec_combs)
vsums = []

for i in range(num):
    dum = np.zeros(l, dtype = object)
    for j in range(len(sym_vec_combs[i])):
        new_vec = rvecs[tuple(sym_vec_combs[i][j])]
        dum = dum + new_vec #np.matrix.sum(dum,new_vec)
    vsums.append(dum)
vsums[99]


# In[96]:


TT = np.array(T,dtype=object)
print(TT.shape)
print(Uprod[0].shape)


# In[97]:


TT


# In[98]:


T = matrix(make_T('B',3))
IT = T.inverse()


# In[105]:


matrix(Uprod[0]).inverse()


# In[114]:


"""
Create a list of all possible B-orbits (potentially with repeats)
"""
B_list = []

for i in range(len(vsums)):
    B = T @ Uprod[i] @ vsums[i] @ matrix(Uprod[i]).inverse() @ IT
    B = simplify(B)
    B_list.append(B)
    print(i)


# In[115]:


B_list[0]


# # Scratch

# In[23]:


num = len(xUlist_mat)
Uprod = np.empty(num)

for i in range(num):
    Uprod[num] = np.matmul(xUlist_mat[i])


# In[ ]:


"""
Save the products of U_alpha as matrices (This is important) in a DICTIONARY (DO THIS) as values 
where the keys are the sums of root vectors as symbols
"""
# First save the products of U_alpha as matrices in the list Umatlist
Uprod = []

for i in range(len(xUlist)):
    dum = 1
    for j in range(len(xUlist[i])):
        dum *= xUlist[i][j]
    Uprod.append(dum)
print(xUlist[1])
Uprod[10]


# In[ ]:


print(type(sym_vec_combs[0]))


# In[ ]:


print(type(xUlist[0][0]))


# In[ ]:


import itertools


# In[ ]:


print(type(xUlist[0][0]))


# In[ ]:


np.array(U[tuple(xUlist_mat[0][0])


# In[ ]:


for i in range(len(Ulist)):
    dum = 1
    for j in range(len(Ulist[i])):
        dum *= Ulist[i][j]
    Uprod.append(dum)
print(Ulist[1])
Uprod[1]


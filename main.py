import numpy as np
#import scipy.sparse as sparse
import scipy as sp
from scipy import linalg
pi = np.pi;

def matmult(a,b):
    zip_b = zip(*b)
    # uncomment next line if python 3 :
    # zip_b = list(zip_b)
    return [[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b))
             for col_b in zip_b] for row_a in a]

def formStiffnessMassTimoshenkoBeam(GDof,numberElements, elementNodes,numberNodes,xx,C,P,rho,I,A):
    stiffness=np.zeros([GDof,GDof])
    mass=np.zeros([GDof,GDof])
    force=np.zeros([GDof,1])
    # stiffness matrix
    gaussLocations=np.array([0.577350269189626, -0.577350269189626])
    gaussWeights=np.ones([2,1])
    # bending contribution for stiffness matrix
    for e in range(0,numberElements):
        indice= elementNodes[e,:]
        elementDof =np.append( indice, indice+numberNodes)
        indiceMass=indice+numberNodes
        ndof=np.size(indice)
        length_element=xx[indice[1]]-xx[indice[0]]
        detJacobian=length_element/2.0
        invJacobian=1.0/detJacobian;
        for q in range(0,np.size(gaussWeights)):
            #print "q=", q
            pt=gaussLocations[q];
            [shape,naturalDerivatives]=shapeFunctionL2(pt);##################
            Xderivatives=naturalDerivatives*invJacobian;###########
            # B matrix
            B=np.zeros([2,2*ndof]);
            for g in range(ndof,2*ndof):
                B[0,g] = Xderivatives[g-ndof];##################
            #K
            ######stiffness[elementDof,elementDof] = stiffness[elementDof,elementDof]+ matmult(B.transpose(),B)*gaussWeights[q]*detJacobian*C[0,0];
            tmp = matmult(B.transpose(),B)*gaussWeights[q]*detJacobian*C[0,0];
            for i in range(np.size(elementDof)):
                for j in range(np.size(elementDof)):
                    stiffness[elementDof[i],elementDof[j]] += tmp[i,j];
            #############
            #############
            # force[indice]=force[indice]+shape*P*detJacobian*gaussWeights[q];
            tmp = shape*P*detJacobian*gaussWeights[q];

            for i in range(np.size(indice)):
                force[indice[i]]+= tmp[0,i];
            #############
            # #############mass[indiceMass,indiceMass]=mass[indiceMass,indiceMass]+matmult(shape.transpose(),shape)*gaussWeights[q]*I*rho*detJacobian;
            tmp = matmult(shape.transpose(),shape)*gaussWeights[q]*I*rho*detJacobian;
            for i in range(np.size(indiceMass)):
                for j in range(np.size(indiceMass)):
                    mass[indiceMass[i],indiceMass[j]] += tmp[i,j];

            #############mass[indice,indice]=mass[indice,indice]+matmult(shape.transpose(),shape)*gaussWeights[q]*A*rho*detJacobian;
            tmp = matmult(shape.transpose(),shape)*gaussWeights[q]*A*rho*detJacobian;
            for i in range(np.size(indice)):
                for j in range(np.size(indice)):
                    mass[indice[i],indice[j]] += tmp[i,j];
            #############
    #  shear contribution for stiffness matrix

    gaussLocations=np.array([ 0 ]);
    gaussWeights = np.ones([1,1])
    gaussWeights[0, 0]=2
    for e in range(0,numberElements):
        indice= elementNodes[e,:]
        # elementDof=[ indice indice+numberNodes];
        elementDof =np.append( indice, indice+numberNodes)
        ndof=np.size(indice)

       # ndof=length(indice);
        length_element=xx[indice[1]]-xx[indice[0]]
        detJ0=length_element/2.0
        invJ0=1.0/detJ0
        for q in range(np.size(gaussWeights)):
            pt=gaussLocations[q]
            [shape,naturalDerivatives]=shapeFunctionL2(pt)
            Xderivatives=naturalDerivatives*invJacobian;


            #B
            B=np.zeros([2, 2*ndof])
            print "Xderivatives = ", Xderivatives
            for i in range(ndof):
                B[1,i] = Xderivatives[i]
                B[1,ndof+i] = shape[0,i]

            #K
            tmp = matmult(B.transpose(), B) * gaussWeights[q]*detJacobian*  C[1, 1]
            #stiffness(elementDof,elementDof)=...
             #   stiffness(elementDof,elementDof)+ B'*B*gaussWeights(q)*detJacobian*C(2,2);
            for i in range(np.size(elementDof)):
                for j in range(np.size(elementDof)):
                    stiffness[elementDof[i], elementDof[j]] += tmp[i, j]

    return stiffness,force,mass


def shapeFunctionL2(xi):
    # shape function and derivatives for L2 elements
    # shape : Shape functions
    # naturalDerivatives: derivatives w.r.t. xi
    # xi: natural coordinates (-1 ... +1)
    shape=np.array([[1-xi,1+xi]])  *0.5
    naturalDerivatives=np.array([-0.5,0.5]) #[-1;1]/2;
    return shape,naturalDerivatives

def CalcMCK_FEM(nodes, numberElements,r, rho, E, G ,Ks):
    # [mass,stiffness,xx,V,D,M,K]
    #................................................................
    # MATLAB codes for Finite Element Analysis
    # problem16.m
    # Timoshenko beam in bending
    # antonio ferreira 2008
    # clear memory
    # clear all
    # close all
    # E; modulus of elasticity
    # G; shear modulus
    # I: second moments of area
    # L: length of beam
    # thickness: thickness of beam
    # E=10e7;
    # poisson = 0.30;
    L = nodes[1]-nodes[0];

    # A = pi*r.^2;
    # dA = 2*pi*r * m;

    I = pi*r*r*r*r /4; # for circular intersection

    A = pi*r*r;
    EI=E*I;

    P = -1; # uniform pressure
    P = 0;
    # constitutive matrix
  # G=E/2/(1+poisson);
    #C=[ EI 0; 0 Ks*A*G];
    C = np.array([[EI, 0],[ 0, Ks*A*G]]);
    # mesh
    # numberElements = 100;

    nodeCoordinates = np.linspace(0,L,numberElements+1);
    xx=nodeCoordinates;
    elementNodes = np.zeros([numberElements,2]);
    for i in range(0, numberElements ):
        elementNodes[i,0]=i;####
        elementNodes[i,1]=i+1;####
    elementNodes = np.int_(elementNodes)
    # generation of coordinates and connectivities
    #numberNodes=np.size(xx,2);
    numberNodes=np.size(xx);
    # GDof: global number of degrees of freedom
    GDof=2*numberNodes;
    # computation of the system stiffness matrix
    [stiffness,force, mass]=formStiffnessMassTimoshenkoBeam(GDof,numberElements,elementNodes,numberNodes,xx,C,P, rho, I, A);
    # boundary conditions (simply-supported at both bords)
    #fixedNodeW =[1 ; numberNodes];
    #fixedNodeTX=[];
    # boundary conditions (clamped at both bords)
    fixedNodeW = np.array([0 , numberNodes-1]);
    fixedNodeTX=fixedNodeW;
    # boundary conditions (cantilever)
    fixedNodeW =  np.array([ numberNodes - 1 ]);#[1]
    fixedNodeTX = np.array([ numberNodes - 1 ]);#[1]
    prescribedDof = np.array([fixedNodeW , fixedNodeTX+numberNodes]);
    print "prescribedDof = ",prescribedDof
    # solution
    ''''
    displacements=solution(GDof,prescribedDof,stiffness,force);
    # output displacements/reactions
    # outputDisplacementsReactions(displacements,stiffness,...
    # GDof,prescribedDof)
    U=displacements;
    ws=1:numberNodes;


    '''
 #   a= np.array(range(0,10))
   # print a
    #activeDof=setdiff([1:GDof],[prescribedDof]);
    #activeDof=np.setdiff1d( np.array(range(0,GDof)), [prescribedDof]);
   # K = stiffness(activeDof,activeDof);
    K = stiffness #[activeDof,activeDof]
    #print "K = ",K
  #  M = mass(activeDof,activeDof);
    M = mass;
    for k in range(np.size(prescribedDof)):
        for i in range(GDof):
            for j in range(GDof):
                if( i == prescribedDof[k] or j == prescribedDof[k]):
                    M[i,j] = 0
                    K[i,j] = 0

    for k in range(np.size(prescribedDof)):
        M[prescribedDof[k],prescribedDof[k]] = 1
        K[prescribedDof[k],prescribedDof[k]] = 1

    np.savetxt("K.txt",K)
    np.savetxt("M.txt",M)
  #  [V,D]=eig(stiffness(activeDof,activeDof), mass(activeDof,activeDof));
    print "hello"
    invMK = linalg.inv(M)*K;
    np.savetxt("invMKpy.txt",invMK)
    V,D = np.linalg.eig(invMK)
    V,D = sp.linalg.eig(K,M)
#    V,D = sparse.linalg.eigs(invMK)
 #   print "hello"

#    print np.size(V)
    print "freq(Hz) = ",np.sqrt(V)/2.0/pi
 #   print "D = ",D
#    print "K eig = ", np.linalg.eig(K)
    #eigvals()




    return mass,stiffness,xx,V,D,M,K


ap = 0.002;
L1 = 0;
L2 = ap;
L3 = 0.06;
nodes = np.array([L1, L3]);


N_sampling = 10;192;
N = N_sampling;

rho = 7850;
E = 200e9; # 200 GPa
G = 79e9;
v = E/2/G-1;
Ks = 6*(1+v)/(7+6*v); # for circular cross-section
r = 0.006;
a = np.zeros([2,6]);
a = np.array([[2,3,5],[5,6,8]]);
b = np.array([2,3,4]);

print "***"
print a
print a.transpose()
print nodes
print np.transpose(nodes)
print a

[mass,stiffness,xx,V,D,M,K] = CalcMCK_FEM(nodes, N_sampling,r, rho, E, G ,Ks)


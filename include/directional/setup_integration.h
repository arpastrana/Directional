// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_SETUP_INTEGRATION_H
#define DIRECTIONAL_SETUP_INTEGRATION_H

#include <queue>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_topology.h>
#include <directional/tree.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/dcel.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/combing.h>

namespace directional
{
  struct IntegrationData
  {
    int N;  //num of parametric functions
    int n;  //actual dimension of problem
    Eigen::MatrixXi symmFunc; //Symmetry function tying the n dofs to the full N
    Eigen::MatrixXi intFunc;  //function spanning integers
    Eigen::SparseMatrix<double> vertexTrans2CutMat;
    Eigen::SparseMatrix<double> constraintMat;
    Eigen::SparseMatrix<double> symmMat;  //in fact, using it for the general reduction of degrees of freedom
    Eigen::SparseMatrix<double> intSpanMat;  //in case some integers are constrained to lie on a lattice
    Eigen::SparseMatrix<double> singIntSpanMat;  //the layer for the singularities
    Eigen::VectorXi constrainedVertices;
    Eigen::VectorXi integerVars;
    Eigen::MatrixXi face2cut;
    
    Eigen::VectorXi fixedIndices;  //the translation fixing indices
    Eigen::VectorXd fixedValues;   //translation fixed values
    Eigen::VectorXi singularIndices;   //the singular-vertex indices
    
    //integer versions, for pure seamless parameterizations
    Eigen::SparseMatrix<int> vertexTrans2CutMatInteger;
    Eigen::SparseMatrix<int> constraintMatInteger;
    Eigen::SparseMatrix<int> symmMatInteger;  //in fact, using it for the general reduction of degrees of freedom
    Eigen::SparseMatrix<int> intSpanMatInteger;  //in case some integers are constrained to lie on a lattice
    Eigen::SparseMatrix<int> singIntSpanMatInteger;  //the layer for the singularities
    
    double lengthRatio;     //global scaling of
    bool integralSeamless;  //do not do translational seamless.
    bool roundSeams;//do not do translational seamless.
    bool verbose;
    bool localInjectivity;
    
    IntegrationData():lengthRatio(0.02), integralSeamless(false), roundSeams(true), verbose(false), localInjectivity(false){}
    ~IntegrationData(){}
  };
  
  IGL_INLINE Eigen::MatrixXi sign_symmetry(int N){
    assert(N%2==0);
    Eigen::MatrixXi symmFunc(N,N/2);
    symmFunc<<Eigen::MatrixXi::Identity(N/2,N/2),-Eigen::MatrixXi::Identity(N/2,N/2);
    std::cout<<"symmFunc: "<<symmFunc<<std::endl;
    return symmFunc;
  }

  IGL_INLINE Eigen::MatrixXi default_period_jumps(int n){
    return Eigen::MatrixXi::Identity(n,n);
  }
  
  

  
  
  // Setting up the seamless integration algorithm
  // Input:
  // symmFunc:    #N x n matrix that describes the transformation between the physical "n"" dofs and the N functions integration into. This should not depend on the mesh.
   // intFunc:    n x n matrix that describes any special relation between the translational jumps. If you don't know what that means, just give Eigen::Identity(d)
  //  wholeV:     #V x 3 vertex coordinates
  //  wholeF:     #F x 3 face vertex indices
  //  EV:         #E x 2 edges to vertices indices
  //  EF:         #E x 2 edges to faces indices
  //  FE:         #F x 3 faces to edges tindices
  // matching:    #E matching function, where vector k in EF(i,0) matches to vector (k+matching(k))%N in EF(i,1). In case of boundary, there is a -1. Most matching should be zero due to prior combing.
  // singVertices: list of singular vertices in wholeV.
  // Output:
  //  intData:         integration data subsequently used in directional::integrate();
  //  cutV:       the Vertices of the cut mesh.
  //  cutF:       The Faces of the cut mesh.
  
  IGL_INLINE void setup_integration(const Eigen::MatrixXi& symmFunc,
                                    const Eigen::MatrixXi& intFunc,
                                    const Eigen::MatrixXd& wholeV,
                                    const Eigen::MatrixXi& wholeF,
                                    const Eigen::MatrixXi& EV,
                                    const Eigen::MatrixXi& EF,
                                    const Eigen::MatrixXi& FE,
                                    const Eigen::MatrixXd& rawField,
                                    const Eigen::VectorXi& matching,
                                    const Eigen::VectorXi& singVertices,
                                    IntegrationData& intData,
                                    Eigen::MatrixXd& cutV,
                                    Eigen::MatrixXi& cutF,
                                    Eigen::MatrixXd& combedField,
                                    Eigen::VectorXi& combedMatching)
  {
    
    using namespace Eigen;
    using namespace std;
    
    //cutting mesh and combing field.
    cut_mesh_with_singularities(wholeV, wholeF, singVertices, intData.face2cut);
    combing(wholeV,wholeF, EV, EF, FE, intData.face2cut, rawField, matching, combedField, combedMatching);
   
    int N = symmFunc.rows();
    int n = symmFunc.cols();
    intData.N=N;
    intData.n=n;
    intData.symmFunc = symmFunc;
    intData.intFunc = intFunc;
    
    MatrixXi EFi,EH, FH;
    MatrixXd FEs;
    VectorXi VH, HV, HE, HF, nextH, prevH, twinH, innerEdges;
    
    // it stores number of edges per face, for now only tirangular
    VectorXi D = VectorXi::Constant(wholeF.rows(), 3);
    
    // mark vertices as being a singularity vertex of the vector field
    VectorXi isSingular = VectorXi::Zero(wholeV.rows());
    for (int i = 0; i < singVertices.size(); i++)
      isSingular(singVertices(i)) = 1;
    
    //cout<<"singVertices: "<<singVertices<<endl;
    
    intData.constrainedVertices = VectorXi::Zero(wholeV.rows());
    
    //computing extra topological information
    std::vector<int> innerEdgesVec; // collects ids of inner edges
    EFi = Eigen::MatrixXi::Constant(EF.rows(), 2, -1); // number of an edge inside the face
    
    /* used later for internal edges there is 1 or  -1 ie if two faces are adjacent then for a given edge we
     * will have 1 in the frst face and -1 in the second
     */
    FEs = Eigen::MatrixXd::Zero(FE.rows(), FE.cols());
    
    /*
     * here we collect information about position of an edge inside each face containing it. Each triangular face
     * has three edges of ids 0, 1, 2. So EFi(i, k) = j means that the edge i is inside the face k \in [0,1]
     * at the position j.
     */
    for(int i = 0; i < EF.rows(); i++)
    {
      for (int k = 0; k < 2; k++)
      {
        if (EF(i, k) == -1)
          continue;
        for (int j = 0; j < D(EF(i, k)); j++)
          if (FE(EF(i, k), j) == i)
            EFi(i, k) = j;
      }
    }
    
    // collect information about inner edges
    for(int i = 0; i < EF.rows(); i++)
    {
      if(EFi(i, 0) != -1)
        FEs(EF(i, 0), EFi(i, 0)) = 1.0;
      if(EFi(i,1) != -1)
        FEs(EF(i, 1), EFi(i, 1)) = -1.0;
      if ((EF(i, 0) !=-1) && (EF(i,1)!=-1))
        innerEdgesVec.push_back(i);
    }
    
    // copy the information into  Eigen vector
    innerEdges.resize(innerEdgesVec.size());
    for (int i = 0; i < innerEdgesVec.size(); i++)
      innerEdges(i) = innerEdgesVec[i];
    
    // compute the half-edge representation
    hedra::dcel(D, wholeF, EV, EF, EFi, innerEdges, VH, EH, FH, HV, HE, HF, nextH, prevH, twinH);
    
    // find boundary vertices and mark them
    VectorXi isBoundary = VectorXi::Zero(wholeV.rows());
    for (int i = 0; i < HV.rows(); i++)
      if (twinH(i) == -1){
        isBoundary(HV(i)) = 1;
        isSingular(HV(i)) = 0; //boundary vertices cannot be singular
      }
    
    
    
    /*for (int i=0;i<wholeV.rows();i++){
     if (isSingular(i))
     cout<<"vertex "<<i<<" is singular "<<endl;
     }*/
    
    // here we compute a permutation matrix
    vector<MatrixXi> constParmMatrices(N);
    MatrixXi unitPermMatrix = MatrixXi::Zero(N, N);
    for (int i = 0; i < N; i++)
      unitPermMatrix((i + 1) % N, i) = 1;
    
    // generate all the members of the permutation group
    constParmMatrices[0] = MatrixXi::Identity(N, N);
    for (int i = 1; i < N; i++)
      constParmMatrices[i] = unitPermMatrix * constParmMatrices[i - 1];
    
    // each edge which is on the cut seam is marked by 1 and 0 otherwise
    VectorXi isSeam = VectorXi::Zero(EV.rows());
    for(int i = 0; i < FE.rows(); i++)
    {
      for (int j = 0; j < 3; j++)
        if (intData.face2cut(i, j)) // face2cut is initalized by directional::cut_mesh_with_singularities
          isSeam(FE(i, j)) = 1;
    }
    
    // do the same for the half-edges, mark edges which correspond to the cut seam
    VectorXi isHEcut = VectorXi::Zero(HE.rows());
    for(int i = 0; i < wholeF.rows(); i++)
    {
      for (int j = 0; j < 3; j++)
        if (intData.face2cut(i, j)) // face2cut is initalized by directional::cut_mesh_with_singularities
          isHEcut(FH(i, j)) = 1; // FH is face to half-edge mapping
    }
    
    // calculate valency of the vertices which lay on the seam
    VectorXi cutValence = VectorXi::Zero(wholeV.rows());
    for(int i = 0; i < EV.rows(); i++)
    {
      if (isSeam(i))
      {
        cutValence(EV(i, 0))++;
        cutValence(EV(i, 1))++;
      }
    }
    
    
    //establishing transition variables by tracing cut curves
    VectorXi Halfedge2TransitionIndices = VectorXi::Constant(HE.rows(), 32767);
    VectorXi Halfedge2Matching(HE.rows());
    VectorXi isHEClaimed = VectorXi::Zero(HE.rows());
    
    // here we convert the matching that was calculated for the vector field over edges to half-edges
    for (int i = 0; i < HE.rows(); i++)
    {
      // HE is a map between half-edges to edges, but it does not carry the direction
      // EH edge to half-edge mapping
      Halfedge2Matching(i) = (EH(HE(i), 0) == i ? -combedMatching(HE(i)) : combedMatching(HE(i)));
      if(Halfedge2Matching(i) < 0)
        Halfedge2Matching(i) = (N + (Halfedge2Matching(i) % N)) % N;
    }
    
    int currTransition = 1;
    
    /*
     * Next steps: cutting mesh and creating map between wholeF and cutF
     */
    
    //cutting the mesh
    vector<int> cut2whole;
    vector<RowVector3d> cutVlist;
    cutF.conservativeResize(wholeF.rows(),3);
    for (int i = 0; i < VH.rows(); i++)
    {
      //creating corners whereever we have non-trivial matching
      int beginH = VH(i);
      int currH = beginH;
      
      //reseting to first cut or first boundary, if exists
      if (!isBoundary(i))
      {
        do
        {
          if (isHEcut(currH)!=0)
            break;
          currH=nextH(twinH(currH));
        } while (beginH!=currH);
      }
      else
      {
        do
        {
          if (twinH(currH)==-1)
            break;
          currH=nextH(twinH(currH));
        } while(twinH(currH)!=-1);
      }
      
      beginH = currH;
      
      do
      {
        if ((isHEcut(currH) != 0) || (beginH == currH))
        {
          cut2whole.push_back(i);
          cutVlist.push_back(wholeV.row(i));
        }
        
        for (int j = 0; j < 3; j++)
          if (wholeF(HF(currH), j) == i)
            cutF(HF(currH), j) = cut2whole.size() - 1;
        currH = twinH(prevH(currH));
      } while((beginH != currH) && (currH != -1));
    }
    
    cutV.conservativeResize(cutVlist.size(), 3);
    for(int i = 0; i < cutVlist.size(); i++)
      cutV.row(i) = cutVlist[i];
    
    //starting from each cut-graph node, we trace cut curves
    for(int i = 0;  i < wholeV.rows(); i++)
    {
      if (((cutValence(i) == 2) && (!isSingular(i))) || (cutValence(i) == 0))
        continue;  //either mid-cut curve or non at all
      
      //tracing curves until next node, if not already filled
      int beginH = VH(i);
      
      //reseting to first boundary
      int currH = beginH;
      
      if (isBoundary(i))
      {
        do
        {
          if (twinH(currH) == -1)
            break;
          currH = nextH(twinH(currH));
        } while(twinH(currH) != -1);
      }
      
      beginH = currH;
      
      int nextHalfedgeInCut = -1;
      do
      {
        //unclaimed inner halfedge
        if ((isHEcut(currH) != 0) && (isHEClaimed(currH) == 0) && (twinH(currH) != -1))
        {
          nextHalfedgeInCut = currH;
          Halfedge2TransitionIndices(nextHalfedgeInCut) = currTransition;
          Halfedge2TransitionIndices(twinH(nextHalfedgeInCut)) = -currTransition;
          isHEClaimed(nextHalfedgeInCut) = 1;
          isHEClaimed(twinH(nextHalfedgeInCut)) = 1;
          int nextCutVertex=HV(nextH(nextHalfedgeInCut));
          //advancing on the cut until next node
          while ((cutValence(nextCutVertex) == 2) && (!isSingular(nextCutVertex)) && (!isBoundary(nextCutVertex)))
          {
            int beginH = VH(nextCutVertex);
            int currH = beginH;
            int nextHalfedgeInCut = -1;
            do
            {
              //unclaimed cut halfedge
              if ((isHEcut(currH) != 0) && (isHEClaimed(currH) == 0))
              {
                nextHalfedgeInCut = currH;
                break;
              }
              currH=twinH(prevH(currH));
            } while (beginH != currH);
            Halfedge2TransitionIndices(nextHalfedgeInCut) = currTransition;
            Halfedge2TransitionIndices(twinH(nextHalfedgeInCut)) = -currTransition;
            isHEClaimed(nextHalfedgeInCut) = 1;
            isHEClaimed(twinH(nextHalfedgeInCut)) = 1;
            nextCutVertex = HV(nextH(nextHalfedgeInCut));
          }
          currTransition++;
        }
        currH = twinH(prevH(currH));
      } while((beginH != currH) && (currH != -1));
    }
    // end of cutting
    
    int numTransitions = currTransition - 1;
    //cout<<"numtransitions: "<<numTransitions<<endl;
    vector<Triplet<double> > vertexTrans2CutTriplets, constTriplets;
    vector<Triplet<int> > vertexTrans2CutTripletsInteger, constTripletsInteger;
    //forming the constraints and the singularity positions
    int currConst = 0;
    // this loop set up the transtions (vector field matching) across the cuts
    for (int i = 0; i < VH.rows(); i++)
    {
      std::vector<MatrixXi> permMatrices;
      std::vector<int> permIndices;  //in the space #V + #transitions
      //The initial corner gets the identity without any transition
      permMatrices.push_back(MatrixXi::Identity(N, N));
      permIndices.push_back(i);
      
      int beginH = VH(i);
      int currH = beginH;
      
      //reseting to first cut or boundary, if exists
      if (!isBoundary(i))
      {
        // travel throu the start of the vertex and stop once the edge on the cut is found
        do
        {
          if (isHEcut(currH) != 0)
            break;
          currH = nextH(twinH(currH));
        } while(beginH != currH);
      }
      else
      {
        do
        {
          // travel until an edge without a twin is found, i.e., boundary
          if (twinH(currH) == -1)
            break;
          currH = nextH(twinH(currH));
        } while(twinH(currH) != -1);
      }
      
      // set the beginning to the edge on the cut or on the boundary
      beginH = currH;
      
      int currCutVertex = -1;
      do
      {
        int currFace = HF(currH); // face containing the half-edge
        int newCutVertex = -1;
        //find position of the vertex i in the face of the initial mesh
        for (int j = 0; j < 3; j++)
        {
          if (wholeF(currFace, j) == i)
            newCutVertex = cutF(currFace, j);
        }
        
        //currCorner gets the permutations so far
        if (newCutVertex != currCutVertex)
        {
          currCutVertex = newCutVertex;
          for(int i = 0; i < permIndices.size(); i++)
          {
            // place the perumtation matrix in a bigger matrix, we need to know how things are connected along the cut, no?
            for(int j = 0; j < N; j++)
              for(int k = 0; k < N; k++){
                vertexTrans2CutTriplets.emplace_back(N * currCutVertex + j, N * permIndices[i] + k, (double) permMatrices[i](j, k));
                vertexTrans2CutTripletsInteger.emplace_back(N * currCutVertex + j, N * permIndices[i] + k, permMatrices[i](j, k));
              }
          }
        }
        
        //updating the matrices for the next corner
        int nextHalfedge = twinH(prevH(currH));
        //reached a boundary
        if(nextHalfedge == -1)
        {
          currH = nextHalfedge;
          continue;
        }
        
        // constParmMatrices contains all the members of the permutation group
        MatrixXi nextPermMatrix = constParmMatrices[Halfedge2Matching(nextHalfedge) % N];
        //no update needed
        if(isHEcut(nextHalfedge) == 0)
        {
          currH = nextHalfedge;
          continue;
        }
        
        //otherwise, updating matrices with transition
        int nextTransition = Halfedge2TransitionIndices(nextHalfedge);
        //Pe*f + Je
        if(nextTransition > 0)
        {
          for(int j = 0; j < permMatrices.size(); j++)
            permMatrices[j] = nextPermMatrix * permMatrices[j];
          
          //and identity on the fresh transition
          permMatrices.push_back(MatrixXi::Identity(N, N));
          permIndices.push_back(wholeV.rows() + nextTransition - 1);
        }
        // (Pe*(f-Je))  matrix is already inverse since halfedge matching is minused
        else
        {
          //reverse order
          permMatrices.push_back(-MatrixXi::Identity(N, N));
          permIndices.push_back(wholeV.rows() - nextTransition - 1);
          
          for(int j = 0; j < permMatrices.size(); j++)
            permMatrices[j] = nextPermMatrix * permMatrices[j];
        }
        currH = nextHalfedge;
      } while((currH != beginH) && (currH != -1));
      
      //cleaning parmMatrices and permIndices to see if there is a constraint or reveal singularity-from-transition
      std::set<int> cleanPermIndicesSet(permIndices.begin(), permIndices.end());
      std::vector<int> cleanPermIndices(cleanPermIndicesSet.begin(), cleanPermIndicesSet.end());
      std::vector<MatrixXi> cleanPermMatrices(cleanPermIndices.size());
      
      for (int j = 0; j < cleanPermIndices.size(); j++)
      {
        cleanPermMatrices[j] = MatrixXi::Zero(N, N);
        for(int k = 0;k < permIndices.size(); k++)
          if(cleanPermIndices[j] == permIndices[k])
            cleanPermMatrices[j] += permMatrices[k];
        if(cleanPermIndices[j] == i)
          cleanPermMatrices[j] -= MatrixXi::Identity(N, N);
      }
      
      //if not all matrices are zero, there is a constraint
      bool isConstraint = false;
      for(int j = 0; j < cleanPermMatrices.size(); j++)
        if (cleanPermMatrices[j].cwiseAbs().maxCoeff() != 0)
          isConstraint = true;
      
      if((isConstraint) && (!isBoundary(i)))
      {
        for(int j = 0; j < cleanPermMatrices.size(); j++)
        {
          for(int k = 0; k < N; k++)
            for(int l = 0; l < N; l++){
              constTriplets.emplace_back(N * currConst + k, N * cleanPermIndices[j] + l, (double) cleanPermMatrices[j](k, l));
              constTripletsInteger.emplace_back(N * currConst + k, N * cleanPermIndices[j] + l, cleanPermMatrices[j](k, l));
            }
        }
        currConst++;
        intData.constrainedVertices(i) = 1;
      }
    }
    
    vector< Triplet< double > > cleanTriplets;
    vector< Triplet< int > > cleanTripletsInteger;
    
    intData.vertexTrans2CutMat.conservativeResize(N * cutV.rows(), N * (wholeV.rows() + numTransitions));
    intData.vertexTrans2CutMatInteger.conservativeResize(N * cutV.rows(), N * (wholeV.rows() + numTransitions));
    cleanTriplets.clear();
    cleanTripletsInteger.clear();
    for(int i = 0; i < vertexTrans2CutTriplets.size(); i++){
      if(vertexTrans2CutTripletsInteger[i].value() != 0){
        cleanTripletsInteger.push_back(vertexTrans2CutTripletsInteger[i]);
        cleanTriplets.push_back(vertexTrans2CutTriplets[i]);
      }
      // if(std::abs((float)vertexTrans2CutTriplets[i].value())>10e-7)
    }
    intData.vertexTrans2CutMat.setFromTriplets(cleanTriplets.begin(), cleanTriplets.end());
    intData.vertexTrans2CutMatInteger.setFromTriplets(cleanTripletsInteger.begin(), cleanTripletsInteger.end());
    
    //
    
    intData.constraintMat.conservativeResize(N * currConst, N * (wholeV.rows() + numTransitions));
    intData.constraintMatInteger.conservativeResize(N * currConst, N * (wholeV.rows() + numTransitions));
    cleanTriplets.clear();
    cleanTripletsInteger.clear();
    for(int i = 0; i < constTriplets.size(); i++){
      if(constTripletsInteger[i].value() != 0){
        cleanTripletsInteger.push_back(constTripletsInteger[i]);
        cleanTriplets.push_back(constTriplets[i]);
      }
    }
    intData.constraintMat.setFromTriplets(cleanTriplets.begin(), cleanTriplets.end());
    intData.constraintMatInteger.setFromTriplets(cleanTripletsInteger.begin(), cleanTripletsInteger.end());
    
    //doing the integer spanning matrix
    intData.intSpanMat.conservativeResize(n * (wholeV.rows() + numTransitions), n * (wholeV.rows() + numTransitions));
    intData.intSpanMatInteger.conservativeResize(n * (wholeV.rows() + numTransitions), n * (wholeV.rows() + numTransitions));
    vector<Triplet<double> > intSpanMatTriplets;
    vector<Triplet<int> > intSpanMatTripletsInteger;
    for (int i=0;i<n*numTransitions;i+=n){
      for(int k = 0; k < n; k++)
        for(int l = 0; l < n; l++){
          if (intFunc(k,l)!=0){
            intSpanMatTriplets.emplace_back(n * wholeV.rows()+i+k, n * wholeV.rows()+i+l, (double)intFunc(k,l));
            intSpanMatTripletsInteger.emplace_back(n * wholeV.rows()+i+k, n * wholeV.rows()+i+l, intFunc(k,l));
          }
        }
    }
    for (int i=0;i<n * wholeV.rows();i++){
      intSpanMatTriplets.emplace_back(i,i,1.0);
      intSpanMatTripletsInteger.emplace_back(i,i,1);
    }
    
    intData.intSpanMat.setFromTriplets(intSpanMatTriplets.begin(), intSpanMatTriplets.end());
    intData.intSpanMatInteger.setFromTriplets(intSpanMatTripletsInteger.begin(), intSpanMatTripletsInteger.end());
    
    //filtering out barycentric symmetry, including sign symmetry. The parameterization should always only include d dof for the surface
    //TODO: this assumes d divides N!
    intData.symmMat.conservativeResize(N * (wholeV.rows() + numTransitions), n * (wholeV.rows() + numTransitions));
    intData.symmMatInteger.conservativeResize(N * (wholeV.rows() + numTransitions), n * (wholeV.rows() + numTransitions));
    vector<Triplet<double> > symmMatTriplets;
    vector<Triplet<int> > symmMatTripletsInteger;
    for(int i = 0; i < N*(wholeV.rows() + numTransitions); i +=N)
      for(int k = 0; k < N; k++)
        for(int l = 0; l < n; l++){
          if (symmFunc(k,l)!=0){
            symmMatTriplets.emplace_back(i + k, i*n/N + l, (double)symmFunc(k,l));
            symmMatTripletsInteger.emplace_back(i + k, i*n/N + l, symmFunc(k,l));
          }
        }
    
    intData.symmMat.setFromTriplets(symmMatTriplets.begin(), symmMatTriplets.end());
    intData.symmMatInteger.setFromTriplets(symmMatTripletsInteger.begin(), symmMatTripletsInteger.end());
    
    //integer variables are per single "d" packet, and the rounding is done for the N functions with projection over symmFunc
    intData.integerVars.conservativeResize(numTransitions);
    intData.integerVars.setZero();
    for(int i = 0; i < numTransitions; i++)
      intData.integerVars(i) = wholeV.rows() + i;
    
    //fixed values
    intData.fixedIndices.resize(intData.n);
    if (isSingular.sum()==0){  //no inner singular vertices; vertex 0 is set to (0....0)
      for (int j=0;j<intData.n;j++)
        intData.fixedIndices(j)=j;
    }else {  //fixing first singularity to (0.5,....0.5)
      int firstSing;
      for (firstSing=0;firstSing<isSingular.size();firstSing++)
        if (isSingular(firstSing))
          break;
      //firstSing=0; //like before
      for (int j=0;j<intData.n;j++)
        intData.fixedIndices(j)=intData.n*firstSing+j;
    }
    
    //creating list of singular corners and singular integer matrix
    VectorXi singularIndices(intData.n * isSingular.sum());
    int counter=0;
    for (int i=0;i<isSingular.size();i++){
      if (isSingular(i))
        for (int j=0;j<intData.n;j++)
          singularIndices(counter++)=intData.n*i+j;
    }
    
    //doing the integer spanning matrix
    intData.singIntSpanMat.conservativeResize(n * (wholeV.rows() + numTransitions), n * (wholeV.rows() + numTransitions));
    intData.singIntSpanMatInteger.conservativeResize(n * (wholeV.rows() + numTransitions), n * (wholeV.rows() + numTransitions));
    vector<Triplet<double> > singIntSpanMatTriplets;
    vector<Triplet<int> > singIntSpanMatTripletsInteger;
    for (int i=0;i<isSingular.size();i++){
      if (!isSingular(i)){
        for (int j=0;j<n;j++){
          singIntSpanMatTriplets.emplace_back(n*i+j,n*i+j,1.0);
          singIntSpanMatTripletsInteger.emplace_back(n*i+j,n*i+j,1);
        }
      } else {
        for(int k = 0; k < n; k++)
          for(int l = 0; l < n; l++){
            if (intFunc(k,l)!=0){
              singIntSpanMatTriplets.emplace_back(n*i+k, n*i+l, (double)intFunc(k,l));
              singIntSpanMatTripletsInteger.emplace_back(n*i+k, n*i+l, intFunc(k,l));
            }
          }
      }
    }
    
    for (int i=n * wholeV.rows() ; i<n*(wholeV.rows()+numTransitions);i++){
      singIntSpanMatTriplets.emplace_back(i,i,1.0);
      singIntSpanMatTripletsInteger.emplace_back(i,i,1);
    }
    
    intData.singIntSpanMat.setFromTriplets(singIntSpanMatTriplets.begin(), singIntSpanMatTriplets.end());
    intData.singIntSpanMatInteger.setFromTriplets(singIntSpanMatTripletsInteger.begin(), singIntSpanMatTripletsInteger.end());
    
    intData.singularIndices=singularIndices;
    intData.fixedValues.resize(intData.n);
    intData.fixedValues.setConstant(0);
    
  }
  }

#endif



#ifndef NumericalKinematicConstrainedFitter_H
#define NumericalKinematicConstrainedFitter_H

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoTauTag/KinematicTau/interface/MultiTrackNumericalKinematicConstraint.h"
#include "RecoVertex/VertexTools/interface/LinearizationPointFinder.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexUpdator.h"
#include "RecoVertex/KinematicFit/interface/VertexKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/ConstrainedTreeBuilder.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class NumericalKinematicConstrainedFitter{

public:

 NumericalKinematicConstrainedFitter();
 NumericalKinematicConstrainedFitter(const LinearizationPointFinder& fnd);
  
 ~NumericalKinematicConstrainedFitter();

 void setParameters(const edm::ParameterSet& pSet);

 RefCountedKinematicTree fit(std::vector<RefCountedKinematicParticle> part) {
   return fit(part, 0, 0);
 }

 RefCountedKinematicTree fit(std::vector<RefCountedKinematicParticle> part,
			     MultiTrackNumericalKinematicConstraint * cs) {
   return fit(part, cs, 0);
 };
 
 RefCountedKinematicTree fit(std::vector<RefCountedKinematicParticle> part,
			     MultiTrackNumericalKinematicConstraint * cs,
			     GlobalPoint * pt);

//return the number of iterations
 int getNit() const;
//return the value of the constraint equation
 float getCSum() const;

private:

 void defaultParameters();

 float theMaxDelta; //maximum (delta parameter)^2/(sigma parameter)^2 per iteration for convergence
 int theMaxStep; 				       
 float theMaxReducedChiSq; //max of initial (after 2 iterations) chisq/dof value
 float theMinChiSqImprovement; //minimum required improvement in chisq to avoid fit termination for cases exceeding theMaxReducedChiSq
 LinearizationPointFinder * finder;				       
 KinematicConstrainedVertexUpdator * updator;
 VertexKinematicConstraint * vCons;
 ConstrainedTreeBuilder * tBuilder;
 int iterations;
 float csum;
};

#endif

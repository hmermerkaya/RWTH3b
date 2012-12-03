#include "RecoTauTag/KinematicTau/interface/NumericalKinematicConstrainedFitter.h"
#include "RecoVertex/KinematicFit/interface/InputSort.h"
#include "RecoVertex/LinearizationPointFinders/interface/DefaultLinearizationPointFinder.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertexFactory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"



NumericalKinematicConstrainedFitter::NumericalKinematicConstrainedFitter()
{
 finder = new DefaultLinearizationPointFinder();
 tBuilder = new ConstrainedTreeBuilder;
 defaultParameters();
 iterations = -1;
 csum = -1000.0;
}

NumericalKinematicConstrainedFitter::NumericalKinematicConstrainedFitter(const LinearizationPointFinder& fnd)
{
 finder = fnd.clone();
 tBuilder = new ConstrainedTreeBuilder;
 defaultParameters();
 iterations = -1;
 csum = -1000.0;
}

NumericalKinematicConstrainedFitter::~NumericalKinematicConstrainedFitter()
{
 delete finder;
 delete tBuilder;
}

void NumericalKinematicConstrainedFitter::setParameters(const edm::ParameterSet& pSet)
{
  theMaxDelta = pSet.getParameter<double>("maxDelta");
  theMaxStep = pSet.getParameter<int>("maxNbrOfIterations");
  theMaxReducedChiSq = pSet.getParameter<double>("maxReducedChiSq");
  theMinChiSqImprovement = pSet.getParameter<double>("minChiSqImprovement");
}

void NumericalKinematicConstrainedFitter::defaultParameters()
{
  theMaxDelta = 0.01;
  theMaxStep = 1000;
  theMaxReducedChiSq = 225.;
  theMinChiSqImprovement = 50.;
}

RefCountedKinematicTree NumericalKinematicConstrainedFitter::fit(std::vector<RefCountedKinematicParticle> part,
                                                             MultiTrackNumericalKinematicConstraint * cs,
                                                             GlobalPoint * pt){

  std::cout << "RefCountedKinematicTree NumericalKinematicConstrainedFitter::fit Start" << std::endl;
  if(part.size()<2) throw VertexException("NumericalKinematicConstrainedFitter::input states are less than 2");
  
  //sorting out the input particles
  InputSort iSort;
  std::pair<std::vector<RefCountedKinematicParticle>, std::vector<FreeTrajectoryState> > input = iSort.sort(part);
  const std::vector<RefCountedKinematicParticle> & particles  = input.first;
  const std::vector<FreeTrajectoryState> & fStates = input.second;
  
  // linearization point
  // (only compute it using the linearization point finder if no point was passed to the fit function):
  GlobalPoint linPoint;
  if (pt!=0) {
    linPoint  = *pt;
  }
  else {
    linPoint = finder->getLinearizationPoint(fStates);
  }
  GlobalPoint lPoint  = linPoint;
  RefCountedKinematicVertex rVtx;
  AlgebraicMatrix refCCov;

  int vSize = particles.size();
  AlgebraicVector inPar(3 + 7*vSize,0);
  AlgebraicVector finPar(3 + 7*vSize,0);
  AlgebraicMatrix inCov(3 + 7*vSize,3 + 7*vSize,0);
  
  //making initial vector of parameters and initial particle-related covariance
  int nSt = 0;
  std::vector<KinematicState> inStates;
  std::vector<KinematicState> lStates = inStates;
  for(std::vector<RefCountedKinematicParticle>::const_iterator i = particles.begin(); i!=particles.end(); i++){
    KinematicState state = (*i)->stateAtPoint(linPoint);
    if (!state.isValid()) {
      LogDebug("KinematicConstrainedVertexFitter")
	<< "State is invalid at point: "<<linPoint<<std::endl;
      return ReferenceCountingPointer<KinematicTree>(new KinematicTree());
    }
    AlgebraicVector prPar = asHepVector<7>(state.kinematicParameters().vector());
    for(int j = 1; j<8; j++){inPar(3 + 7*nSt + j) = prPar(j);}
    AlgebraicSymMatrix l_cov  = asHepMatrix<7>(state.kinematicParametersError().matrix());
    inCov.sub(4 + 7*nSt,4 + 7*nSt ,l_cov);
    inStates.push_back(state);
    ++nSt;
  }
  
  //iterarions over the updator: each time updated parameters
  //are taken as new linearization point
  cs->ConfigureIntialState(inStates,linPoint);
  int nit(0);
  double chi2(1e6),delta(0);
  for(nit=0;nit<=theMaxStep;nit++){
    bool passed=cs->ApplyLagrangianConstraints(chi2,delta);
    /*if (!passed || (nit==theMaxStep && delta>=4.0*theMaxDelta)) {
      std::cout << "RefCountedKinematicTree NumericalKinematicConstrainedFitter::fit end" << std::endl;
      return ReferenceCountingPointer<KinematicTree>(new KinematicTree());
      } */
    std::cout << "Loop Iteration " << nit << " chi2 " << chi2 << " delta " << delta << std::endl; 
    if(cs->isConverged(chi2,delta))break;
  }
  std::pair< std::pair< std::vector<KinematicState>, AlgebraicMatrix >,RefCountedKinematicVertex> lRes =cs->ConvertStateToParameters(inStates,linPoint);  
  const std::vector<KinematicState> &newStates = lRes.first.first;
  lStates = newStates;
  rVtx = lRes.second;
  refCCov = lRes.first.second;
  iterations = nit;
  csum = cs->numberOfEquations();
  std::cout << "RefCountedKinematicTree NumericalKinematicConstrainedFitter::fit end 1" << std::endl;
  return  tBuilder->buildTree(particles, lStates, rVtx, refCCov);
}

int NumericalKinematicConstrainedFitter::getNit() const {
    return iterations;
}

float NumericalKinematicConstrainedFitter::getCSum() const {
    return csum;
}

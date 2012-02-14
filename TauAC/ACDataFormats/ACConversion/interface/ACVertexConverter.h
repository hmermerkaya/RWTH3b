#ifndef ACVertexConverter_h
#define ACVertexConverter_h
/**
 Vertex helper class to convert from CMSSW objects to ACCommon objects. In general, helper classes inherit from the corresponding tree object, expand the constructor. They do not own additional private members.

 @author Lars Perchalla & Philip Sauerland
 @date 2011
 */

#include "../../ACGeneral/interface/ACVertex.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"

class ACVertexConverter : public ACVertex {
public:
    virtual std::string classname() const { return "ACVertexConverter"; }
    /// constructor from reco::Vertex to initiate an ACVertex
    ACVertexConverter(const reco::Vertex& vertex);
    virtual ~ACVertexConverter() {};
};

#endif

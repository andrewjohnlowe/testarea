#ifndef JETSHAPES_H
#define JETSHAPES_H

// #include "JetEvent/Jet.h"
// #include "JetEvent/AliPicoTrack.h"
// #include "JetEvent/JetAssociationGeneric.h" // To get subjets
// #include "JetUtils/TClonesArrayHelper.h"
// #include "JetTagInfo/TruthInfo.h" // For MC truth tagging only
// #include "JetUtils/JetSorters.h" // For sorting (sub)jets

#include "AliEmcalJet.h"

#include "AliVCluster.h"
#include "AliPicoTrack.h"
#include "AliPicoTrack.h"
#include "AliVParticle.h"

#include "CLHEP/Vector/TwoVector.h" // Used in some jet shapes
#include "CLHEP/Vector/ThreeVector.h" // Common to most variables
#include "CLHEP/Units/PhysicalConstants.h" // Used for angularities (tau)
#include "CLHEP/Vector/LorentzVector.h"

#include "/home/andy/testarea/Hep2Matrix.h" // Pull vectors and planar flow
#include "/home/andy/testarea/Numeric.h" // Common numeric stuff
#include "/home/andy/testarea/JetFlags.h" // Common jet reco flags
#include <vector> // For holding jets temporarily
#include <utility> // For pairs returned by tau functor(s)
#include <string> // To determine jet inputs and algorithm
#include <map> // An elegant kludge to workaround AliPicoTrack having no copy CTOR beyond shallow copy 
 
// debugging macros so we can pin down message provenance at a glance
#include <iostream>
#define DEBUG(x)							\
  std::cout << "DEBUG: " << "\t" << #x << "\t" << x << "\t"		\
  << __FUNCTION__ << "\t" << __FILE__ << "\t" << __LINE__ << std::endl;
#define DEBUG(x)

using namespace Numeric;

//typedef AliPicoTrack AliPicoTrack;

class JetMisc {
public:
  JetMisc(const AliEmcalJet*& j): _j(j) {
  //_j(j), _association(0), _name(""), _subjetColl(0), _numSubjets(0) {
    TLorentzVector *hlv = new TLorentzVector;
    j->GetMom(*hlv);
    _theHLV = new CLHEP::HepLorentzVector(hlv->X(), hlv->Y(), hlv->Z(), hlv->T());
    // if(j != 0 && (theHLV.restMass2() <= m2min || theHLV.findBoostToCM().mag2() >= b2max)) { // !ZMxpvTachyonic
    //   i.setX(NaN);
    //   i.setY(NaN);
    //   i.setZ(NaN);
    //   i.setT(NaN);
    //   j = 0; // signal that we are to not boost, pass back NaN HLV instead
    // }
    // _association = _j->getAssociation<JetAssociationGeneric<TClonesArray> >("CoreSubJets"); DEBUG(_association);
    // if(_association) {
    //   _name = _association->name(); DEBUG(_name);
    //   _subjetColl = _association->get(); DEBUG(_subjetColl);
    //   if(_subjetColl) _numSubjets = _subjetColl->size(); DEBUG(_numSubjets);
    // }
  }
  int/* jStr_indexJ */index() {}
  int/* jStr_jetMiscJ_*/numConstituents() { // PUBLISH_JET_SHAPE
    return _j->Nch();
  }
  double/* jStr_jetMiscJ_*/jetM() { // PUBLISH_JET_SHAPE
    return _j->M();
  }
  double/* jStr_jetMiscJ_*/jetMt() { // PUBLISH_JET_SHAPE
    //return _j->hlv().mt();
  return _theHLV->mt();
  }
  double/* jStr_jetMiscJ_*/jetE() { // PUBLISH_JET_SHAPE
    //return _j->e();
  return _theHLV->e();
  }
  double/* jStr_jetMiscJ_*/jetP() { // PUBLISH_JET_SHAPE
    //return _j->p();
  return _theHLV->vect().mag();
  }
  double/* jStr_jetMiscJ_*/jetEt() { // PUBLISH_JET_SHAPE
  //return _j->et();
  return _theHLV->et();
  }
  double/* jStr_jetMiscJ_*/jetPt() { // PUBLISH_JET_SHAPE
    //return _j->pt();
  return _theHLV->perp();
  }
  double/* jStr_jetMiscJ_*/jetPhi() { // PUBLISH_JET_SHAPE
    //return _j->phi();
  return _theHLV->phi();
  }
  double/* jStr_jetMiscJ_*/jetEta() { // PUBLISH_JET_SHAPE
    //return _j->eta();
  return _theHLV->eta();
  }
  double/* jStr_jetMiscJ_*/jetRapidity() { // PUBLISH_JET_SHAPE
    //return _j->rapidity();
  return _theHLV->rapidity();
  }
  // http://arxiv.org/pdf/1101.1335v1
  double/* jStr_jetMiscJ_*/xJ() { // PUBLISH_JET_SHAPE
    double m = jetM();
    double pt = jetPt();
    return Divide(m * m, pt * pt);
  }
  // http://arxiv.org/pdf/1101.1335v1
  double/* jStr_jetMiscJ_*/gamma() { // PUBLISH_JET_SHAPE
    return Sqrt(Divide(1., xJ()) + 1.);
  }
  double/* jStr_jetMiscJ_*/R() { // PUBLISH_JET_SHAPE
    RadialParameter R;
    return 0.1 * R(_j);
  }
  double/* jStr_jetMiscJ_*/Cbar() { // PUBLISH_JET_SHAPE
    return Divide(jetM(), jetPt() * R());
  }
  const TClonesArray* subjets() {
    //return _subjetColl.get();
    return 0;//_subjetColl;
  }
  int/* jStr_jetMiscJ_*/numSubjets() { // PUBLISH_JET_SHAPE
    return -1;//_numSubjets;
  }
  std::string name() {
    return "";//_name;
  }

~JetMisc() {
  delete _theHLV;
  }

private:
  const AliEmcalJet*& _j;
  CLHEP::HepLorentzVector* _theHLV;
  //const TObject* _association;
  std::string _name;
  //const TClonesArray* _subjetColl;
  //unsigned int _numSubjets;
};

// arXiv: 1001.5027v3 [hep-ph] // changed pt --> et
class Pull {
public:
  Pull(const AliEmcalJet*& j, double kappa = 0.3): _jet(j), _kappa(kappa), _G(0), _halfG(0), _QG(0), _sumPxNmag(0.), _sumPmag(0.) {
    TLorentzVector *hlv = new TLorentzVector;
    j->GetMom(*hlv);
    _theHLV = new CLHEP::HepLorentzVector(hlv->X(), hlv->Y(), hlv->Z(), hlv->T());
    _sumPt = 0.;
    _sumPt2 = 0.;
    _sumQ = 0.;
    _sumBQ = 0.;
    _sumBPt = 0.;
    _momPhi2 = 0.;
    _momY2 = 0.;
    if(_jet->Nch() <= 1) {
      _G = NaN;
      _halfG = NaN;
      _QG = NaN;
      _sumPxNmag = NaN;
      _sumPmag = NaN;
      _momPhi2 = NaN;
      _momY2 = NaN;
      _t = CLHEP::Hep2Vector(NaN, NaN);
      _C[0][0] = _C[0][1] = _C[1][0] = _C[1][1] = NaN;
    }
    else {
      _t = CLHEP::Hep2Vector(0, 0);
      _C = Hep2Matrix(0);
    }
  }
  
  ~Pull() {
    delete _theHLV;
  }

  void operator()(const AliPicoTrack i) {
    term(i);
  }

  void operator()(const AliPicoTrack* i) {
    term(*i);
  }
  
  const AliEmcalJet* jet() const { return _jet; } 
  const CLHEP::Hep2Vector* t() const { return &_t; } // PUBLISH_JET_SHAPE
  const Hep2Matrix* C() const { return &_C; } // PUBLISH_JET_SHAPE

  // Jason Gallicchio http://frank.harvard.edu/~jason/physics/qvsg/qvsg_beamer_v02.pdf
  double/* jStr_pullJ_*/det() const { return _C.determinant(); } // PUBLISH_JET_SHAPE

  // Jason Gallicchio http://frank.harvard.edu/~jason/physics/qvsg/qvsg_beamer_v02.pdf
  double/* jStr_pullJ_*/ratio() const { // PUBLISH_JET_SHAPE
    double e1 = _C.eigenvalues().first;
    double e2 = _C.eigenvalues().second;
    return Divide(e2, e1);
  }

  // Jason Gallicchio http://frank.harvard.edu/~jason/physics/qvsg/qvsg_beamer_v02.pdf
  double/* jStr_pullJ_*/pullPf() const { // PUBLISH_JET_SHAPE
    double e1 = _C.eigenvalues().first;
    double e2 = _C.eigenvalues().second;
    return Divide(4.*e1*e2, (e1 + e2)*(e1 + e2));
  }

  // http://en.wikipedia.org/wiki/Angular_eccentricity
  double/* jStr_pullJ_*/angularEccentricity() const { // PUBLISH_JET_SHAPE
    double e1 = _C.eigenvalues().first;
    double e2 = _C.eigenvalues().second;
    return Arccos(Divide(e2, e1));
  }

  // Jason Gallicchio http://frank.harvard.edu/~jason/physics/qvsg/qvsg_beamer_v02.pdf
  double/* jStr_pullJ_*/orientation() const { // PUBLISH_JET_SHAPE
    return (-1.* _C.eigenvectors().first).phi();
  }

  double/* jStr_pullJ_*/girth() const { return _G; } // PUBLISH_JET_SHAPE

  // arxiv.org/pdf/1211.703 (eqn 8.5)
  double/* jStr_pullJ_*/halfgirth() const { return _halfG; } // PUBLISH_JET_SHAPE

  // My own invention; try to get insight on charge distribution inside jet
  double/* jStr_pullJ_*/Qgirth() const { return _QG; } // PUBLISH_JET_SHAPE

  // http://www-cdf.fnal.gov/physics/new/hdg/results/vhqqbb_oct09/AllHadronicHiggsSearch/Jet_Shapes_Studies.html
  double/* jStr_pullJ_*/phiMom() const { return sqrt(_momPhi2); } // PUBLISH_JET_SHAPE
  double/* jStr_pullJ_*/yMom() const { return sqrt(_momY2); } // PUBLISH_JET_SHAPE

  double/* jStr_pullJ_*/Cbar() const { return Divide(_theHLV->m(), _theHLV->perp() * _G); } // PUBLISH_JET_SHAPE

  // size
  double/* jStr_pullJ_*/g() const { // PUBLISH_JET_SHAPE
    double e1 = _C.eigenvalues().first;
    double e2 = _C.eigenvalues().second;
    return Sqrt((e1 * e1) + (e2 * e2));
  }

  // != linear eccentricity (c)
  double/* jStr_pullJ_*/e() const { // PUBLISH_JET_SHAPE
    double e1 = _C.eigenvalues().first;
    double e2 = _C.eigenvalues().second;
    return Sqrt(Divide((e1 * e1) - (e2 * e2), e1));
  }

  // Broadening, http://jets.physics.harvard.edu/qvg/
  double/* jStr_pullJ_*/B() const { // PUBLISH_JET_SHAPE
    return Divide(_sumPxNmag, _sumPmag);
  }

  double/* jStr_pullJ_*/logB() const { // PUBLISH_JET_SHAPE
    return Ln(Divide(_sumPxNmag, _sumPmag));
  }

  double/* jStr_pullJ_*/pullTheta() const { return _t.phi(); } // PUBLISH_JET_SHAPE
  double/* jStr_pullJ_*/pullMag() const { return _t.mag(); } // PUBLISH_JET_SHAPE
  double/* jStr_pullJ_*/JetW() const { return sqrt(_sumPt2/_jet->Nch()) / _jet->Pt(); }
  double/* jStr_pullJ_*/JetQ() const { return _sumQ * (1./pow(_jet->Pt(),_kappa)); }
  // A study of b-jet charge structure, Robert Astalos:
  double/* jStr_pullJ_*/JetB() const { return _sumBQ / _sumBPt; }


private:
  void term(const AliPicoTrack i) {
    _sumQ += ((double)i.Charge()) * pow(i.Pt(),_kappa); DEBUG(i.Charge()); DEBUG(i.Pt()); DEBUG(pow(i.Pt(),_kappa));
    _sumPt += i.Pt();
    _sumPt2 += i.Pt() * i.Pt();

    CLHEP::Hep3Vector P_i = CLHEP::Hep3Vector(i.Px(), i.Py(), i.Pz());
    CLHEP::Hep3Vector JJ = CLHEP::Hep3Vector(_theHLV->px(), _theHLV->py(), _theHLV->pz());

    _sumBPt += pow(fabs(JJ.dot(P_i)), _kappa);
    _sumBQ += ((double)i.Charge()) * pow(fabs(JJ.dot(P_i)), _kappa);

    if(_jet->Nch() <= 1) return;
    CLHEP::Hep2Vector ri;
    Hep2Matrix Vi(0);
    
    CLHEP::Hep3Vector N = (CLHEP::Hep3Vector(_theHLV->px(), _theHLV->py(), _theHLV->pz())).unit();
    CLHEP::Hep3Vector P = CLHEP::Hep3Vector(i.Px(), i.Py(), i.Pz());
    double PxNmag = P.cross(N).mag();
    double Pmag = P.mag();

    _sumPxNmag += PxNmag;
    _sumPmag += Pmag;

    CLHEP::HepLorentzVector t = CLHEP::HepLorentzVector(i.Px(), i.Py(), i.Pz(), i.E());

    double dy = t.rapidity() - _theHLV->rapidity();
    double dphi = deltaPhi(t.phi(), _theHLV->phi());
    ri = CLHEP::Hep2Vector(dy, dphi);

    Vi[0][0] = ri.x() * ri.x();
    Vi[0][1] = ri.x() * ri.y();
    Vi[1][0] = Vi[0][1];
    Vi[1][1] = ri.y() * ri.y();

    double et_J = _theHLV->et();

  // http://www-cdf.fnal.gov/physics/new/hdg/results/vhqqbb_oct09/AllHadronicHiggsSearch/Jet_Shapes_Studies.html
    _momY2 += Divide((t.et() * dy*dy), et_J);
    _momPhi2 += Divide((t.et() * dphi*dphi), et_J);

    _G += Divide((t.et() * ri.mag()), et_J); // arXiv:1012.2077
    _halfG += Divide((t.et() * sqrt(ri.mag())), et_J); // arXiv:1211.7038
    _QG += Divide(((double)i.Charge()) * (t.et() * sqrt(ri.mag())), et_J);
    _t += Divide((t.et() * ri.mag()), et_J) * ri;
    _C += Divide((t.et() * ri.mag()), et_J) * Vi;
  }
  
  CLHEP::HepLorentzVector* _theHLV;
  PhiCorr phiCorr;
  DeltaPhi deltaPhi;
  const AliEmcalJet* _jet;
  CLHEP::Hep2Vector _t;
  Hep2Matrix _C;
  double _G, _halfG, _QG;
  double _sumPxNmag, _sumPmag;
  double _sumQ;
  double _sumPt;
  double _sumBQ;
  double _sumBPt;
  double _sumPt2;
  double _momPhi2;
  double _momY2;
  double _kappa;
};

// // You can't do any of this until you can extract subjets from a jet ala BDRS
// // arXiv:1006.1650v1 [hep-ph] (Unburied Higgs)
// class UnburiedHiggs {
// public:
//   UnburiedHiggs(const AliEmcalJet*& j):
//   _jet(j), _sortedSubJets(0), _etsum(0.), _ptsum(0.) {
//     JetMisc misc(_jet);
//     _subjets = misc.subjets();
//     _numSubjets = misc.numSubjets();
//     _subjetsV.clear();
//     if(_subjets) {
//       TClonesArrayHelper::TClonesArray_t sortedSubJets(_subjets->begin(), _subjets->end());
//       TClonesArrayHelper::sort_jets(sortedSubJets, JetSorters::sortJetByEtDown());
//       _sortedSubJets = sortedSubJets;
//       TClonesArrayHelper::TClonesArray_t::iterator firstSubJet = _sortedSubJets.begin();
//       TClonesArrayHelper::TClonesArray_t::iterator lastSubJet  = _sortedSubJets.end();
//       for(;firstSubJet != lastSubJet; ++firstSubJet) {
// 	if((*firstSubJet) == 0) continue;
// 	const AliEmcalJet* aSubJet = (*firstSubJet)->clone(true, true);
// 	if(aSubJet == 0) continue;
// 	_subjetsV.push_back(aSubJet);
// 	_etsum += aSubJet->et();
// 	_ptsum += aSubJet->pt();
//       }
//     }
//     _firstSubJet = _numSubjets >= 1? _subjetsV[0]: 0;
//     _secondSubJet = _numSubjets >= 2? _subjetsV[1]: 0; 
//     _thirdSubJet = _numSubjets >= 3? _subjetsV[2]: 0;
//   }

//   ~UnburiedHiggs() {
//     for(unsigned int i = 0; i != _subjetsV.size(); i++) delete _subjetsV[i];
//     _subjetsV.clear();
//   }

//   double/* jStr_ubhJ_*/z() const { // PUBLISH_JET_SHAPE, min(et_i)/et_J
//     double result = NaN;
//     if(!_subjetsV.empty()) {
//       double minEt = (_subjetsV[_subjetsV.size() - 1])->et();
//       double et_J = _jet->et();
//       result = Divide(minEt, et_J);
//     }
//     return result;
//   }

//   double/* jStr_ubhJ_*/z2() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 2) {
//       double et1 = _firstSubJet->et();
//       double et2 = _secondSubJet->et();
//       result = Divide(std::min(et1, et2), et1 + et2);
//     }
//     return result;
//   }

//   // arXiv:1101.1335v1 [hep-ph]
//   double/* jStr_ubhJ_*/a1() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 1) result = Divide(_firstSubJet->m(), _jet->m());
//     return result;
//   }
  
//   // arXiv:1101.1335v1 [hep-ph]
//   double/* jStr_ubhJ_*/a2() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 2) result = Divide(_secondSubJet->m(), _jet->m());
//     return result;
//   }

//   // Based on above
//   double/* jStr_ubhJ_*/a3() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 3) result = Divide(_thirdSubJet->m(), _jet->m());
//     return result;
//   }
  
//   // http://arxiv.org/pdf/1010.3698 
//   double/* jStr_ubhJ_*/meanpt() const { 
//     double result = NaN;
//     if(_numSubjets > 0) result = Divide(_ptsum, (double)_numSubjets);
//     return result;
//   }
  
//   double/* jStr_ubhJ_*/meanet() const { 
//     double result = NaN;
//     if(_numSubjets > 0) result = Divide(_etsum, (double)_numSubjets);
//     return result;
//   }

//   // arXiv:1006.1650v1 [hep-ph] (Unburied Higgs)
//   double/* jStr_ubhJ_*/mbar() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 2) {
//       double m1 = _firstSubJet->m();
//       double m2 = _secondSubJet->m();
//       result = Divide(m1 + m2, 2.);
//     }
//     return result;
//   }

//   // arXiv:1006.1650v1 [hep-ph] (Unburied Higgs)
//   double/* jStr_ubhJ_*/massDemocracy() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 2) {
//       double m1 = _firstSubJet->m();
//       double m2 = _secondSubJet->m();
//       result = std::min(Divide(m1, m2), Divide(m2, m1));
//     }
//     return result;
//   }

//   // https://indico.cern.ch/getFile.py/access?contribId=10&resId=0&materialId=slides&confId=150186
//   double/* jStr_ubhJ_*/fE1() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 1) {
//       double e1 = _firstSubJet->e();
//       result = Divide(e1, _jet->e());
//     }
//     return result;
//   }

//   double/* jStr_ubhJ_*/fE2() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 2) {
//       double e2 = _secondSubJet->e();
//       result = Divide(e2, _jet->e());
//     }
//     return result;
//   }

//   double/* jStr_ubhJ_*/fE3() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 3) {
//       double e3 = _thirdSubJet->e();
//       result = Divide(e3, _jet->e());
//     }
//     return result;
//   }

//   double/* jStr_ubhJ_*/fET1() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 1) {
//       double et1 = _firstSubJet->et();
//       result = Divide(et1, _jet->et());
//     }
//     return result;
//   }

//   double/* jStr_ubhJ_*/fET2() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 2) {
//       double et2 = _secondSubJet->et();
//       result = Divide(et2, _jet->et());
//     }
//     return result;
//   }

//   double/* jStr_ubhJ_*/fET3() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 3) {
//       double et3 = _thirdSubJet->et();
//       result = Divide(et3, _jet->et());
//     }
//     return result;
//   }

//   // arXiv: 1101.1628v2 [hep-ex] and other CMS SUSY notes; unboosted
//   // versions; not in jet CM frame
//   double/* jStr_ubhJ_*/Alpha() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 2) {
//       HepLorentzVector hlv1 = _firstSubJet->hlv();
//       HepLorentzVector hlv2 = _secondSubJet->hlv();
//       double m12 = (hlv1 + hlv2).m();
//       double et2 = _secondSubJet->et();
//       result = Divide(et2, m12);
//     }
//     return result;
//   }

//   double/* jStr_ubhJ_*/AlphaT() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 2) {
//       HepLorentzVector hlv1 = _firstSubJet->hlv();
//       HepLorentzVector hlv2 = _secondSubJet->hlv();
//       double mt12 = (hlv1 + hlv2).mt();
//       double et2 = _secondSubJet->et();
//       result = Divide(et2, mt12);
//     }
//     return result;
//   }

//   // arXiv:1006.1650v1 [hep-ph] (Unburied Higgs) // changed pt --> et
//   double/* jStr_ubhJ_*/betaflow(double etmin = 0.) const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 3) {
//       double et1 = _firstSubJet->et();
//       double et2 = _secondSubJet->et();
//       double et3 = _thirdSubJet->et();
//       double Q = Divide(et3, (et1 + et2));
//       result = etmin == 0.? Q: // the default
// 	et3 < etmin? 0.: Q; // suppress very soft radiation; suitable etmin: 1 GeV, 5 GeV (paper) 
//     }
//     return result;
//   }

//   // https://indico.cern.ch/getFile.py/access?contribId=6&resId=0&materialId=slides&confId=160270
//   double/* jStr_ubhJ_*/y23(double etmin = 0.) const { // PUBLISH_JET_SHAPE
//     double b = betaflow(etmin);
//     return b * b;
//   }

//   double/* jStr_ubhJ_*/lny23(double etmin = 0.) const { // PUBLISH_JET_SHAPE
//     double b = betaflow(etmin);
//     return Ln(b * b);
//   }

//   double/* jStr_ubhJ_*/subjetAsymmetry() const { // PUBLISH_JET_SHAPE
//     double result = NaN;
//     if(_numSubjets >= 2) {
//       double et1 = _firstSubJet->et();
//       double et2 = _secondSubJet->et();
//       result = Divide(et1 - et2, et1 + et2);
//     }
//     return result;
//   }

// private:
//   const AliEmcalJet* _jet;
//   const TClonesArray* _subjets;
//   int _numSubjets;
//   TClonesArrayHelper::TClonesArray_t _sortedSubJets;
//   double _etsum, _ptsum;
//   std::vector<const Jet*> _subjetsV;
//   const Jet* _firstSubJet;
//   const Jet* _secondSubJet;
//   const Jet* _thirdSubJet;
// };

// // You can't do this either until you've cracked BDRS
// class JetDipolarity {
// public:
//   JetDipolarity(const AliEmcalJet*& j):
//   _jet(j), _sortedSubJets(0) {
//     JetMisc misc(_jet);
//     _subjets = misc.subjets();
//     _numSubjets = misc.numSubjets();
//     _subjetsV.clear();
//     if(_subjets) {
//       TClonesArrayHelper::TClonesArray_t sortedSubJets(_subjets->begin(), _subjets->end());
//       TClonesArrayHelper::sort_jets(sortedSubJets, JetSorters::sortJetByEtDown());
//       _sortedSubJets = sortedSubJets;
//       TClonesArrayHelper::TClonesArray_t::iterator firstSubJet = _sortedSubJets.begin();
//       TClonesArrayHelper::TClonesArray_t::iterator lastSubJet  = _sortedSubJets.end();
//       for(unsigned int i = 0; firstSubJet != lastSubJet && i != 2; ++firstSubJet, i++) {
// 	if((*firstSubJet) == 0) continue;
// 	const AliEmcalJet* aSubJet = (*firstSubJet)->clone(true, true);
// 	if(aSubJet == 0) continue;
// 	_subjetsV.push_back(aSubJet);
//       }
//     }
//     _firstSubJet = _numSubjets >= 1? _subjetsV[0]: 0;
//     _secondSubJet = _numSubjets >= 2? _subjetsV[1]: 0;
//     _sum = NaN;
//     _jetSep2 = NaN;
//     if(_numSubjets >= 2) {
//       _rapidity1 = _firstSubJet->rapidity(); // x1
//       _rapidity2 = _secondSubJet->rapidity(); // x2
//       _phi1 = _firstSubJet->phi(); // y1
//       _phi2 = _secondSubJet->phi(); // y2
//       _dRapidity = _rapidity2 - _rapidity1; // ux = x2 - x1
//       _dPhi = _deltaPhi(_phi2, _phi1); // uy = y2 - y1
//       _jetSep2 = (_dRapidity * _dRapidity) + (_dPhi * _dPhi); // length = ux2 + uy2
//       _jetSep = sqrt(_jetSep2); // sqrt(length)
//       _et_J = _jet->et();
//       _sum = 0.;
//     }
//   }

//   ~JetDipolarity() {
//     for(unsigned int i = 0; i != _subjetsV.size(); i++) delete _subjetsV[i];
//     _subjetsV.clear();
//   }

//   void operator()(const AliPicoTrack i) {
//     term(i);
//   }

//   void operator()(const AliPicoTrack* i) {
//     term(*i);
//   }

//   double/* jStr_dipJ_*/dipolarity() const { // PUBLISH_JET_SHAPE
//     return Divide(_sum, _jetSep2);
//   }
  
// private:
//   // http://www.codeguru.com/forum/showthread.php?t=194400&page=2
//   // http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
//   // http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
//   void term(const AliPicoTrack i) {
//     if(_jet->size() <= 1 || _numSubjets < 2) return;
//     double et_i = i.et();
    
//     double rapidity_i = i.rapidity();
//     double dRapidity_i = _rapidity1 - rapidity_i; // vx = x1 - X
//     double phi_i = i.phi();
//     double dPhi_i = _deltaPhi(_phi1, phi_i); // vy = y1 - Y
    
//     double Ri;
//     double det = (dRapidity_i * _dRapidity) + (dPhi_i * _dPhi); // det = (vx * ux) + (vy * uy) = v.u = |v||u|cos(theta)
    
//     if(det < 0 || det > _jetSep) {
//       _dRapidity = _rapidity2 - rapidity_i; // ux = x2 - X
//       _dPhi = _deltaPhi(_phi2, phi_i); // uy = y2 - Y
//       Ri = sqrt(std::min(
// 			 (dRapidity_i * dRapidity_i) + (dPhi_i * dPhi_i), // vx2 + vy2
// 			 (_dRapidity * _dRapidity) + (_dPhi * _dPhi) // ux2 + uy2
// 			 ));
//     }
//     else {
//       Ri = Divide(fabs(dRapidity_i * _dPhi - _dRapidity * dPhi_i), _jetSep); // |vx * uy - ux * vy| / length
//     }
//     _sum += Divide(et_i, _et_J) * (Ri * Ri);
//   }

//   DeltaPhi _deltaPhi;
//   const AliEmcalJet* _jet;
//   const TClonesArray* _subjets;
//   int _numSubjets;
//   TClonesArrayHelper::TClonesArray_t _sortedSubJets;
//   std::vector<const Jet*> _subjetsV;
//   const AliEmcalJet* _firstSubJet;
//   const AliEmcalJet* _secondSubJet;
//   double _rapidity1;
//   double _rapidity2;
//   double _phi1;
//   double _phi2;
//   double _dRapidity;
//   double _dPhi;
//   double _jetSep2;
//   double _jetSep;
//   double _et_J;
//   double _sum;
// };

// // http://arxiv.org/pdf/1011.2268v3
// // https://indico.cern.ch/getFile.py/access?contribId=24&resId=0&materialId=slides&confId=138809
// // et --> pt 2011-11-22
// class Nsubjettiness {
// public:
//   Nsubjettiness(const Jet*& j, double beta = 1.):
//   _jet(j), _sortedSubJets(0), _beta(beta), _d0(0.)  {
//     if(_numConstituents <= 1) {
//       _d0 = NaN;
//     }
//     _sum.clear();
//     JetMisc misc(_jet);
//     _R = misc.R();
//     _subjets = misc.subjets();
//     _numSubjets = misc.numSubjets();
//     _numConstituents = misc.numConstituents();
//     if(_numSubjets) {
//       _sum.resize(_numSubjets);
//       for(unsigned int i = 0; i != _numSubjets; i++) {
// 	_sum[i] = 0.;
//       }
//     }
//     if(_subjets) {
//       TClonesArrayHelper::TClonesArray_t sortedSubJets(_subjets->begin(), _subjets->end());
//       TClonesArrayHelper::sort_jets(sortedSubJets, JetSorters::sortJetByEtDown());
//       _sortedSubJets = sortedSubJets;
//     }
//   }

//   void operator()(const AliPicoTrack i) {
//     term(i);
//   }

//   void operator()(const AliPicoTrack* i) {
//     term(*i);
//   }

//   const Jet* jet() const { return _jet; } 
//   double tauN(unsigned int N) const { // PUBLISH_JET_SHAPE
//     return N > 0 && _numSubjets >= N? Divide(_sum[N-1], _d0): 0.; // 2011-11-22
//     //return N > 0? Divide(_sum[N-1], _d0): 0.;
//   }

//   double/* jStr_nsubjnessJ_*/tau1() const { return tauN(1); } // PUBLISH_JET_SHAPE 
//   double/* jStr_nsubjnessJ_*/tau2() const { return tauN(2); } // PUBLISH_JET_SHAPE 
//   double/* jStr_nsubjnessJ_*/tau3() const { return tauN(3); } // PUBLISH_JET_SHAPE 
//   double/* jStr_nsubjnessJ_*/tau2tau1() const { return Divide(tau2(), tau1()); } // PUBLISH_JET_SHAPE 
//   double/* jStr_nsubjnessJ_*/tau3tau2() const {return Divide(tau3(), tau2()); } // PUBLISH_JET_SHAPE 
  
// private:
//   void term(const AliPicoTrack i) {
//     if(_numConstituents <= 1 || _numSubjets == 0) return;
  
//     // double et_i = i.et();
//     double pt_i = i.pt();
//     double rapidity_i = i.rapidity();
//     double phi_i = i.phi();
   
//     // _d0 += et_i * _R;
//     _d0 += pt_i * _R;
    
//     std::vector<double> deltaRs;
//     std::vector<double> minDeltaRs;
//     deltaRs.clear();
//     minDeltaRs.clear();
    
//     if(_subjets) {
//       TClonesArrayHelper::TClonesArray_t::iterator firstJet = _sortedSubJets.begin();
//       TClonesArrayHelper::TClonesArray_t::iterator lastJet  = _sortedSubJets.end();
//       for(; firstJet != lastJet; ++firstJet) {
// 	if((*firstJet) == 0) continue;
// 	const Jet* aJet = (*firstJet)->clone(true, true);
// 	if(aJet == 0) continue;
// 	double rapidity_J = aJet->rapidity();
// 	double phi_J = aJet->phi();
// 	double dRapidity_i = rapidity_J - rapidity_i;
// 	double dPhi_i = _deltaPhi(phi_J, phi_i);
// 	double deltaR = sqrt((dRapidity_i * dRapidity_i) + (dPhi_i * dPhi_i));
// 	deltaRs.push_back(_beta == 1.? deltaR: Power(deltaR, _beta));
// 	std::vector<double>::iterator it = std::min_element(deltaRs.begin(), deltaRs.end());
// 	double minDeltaR = *it;
// 	minDeltaRs.push_back(minDeltaR);
// 	delete aJet;
//       }
//     }
//     for(unsigned int j = 0; j != _numSubjets; j++) _sum[j] += minDeltaRs[j] * pt_i; // et_i
//   }

//   DeltaPhi _deltaPhi;
//   const Jet* _jet;
//   double _R;
//   const TClonesArray* _subjets;
//   unsigned int _numSubjets, _numConstituents;
//   TClonesArrayHelper::TClonesArray_t _sortedSubJets;
//   std::vector<double> _sum;
//   double _beta, _d0;
// };

// http://indico.cern.ch/conferenceDisplay.py?confId=114321
class Psi {
public:
  Psi(const AliEmcalJet*& j, double r, double deltaR = 0.1):
  _jet(j), _r(r), _deltaR(deltaR), _psi(0), _rho(0) {
    if(_jet->Nch() <= 1) {
      _psi = NaN;
      _rho = NaN;
    }
    TLorentzVector *hlv = new TLorentzVector;
    j->GetMom(*hlv);
    _theHLV = new CLHEP::HepLorentzVector(hlv->X(), hlv->Y(), hlv->Z(), hlv->T());
  }

  ~Psi() {
    delete _theHLV;
  }

  void operator()(const AliPicoTrack i) {
    term(i);
  }

  void operator()(const AliPicoTrack* i) {
    term(*i);
  }

  const AliEmcalJet* jet() const { return _jet; } 
  double/* jStr_psiJX_*/psi() const { return _psi; } // PUBLISH_JET_SHAPE
  double/* jStr_psiJX_*/rho() const { return _rho / _deltaR; } // PUBLISH_JET_SHAPE

private:
  void term(const AliPicoTrack i) {
    if(_jet->Nch() <= 1) return;
    CLHEP::Hep2Vector ri;
    double core, annulus;

    CLHEP::HepLorentzVector t = CLHEP::HepLorentzVector(i.Px(), i.Py(), i.Pz(), i.E());

    double dy = t.rapidity() - _theHLV->rapidity();
    double dphi = deltaPhi(t.phi(), _theHLV->phi());
    ri = CLHEP::Hep2Vector(dy, dphi);

    core = ri.mag() < _r? t.et(): 0;
    annulus = ri.mag() > (_r - (_deltaR / 2.)) && ri.mag() < (_r + (_deltaR / 2.))? t.et(): 0;

    double et_J = _theHLV->et();

    _psi += Divide(core, et_J);
    _rho += Divide(annulus, et_J);
  }

  CLHEP::HepLorentzVector* _theHLV;
  DeltaPhi deltaPhi;
  const AliEmcalJet* _jet;
  double _r, _deltaR, _psi, _rho;
};

// http://thesis.library.caltech.edu/6221/1/ThesisMain_D.pdf
class PsiRatios {
public:
  PsiRatios(const AliEmcalJet*& j): _psi1(j, 0.1), _psi2(j, 0.2), _psi3(j, 0.3), _psi7(j, 0.7) {}

  void operator()(const AliPicoTrack i) {
    term(i);
  }

  void operator()(const AliPicoTrack* i) {
    term(*i);
  }

  double/* jStr_psiRatiosJ_*/psi1() const { return _psi1.psi(); } // PUBLISH_JET_SHAPE
  double/* jStr_psiRatiosJ_*/psi2() const { return _psi2.psi(); } // PUBLISH_JET_SHAPE
  double/* jStr_psiRatiosJ_*/psi3() const { return _psi3.psi(); } // PUBLISH_JET_SHAPE
  double/* jStr_psiRatiosJ_*/psi7() const { return _psi7.psi(); } // PUBLISH_JET_SHAPE
  double/* jStr_psiRatiosJ_*/psi717() const { return Divide(psi7() - psi1(), psi7()); } // PUBLISH_JET_SHAPE
  double/* jStr_psiRatiosJ_*/psi127() const { return Divide(psi1() - psi2(), psi7()); } // PUBLISH_JET_SHAPE
  double/* jStr_psiRatiosJ_*/psi37() const { return Divide(psi3(), psi7()); } // PUBLISH_JET_SHAPE

private:
  void term(const AliPicoTrack i) {
    _psi1(i);
    _psi2(i);
    _psi3(i);
    _psi7(i);
  }

  Psi _psi1;
  Psi _psi2;
  Psi _psi3;
  Psi _psi7;
};

// YOU ARE HERE!

// class BCF {
// public:
//   BCF(const Jet*& j, int version = 0, unsigned int a = 1):
//   _jet(j),
//   _maxRapidity(1000.),
//   _version(version),
//   _a(a),
//   _sumDeltaY(0.), // orthogonal beam colour flow var
//   _sumDeltaPhi(0.), // original beam colour flow var
//   _sumInDeltaY(0.), _sumOutDeltaY(0.), _sumPosDeltaPhi(0.), _sumNegDeltaPhi(0.),
//   _rapidity1(_jet->rapidity()), // x1
//   _rapidity2(_rapidity1 >= 0.? _maxRapidity: -_maxRapidity), // x2
//   _phi1(_jet->phi()), // y1 (phi2 = phi1 => dPhi = uy = 0)
//   _dRapidity(_rapidity2 - _rapidity1), // ux = x2 - x1
//   _et_J(_jet->et()) { // CTOR body
//     if(_jet->size() <= 1) {
//       _sumDeltaY = _sumDeltaPhi = NaN;
//       _sumInDeltaY = _sumOutDeltaY = _sumPosDeltaPhi = _sumNegDeltaPhi = NaN;
//     }
//   }

//   void operator()(const AliPicoTrack i) {
//     term(i);
//   }

//   void operator()(const AliPicoTrack* i) {
//     term(*i);
//   }

//   const Jet* jet() const { return _jet; } 

//   int/* jStr_bcfJ_XX_*/bcfVersion() const { return _version; }
//   int/* jStr_bcfJ_XX_*/a() const { return _a; }
//   double/* jStr_bcfJ_XX_*/bcfT() const { return _sumDeltaY; } // PUBLISH_JET_SHAPE
//   double/* jStr_bcfJ_XX_*/bcf() const { return _sumDeltaPhi; } // PUBLISH_JET_SHAPE
//   double/* jStr_bcfJ_XX_*/bcfAsymY() const { return Asymmetry(_sumInDeltaY, _sumOutDeltaY); } // PUBLISH_JET_SHAPE
//   double/* jStr_bcfJ_XX_*/bcfAsymPhi() const { return Asymmetry(_sumPosDeltaPhi, _sumNegDeltaPhi); } // PUBLISH_JET_SHAPE
//   double/* jStr_bcfJ_XX_*/bcfAsymYPhi() const { return Asymmetry(_sumDeltaY, _sumDeltaPhi); } // PUBLISH_JET_SHAPE
//   double/* jStr_bcfJ_XX_*/bcfAsymYPhi2() const { return Asymmetry(_sumDeltaPhi, _sumOutDeltaY); } // PUBLISH_JET_SHAPE

// private:
//   void term(const AliPicoTrack i) {
//     if(_jet->size() <= 1) return;
//     double et_i = i.et();
    
//     double rapidity_i = i.rapidity();
//     double dRapidity_i = _rapidity1 - rapidity_i; // vx = x1 - X
//     double phi_i = i.phi();
//     double dPhi_i = _deltaPhi(_phi1, phi_i); // vy = y1 - Y
    
//     double Ri = 0.;
//     double Ti = 0.;
//     double InDeltaY = 0.;
//     double OutDeltaY = 0.;
//     double PosDeltaPhi = 0.;
//     double NegDeltaPhi = 0.;

//     double det = (dRapidity_i * _dRapidity); // det = (vx * ux) + (vy * uy) = v.u = |v||u| cos(theta)    

//     if(det < 0. || _rapidity1 == 0.) { // we're off the end of the line OR can't tell direction of line (+/-)
//       if(_rapidity1 != 0.) InDeltaY = fabs(dRapidity_i); // No contribution if we can't tell +ve from -ve

//       switch(_version) {
//       case 1:
// 	Ri = sqrt((dRapidity_i * dRapidity_i) + (dPhi_i * dPhi_i));
// 	Ti = 0.;
// 	break;
//       case 2:
// 	Ri = 0.;
// 	Ti = 0.;
// 	break;
//       default:
// 	Ri = fabs(dPhi_i); // |vx * uy - ux * vy| / length, originally Divide(fabs(dRapidity * dPhi_i), fabs(dRapidity));
// 	Ti = fabs(dRapidity_i);
// 	break;
//       }
//     }
//     else { // works even if dRapidity_i == 0 because det < 0 always false, regardless of sign of _dRapidity, because det == NaN
//       OutDeltaY = fabs(dRapidity_i);
//       Ri = fabs(dPhi_i); // |vx * uy - ux * vy| / length, originally Divide(fabs(dRapidity * dPhi_i), fabs(dRapidity));
//       Ti = fabs(dRapidity_i);
//     }

//     PosDeltaPhi = dPhi_i > 0.? fabs(dPhi_i): 0.;
//     NegDeltaPhi = dPhi_i < 0.? fabs(dPhi_i): 0.;  

//     _sumDeltaY += Divide(et_i, _et_J) * Power(Ti, _a);
//     _sumDeltaPhi += Divide(et_i, _et_J) * Power(Ri, _a);

//     _sumInDeltaY +=  Divide(et_i, _et_J) * Power(InDeltaY, _a);
//     _sumOutDeltaY +=  Divide(et_i, _et_J) * Power(OutDeltaY, _a);
//     _sumPosDeltaPhi +=  Divide(et_i, _et_J) * Power(PosDeltaPhi, _a);
//     _sumNegDeltaPhi +=  Divide(et_i, _et_J) * Power(NegDeltaPhi, _a);
//   }

//   const Jet* _jet;
//   double _maxRapidity;
//   int _version;
//   unsigned int _a;
//   DeltaPhi _deltaPhi;
//   double _sumDeltaY, _sumDeltaPhi;
//   double _sumInDeltaY, _sumOutDeltaY, _sumPosDeltaPhi, _sumNegDeltaPhi;
//   double _rapidity1;
//   double _rapidity2;
//   double _phi1;
//   double _dRapidity;
//   double _et_J;
// };

// http://prd.aps.org/pdf/PRD/v79/i7/e074017
// arXiv:1012.2077, arXiv:1006.2035v1 [hep-ph], arXiv:1012.5412v1 [hep-ph]
class PlanarFlow {
public:
  PlanarFlow(const AliEmcalJet*& j): _jet(j) {
    _zaxis = CLHEP::Hep3Vector(_jet->Px(), _jet->Py(), _jet->Pz());
    _zaxis *= 1. / _zaxis.mag();
    _zbeam = CLHEP::Hep3Vector(0., 0., 1.);
    _xaxis = _zaxis.cross(_zbeam);
    _xaxis *= 1. / _xaxis.mag();
    _yaxis = _xaxis.cross(_zaxis);
    _yaxis *= 1. / _yaxis.mag();
    _mJ = _jet->M();
    if(_jet->Nch() <= 1) {
      _I[0][0] = _I[0][1] = _I[1][0] = _I[1][1] = NaN;
      _s[0][0] = _s[0][1] = _s[1][0] = _s[1][1] = NaN;    
      _ST[0][0] = _ST[0][1] = _ST[1][0] = _ST[1][1] = NaN;
    }
    else {
      _I[0][0] = _I[0][1] = _I[1][0] = _I[1][1] = 0.;
      _s[0][0] = _s[0][1] = _s[1][0] = _s[1][1] = 0.;
      _ST[0][0] = _ST[0][1] = _ST[1][0] = _ST[1][1] = 0.;
    }
  }
  
  void operator()(const AliPicoTrack i) {
    term(i);
  }

  void operator()(const AliPicoTrack* i) {
    term(*i);
  }

  const AliEmcalJet* jet() const { return _jet; }
  double/* jStr_pfJ_*/pf() const { // PUBLISH_JET_SHAPE
    double D = (_I[0][0] * _I[1][1]) - (_I[0][1] * _I[1][0]);
    double T = _I[0][0] + _I[1][1];
    return Divide(4.* D, T * T);
  }

  // arXiv:0806.0023v2 [hep-ph]
  double/* jStr_pfJ_*/detST() { return ST()->determinant(); } // PUBLISH_JET_SHAPE
  double/* jStr_pfJ_*/lambdaST() { return ST()->eigenvalues().first; } // PUBLISH_JET_SHAPE
  
private:
  void term(const AliPicoTrack i) {
    if(_jet->Nch() <= 1) return;
    double  px, py, pmag, e;
    CLHEP::Hep3Vector p = CLHEP::Hep3Vector(i.Px(), i.Py(), i.Pz());
    px = p.dot(_xaxis);
    py = p.dot(_yaxis);
    e = i.E();

    //double w = e * _mJ;
    double w = e; // http://prd.aps.org/pdf/PRD/v79/i7/e074012

    _I[0][0] += Divide(px * px, w);
    _I[0][1] += Divide(px * py, w);
    _I[1][0] += Divide(py * px, w);
    _I[1][1] += Divide(py * py, w);

    pmag = Sqrt((px * px) + (py * py));

    _pmagsum += pmag;
    _s[0][0] += Divide(px * px, pmag);
    _s[0][1] += Divide(px * py, pmag);
    _s[1][0] += Divide(py * px, pmag);
    _s[1][1] += Divide(py * py, pmag);
  }
  
  const Hep2Matrix* I() const { return &_I; }
  const Hep2Matrix* ST() { DEBUG(_pmagsum);
    _ST[0][0] = _pmagsum? Divide(_s[0][0], _pmagsum): NaN;
    _ST[0][1] = _pmagsum? Divide(_s[0][1], _pmagsum): NaN;
    _ST[1][0] = _pmagsum? Divide(_s[1][0], _pmagsum): NaN;
    _ST[1][1] = _pmagsum? Divide(_s[1][1], _pmagsum): NaN;
    return &_ST;
  }
  
  const AliEmcalJet* _jet;
  Hep2Matrix _I, _s, _ST;
  CLHEP::Hep3Vector _zbeam, _xaxis, _yaxis, _zaxis;
  double _mJ, _px, _py, _pmagsum, _e;
};

// track-based energy-energy correlation angularities
// Light-quark and gluon jet discrimination in pp collisions at sqrt(s) = 7 TeV with the ATLAS detector
// http://arxiv.org/pdf/1405.6583v3.pdf
// http://xxx.tau.ac.il/pdf/1305.0007v2.pdf
class AngEEC { // PUBLISH_JET_SHAPE
public:
  double/* jStr_AngEECJ_*/operator()(const AliEmcalJet*& jet, TClonesArray*& fTracks, double beta) {
  Int_t Nch = jet->Nch(); DEBUG(Nch);     DEBUG("I am here");
  double _sumNumerator = 0.;
  double _sumPt = 0.;
      for(Int_t i = 0; i != Nch; i++) {
        AliVTrack *track_i = static_cast<AliVTrack*>(jet->TrackAt(i, fTracks));
        if(track_i) {
          double iPt = track_i->Pt();
          _sumPt += iPt; 
          for(Int_t j = 0; j != Nch; j++) {
            if(j <= i) continue;
            AliVTrack *track_j = static_cast<AliVTrack*>(jet->TrackAt(j, fTracks));
            if(track_j) {
              double dEta = fabs(track_i->Eta() - track_j->Eta());
              double dPhi = DeltaPhi()(track_i->Phi(), track_j->Phi());
              double dR = Sqrt((dEta * dEta) + (dPhi * dPhi));
              _sumNumerator += track_i->Pt() * track_j->Pt() * pow(dR, beta);
            }
          }
        }
      }
      return Divide(_sumNumerator, _sumPt * _sumPt);
    }
  };

  class TauEEC { // PUBLISH_JET_SHAPE
public:
  double/* jStr_AngEECJ_*/operator()(const AliEmcalJet*& jet, TClonesArray*& fTracks, double beta) {
  Int_t Nch = jet->Nch(); DEBUG(Nch);     DEBUG("I am here");
  double _sumNumerator = 0.;
  double _sumPt = 0.;
      for(Int_t i = 0; i != Nch; i++) {
        AliVTrack *track_i = static_cast<AliVTrack*>(jet->TrackAt(i, fTracks));
        if(track_i) {
          double iPt = track_i->Pt();
          _sumPt += iPt; 
          for(Int_t j = 0; j != Nch; j++) {
            if(j <= i) continue;
            AliVTrack *track_j = static_cast<AliVTrack*>(jet->TrackAt(j, fTracks));
            if(track_j) {
              double dEta = fabs(track_i->Eta() - track_j->Eta());
              double dPhi = DeltaPhi()(track_i->Phi(), track_j->Phi());
              double tau = Arctan2(fabs(dPhi), fabs(dEta));
              //double dR = Sqrt((dEta * dEta) + (dPhi * dPhi));
              _sumNumerator += track_i->Pt() * track_j->Pt() * pow(tau, beta);
            }
          }
        }
      }
      return Divide(_sumNumerator, _sumPt * _sumPt);
    }
  };

// ArXiv:0909.3855, arXiv:0806.0023v2 [hep-ph], also arXiv:0903.5081v4 [hep-ph] is similar (pt instead of et)
class ZminJ { // PUBLISH_JET_SHAPE
public:
  double/* jStr_zminJ_*/operator()(const AliEmcalJet*& j, TClonesArray*& fTracks) {
  Int_t Nch = j->Nch(); DEBUG(Nch);     DEBUG("I am here");
  double E, minE = std::numeric_limits<double>::infinity();
      for(Int_t i = 0; i != Nch; i++) {
        AliVTrack *track = static_cast<AliVTrack*>(j->TrackAt(i, fTracks));
        if(track) {
          E = track->E();
          if(E < minE) minE = E;
        }
      }
    return Divide(minE, j->E());
  }
};

// Idea loosely based on ArXiv:0909.3855
class ZmaxJ { // PUBLISH_JET_SHAPE
public:
  double/* jStr_zmaxJ_*/operator()(const AliEmcalJet*& j, TClonesArray*& fTracks) {
  Int_t Nch = j->Nch(); DEBUG(Nch);     DEBUG("I am here");
  double E, maxE = 0.;
      for(Int_t i = 0; i != Nch; i++) {
        AliVTrack *track = static_cast<AliVTrack*>(j->TrackAt(i, fTracks));
        if(track) {
          E = track->E();
          if(E > maxE) maxE = E;
        }
      }
    return Divide(maxE, j->E());
  }
};

  class R90 { // PUBLISH_JET_SHAPE
  public:
    R90(const AliEmcalJet*& j): 
    _eta(j->Eta()),
    _phi(j->Phi()),
    _jetP(j->P()),
    _jetPt(j->P()),
    _jetE(j->P()),
    _threshold(0.90 * _jetP) // This defines R90 to be 90% of jet's momentum
  {/*empty*/}

  // You should clean this up, but it's fine for some quick and dirty numbers...
  double/* jStr_R90_*/operator()(const AliEmcalJet*& j, TClonesArray*& fTracks) {
    Int_t Nch = j->Nch(); DEBUG(Nch);
    if(Nch == 1) return NaN; // Single track --> is R90 meaningful anymore?
  double P = 0, Pt = 0, E = 0, dR = 0;
  std::map<double, Int_t> sort_index;
  sort_index.clear();

  for(Int_t i = 0; i != Nch; i++) {
    AliPicoTrack* track = static_cast<AliPicoTrack*>(j->TrackAt(i, fTracks));
    if(track) {
  // How far is this track away from the jet axis?
      double dEta = fabs(track->Eta() - _eta);
      double dPhi = DeltaPhi()(track->Phi(), _phi);
      dR = Sqrt((dEta * dEta) + (dPhi * dPhi));
      // In order of increasing dR, sort the index of the tracks:
      sort_index[dR] = i;
      // There are better ways of doing this, but I have to workaround AliROOT's braindead contraints

    } // got track
  } // 1st loop over tracks

  std::map<double, Int_t>::iterator it; DEBUG("Calculate R90"); // This is dR
  DEBUG(Nch);
  for(it = sort_index.begin(); it != sort_index.end(); ++it) { // Iterate in order of increasing dR
    DEBUG(it->first);
    Int_t k = it->second; DEBUG(k); // This is the track index
    AliPicoTrack* track = static_cast<AliPicoTrack*>(j->TrackAt(k, fTracks)); // Get track
    if(track) { // Got a good track
      P += track->P(); DEBUG(P); // Change this if R90 is to mean something else
      DEBUG(_threshold);
      if(P > _threshold) break; // Stop once we've got 90% of jet momentum
      double dEta = fabs(track->Eta() - _eta);
      double dPhi = DeltaPhi()(track->Phi(), _phi);
      dR = Sqrt((dEta * dEta) + (dPhi * dPhi)); DEBUG(dR); // This is the dR that contains 90%
    } // got the correct track
 } // 2nd loop over tracks
 DEBUG("final dR:")
 DEBUG(dR); 
  return dR;
  }

  double _eta;
  double _phi;
  double _jetP;
  double _jetPt;
  double _jetE;
  double _threshold;
};

//arXiv: 1006.2035v1 [hep-ph] / 0807.0234v1 [hep-ph]
// Default version, uses default pow
class Tau {
public:
  Tau(double a, const AliEmcalJet*& j, double R = CLHEP::halfpi): _a(a), _jet(j), _R(R), _sum(0) {
    _J = CLHEP::Hep2Vector(j->Y(), phiCorr(j->Phi()));
    _mJinv = Divide(1., j->M());
    if(_jet->Nch() <= 1) _sum = NaN;
  }

  // Use standard C++ pow
  void operator()(const AliPicoTrack i) {
    term(i);
  }

  // Use standard C++ pow
  void operator()(const AliPicoTrack* i) {
    term(*i);
  }

  const AliEmcalJet* jet() const { return _jet; }
  double/* jStr_tauJXX_*/a() const { return _a; } 
  double/* jStr_tauJXX_*/tau() const { return _mJinv * _sum; } // PUBLISH_JET_SHAPE
  std::pair<double, double> tauPair() const { return std::make_pair(_a, tau()); } // PUBLISH_JET_SHAPE

private:
  // Use standard C++ pow
  void term(const AliPicoTrack i) {
    if(_jet->Nch() <= 1) return; //DEBUG(_jet->size());
    CLHEP::Hep2Vector _ci = CLHEP::Hep2Vector(i.Y(), phiCorr(i.Phi())); //DEBUG(0);
    double theta = Arccos(Divide(_ci.dot(_J), _ci.mag() * _J.mag()));  //DEBUG(theta);
    double ei = i.E(); //DEBUG(ei);
    double sin = 0.; //DEBUG(0);
    double cos = 0.; //DEBUG(0);
    sincos(((CLHEP::halfpi * theta) / _R), &sin, &cos); //DEBUG(0);
    //DEBUG(sin); //DEBUG(_a);
    //DEBUG(1. - cos); //DEBUG( 1. - _a);
    //DEBUG(Power(sin, _a));
    //DEBUG(Power(1. - cos, 1. - _a));
    _sum += ei * Multiply(Power(sin, _a), Power(1. - cos, 1. - _a)); //DEBUG(_sum);
  }  

  const double _a;
  PhiCorr phiCorr;
  const AliEmcalJet* _jet;
  const double _R;
  double _mJinv, _sum;
  CLHEP::Hep2Vector _J;
};

// Fast version for integer powers only
// template<int a = -2>
// class iTau {
// public:
//   iTau(const Jet*& j, double R = CLHEP::halfpi): _jet(j), _R(R), _sum(0) {
//     _J = Hep2Vector(j->rapidity(), phiCorr(j->phi()));
//     _mJinv = Divide(1., j->m());
//     if(_jet->size() <= 1) _sum = NaN;
//   }

//   void operator()(const AliPicoTrack i) {
//     term(i);
//   }

//   void operator()(const AliPicoTrack* i) {
//     term(*i);
//   }

//   const Jet* jet() const { return _jet; } 
//   double tau() const { return _mJinv * _sum; } // PUBLISH_JET_SHAPE

// private:
//   void term(const AliPicoTrack i) {
//     if(_jet->size() <= 1) return;
//     CLHEP::Hep2Vector _ci = Hep2Vector(i.rapidity(), phiCorr(i.phi()));
//     double theta = Arccos(Divide(_ci.dot(_J), _ci.mag() * _J.mag()));
//     double ei = i.e();
//     double sin = 0.;
//     double cos = 0.;
//     sincos(((CLHEP::halfpi * theta) / _R), &sin, &cos);
//     _sum += ei * intpow<a>(sin) * intpow<1 - a>(1. - cos);
//   }

//   PhiCorr phiCorr;
//   const Jet* _jet;
//   const double _R;
//   double _mJinv, _sum;
//   CLHEP::Hep2Vector _J;
// };

// Fast version for non-integer powers:
// template<double *a>
// class dTau {
// public:
//   dTau(const Jet*& j, double R = CLHEP::halfpi): _jet(j), _R(R), _sum(0) {
//     _J = Hep2Vector(j->rapidity(), phiCorr(j->phi()));
//     _mJinv = Divide(1., j->m());
//     if(_jet->size() <= 1) _sum = NaN;
//   }

//   void operator()(const AliPicoTrack i) {
//     term(i);
//   }

//   void operator()(const AliPicoTrack* i) {
//     term(*i);
//   }

//   const Jet* jet() const { return _jet; } 
//   double tau() const { return _mJinv * _sum; } // PUBLISH_JET_SHAPE

// private:
//   void term(const AliPicoTrack i) {
//     if(_jet->size() <= 1) return;
//     CLHEP::Hep2Vector _ci = Hep2Vector(i.rapidity(), phiCorr(i.phi()));
//     double theta = Arccos(Divide(_ci.dot(_J), _ci.mag() * _J.mag()));
//     double ei = i.e();
//     double sin = 0.;
//     double cos = 0.;
//     sincos(((CLHEP::halfpi * theta) / _R), &sin, &cos);
//     _sum += ei * FASTPOW(sin, *a) * FASTPOW(1. - cos, 1 - *a);
//   }

//   PhiCorr phiCorr;
//   const Jet* _jet;
//   const double _R;
//   double _mJinv, _sum;
//   CLHEP::Hep2Vector _J;
// };

#endif


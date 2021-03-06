#ifndef SPHERICITYSPHEROCITY_H
#define SPHERICITYSPHEROCITY_H

#include "AliEmcalJet.h"

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliVParticle.h"

//#include "JetEvent/Jet.h"
//#include "JetEvent/JetConstituentIterator.h" // For Fox-Wolfram in CM frame exclusively
#include "CLHEP/Vector/ThreeVector.h" // Common to most variables
//#include "CLHEP/Vector/TwoVector.h" // For transverse spherocity
#include "CLHEP/Vector/LorentzVector.h" // For event shapes only
#include "CLHEP/Matrix/Matrix.h" // For SphericitySpherocity

#include "/home/andy/testarea/Hep2Matrix.h" // Linearised momentum tensor
#include "/home/andy/testarea/EigensolverSym3x3.h" // For SphericitySpherocity
#include "/home/andy/testarea/Numeric.h" // Common numeric stuff
#include <algorithm> // for std::min/max_element
#include <cmath> // For isnan(x)

// debugging macros so we can pin down message provenance at a glance
#include <iostream>
#define DEBUG(x)							\
  std::cout << "DEBUG:\t" << "\t" << #x << "\t" << x << "\t"		\
  << __FUNCTION__ << "\t" << __FILE__ << "\t" << __LINE__ << std::endl;
#define DEBUG(x)

#define SPH_NOBOOST true

using namespace Numeric;

// arXiv:1010.3698v1 [hep-ph]
class SphericitySpherocity {
public:
  SphericitySpherocity(const TClonesArray* fJets, TClonesArray*& fTracks, char* instName = "", const AliEmcalJet* parentJet = 0, bool useConstituents = false): 
  _theJetCollection(fJets),
  _qtVect(0), _3momenta(0), _qtDotN(0), _qtCrossN(0), _thrusts(0), _thrustAxes(0),
  _aco(0.), _pmag(0), _pmag2(0), _pmagSum(0), _pmag2Sum(0), _pTsum(0), _pT2sum(0), _sumI(0), _sumJ(0), _sumB(NaN), _sumQt(0.), _maxpT(0), _nConstituents(0),
  _sphericityLambda1(NaN), _sphericityLambda2(NaN), _sphericityLambda3(NaN),
  _spherocityLambda1(NaN), _spherocityLambda2(NaN), _spherocityLambda3(NaN),
  _minQtCrossN(NaN), _maxQtCrossN(NaN), _minQtDotN(NaN), _maxQtDotN(NaN),
  _fTracks(fTracks),
  _instName(instName) {
    _sphericityTensor = CLHEP::HepMatrix(3, 3, 0.);
    _spherocityTensor = CLHEP::HepMatrix(3, 3, 0.);
    _Mij = CLHEP::HepMatrix(4,4,0.);
    _Nij = CLHEP::HepMatrix(4,4,0.);
    _mu = CLHEP::HepMatrix(4,4,0.);
    _eta = CLHEP::HepMatrix(4,4,0.);
    _nu = CLHEP::HepMatrix(4,4,0.);
    _zeta = CLHEP::HepMatrix(4,4,0.);
    _ThetaHu = _eHu = _I1 = _I2 = _I3 = _I4 = _I5 = _I6 = _I7 = _I8 = NaN;
    _ThetaHu2 = _eHu2 = _J1 = _J2 = _J3 = _J4 = _J5 = _J6 = _J7 = _J8 = NaN;

    _Mlin = Hep2Matrix(0.);
    
    _qtVect = new std::vector<CLHEP::Hep3Vector*>;
    _qtVect->clear(); DEBUG(0);

    _3momenta = new std::vector<CLHEP::Hep3Vector*>;
    _3momenta->clear(); DEBUG(0);

    _thrustAxes = new std::vector<CLHEP::Hep3Vector*>;
    _thrustAxes->clear(); DEBUG(0);

    _thrusts = new std::vector<double>;
    _thrusts->clear(); DEBUG(0);

    _qtDotN = new std::vector<double>;
    _qtDotN->clear(); DEBUG(0);
    _qtCrossN = new std::vector<double>;
    _qtCrossN->clear(); DEBUG(0);

    if(_theJetCollection) { DEBUG(_theJetCollection);
      compute(parentJet, useConstituents); DEBUG(parentJet); DEBUG(useConstituents);
      calcHu();
      calcHu2();
      
      double lambda[3] = {NaN, NaN, NaN};
      if(normalisedEigenvalues(sphericityTensor(), lambda)) {
	_sphericityLambda1 = lambda[0];
	_sphericityLambda2 = lambda[1];
	_sphericityLambda3 = lambda[2];
      }
      
      if(normalisedEigenvalues(spherocityTensor(), lambda)) {
	_spherocityLambda1 = lambda[0];
	_spherocityLambda2 = lambda[1];
	_spherocityLambda3 = lambda[2];
      }
    } // IF we have a jet collection
    else { // no jet collection 
      compute(0, useConstituents); DEBUG(useConstituents);
      calcHu();
      calcHu2();

      double lambda[3] = {NaN, NaN, NaN};
      if(normalisedEigenvalues(sphericityTensor(1), lambda)) {
  _sphericityLambda1 = lambda[0];
  _sphericityLambda2 = lambda[1];
  _sphericityLambda3 = lambda[2];
      }
      
      if(normalisedEigenvalues(spherocityTensor(1), lambda)) {
  _spherocityLambda1 = lambda[0];
  _spherocityLambda2 = lambda[1];
  _spherocityLambda3 = lambda[2];
      }
    } // IF no jet collection
  }
  
  ~SphericitySpherocity() { DEBUG(_qtVect); DEBUG(_instName);
    if(_qtVect) { DEBUG(_qtVect);
      for(unsigned int i = 0; i != _qtVect->size(); i++) { DEBUG(_qtVect->size());
        delete _qtVect->at(i); DEBUG(0);
        _qtVect->at(i) = 0; DEBUG(i);
      }
      _qtVect->clear(); DEBUG(0);
    }
    delete _qtVect; DEBUG(0);
    _qtVect = 0; DEBUG(0);
    
    if(_3momenta) { DEBUG(_3momenta);
      for(unsigned int i = 0; i != _3momenta->size(); i++) { DEBUG(_3momenta->size());
        delete _3momenta->at(i); DEBUG(0);
        _3momenta->at(i) = 0; DEBUG(i);
      }
      _3momenta->clear(); DEBUG(0);
    }
    delete _3momenta; DEBUG(0);
    _3momenta = 0; DEBUG(0);
    
    if(_thrustAxes) { DEBUG(_thrustAxes);
      for(unsigned int i = 0; i != _thrustAxes->size(); i++) { DEBUG(_thrustAxes->size());
        delete _thrustAxes->at(i); DEBUG(0);
        _thrustAxes->at(i) = 0; DEBUG(i);
      }
      _thrustAxes->clear(); DEBUG(0);
    }
    delete _thrustAxes; DEBUG(0);
    _thrustAxes = 0; DEBUG(0);
    
    if(!_thrusts->empty()) _thrusts->clear(); DEBUG(_thrusts->size());
    if(!_qtDotN->empty()) _qtDotN->clear(); DEBUG(_qtDotN->size());
    if(_qtDotN) delete _qtDotN; _qtDotN = 0; DEBUG(_qtDotN);
    if(!_qtCrossN->empty())_qtCrossN->clear(); DEBUG(_qtCrossN->size());
    if(_qtCrossN) delete _qtCrossN; _qtCrossN = 0; DEBUG(_qtCrossN);
  }


  // http://arxiv.org/pdf/0806.0023
  double/* jStr_sphXX_*/detSphericity() const { return sphericityTensor()? sphericityTensor()->determinant(): NaN; }
  double/* jStr_sphXX_*/detSpherocity() const { return spherocityTensor()? spherocityTensor()->determinant(): NaN; }

  double/* jStr_sphXX_*/sphericityLambda1() const { return _sphericityLambda1; }
  double/* jStr_sphXX_*/sphericityLambda2() const { return _sphericityLambda2; }
  double/* jStr_sphXX_*/sphericityLambda3() const { return _sphericityLambda3; }
  double/* jStr_sphXX_*/spherocityLambda1() const { return _spherocityLambda1; } 
  double/* jStr_sphXX_*/spherocityLambda2() const { return _spherocityLambda2; }
  double/* jStr_sphXX_*/spherocityLambda3() const { return _spherocityLambda3; }

  double/* jStr_HuXX_*/HuM00() const { return _Mij[0][0]; }
  double/* jStr_HuXX_*/HuI1() const { return _I1; }
  double/* jStr_HuXX_*/HuI2() const { return _I2; }
  double/* jStr_HuXX_*/HuI3() const { return _I3; }
  double/* jStr_HuXX_*/HuI4() const { return _I4; }
  double/* jStr_HuXX_*/HuI5() const { return _I5; }
  double/* jStr_HuXX_*/HuI6() const { return _I6; }
  double/* jStr_HuXX_*/HuI7() const { return _I7; }
  double/* jStr_HuXX_*/HuI8() const { return _I8; }
  double/* jStr_HuXX_*/HuH1() const { return sgn(_I1) * log(fabs(_I1)); }
  double/* jStr_HuXX_*/HuH2() const { return sgn(_I2) * log(fabs(_I2)); }
  double/* jStr_HuXX_*/HuH3() const { return sgn(_I3) * log(fabs(_I3)); }
  double/* jStr_HuXX_*/HuH4() const { return sgn(_I4) * log(fabs(_I4)); }
  double/* jStr_HuXX_*/HuH5() const { return sgn(_I5) * log(fabs(_I5)); }
  double/* jStr_HuXX_*/HuH6() const { return sgn(_I6) * log(fabs(_I6)); }
  double/* jStr_HuXX_*/HuH7() const { return sgn(_I7) * log(fabs(_I7)); }
  double/* jStr_HuXX_*/HuH8() const { return sgn(_I8) * log(fabs(_I8)); }
  double/* jStr_HuXX_*/HuTheta() const { return _ThetaHu; }
  double/* jStr_HuXX_*/HuEcc() const { return _eHu; }

  double/* jStr_HuXX_*/Humu00() const { return _mu[0][0]; }
  double/* jStr_HuXX_*/Humu11() const { return _mu[1][1]; }
  double/* jStr_HuXX_*/Humu20() const { return _mu[2][0]; }
  double/* jStr_HuXX_*/Humu02() const { return _mu[0][2]; }
  double/* jStr_HuXX_*/Humu21() const { return _mu[2][1]; }
  double/* jStr_HuXX_*/Humu12() const { return _mu[1][2]; }
  double/* jStr_HuXX_*/Humu30() const { return _mu[3][0]; }
  double/* jStr_HuXX_*/Humu03() const { return _mu[0][3]; }

  double/* jStr_HuXX_*/Hueta00() const { return _eta[0][0]; }
  double/* jStr_HuXX_*/Hueta11() const { return _eta[1][1]; }
  double/* jStr_HuXX_*/Hueta20() const { return _eta[2][0]; }
  double/* jStr_HuXX_*/Hueta02() const { return _eta[0][2]; }
  double/* jStr_HuXX_*/Hueta21() const { return _eta[2][1]; }
  double/* jStr_HuXX_*/Hueta12() const { return _eta[1][2]; }
  double/* jStr_HuXX_*/Hueta30() const { return _eta[3][0]; }
  double/* jStr_HuXX_*/Hueta03() const { return _eta[0][3]; }

  double/* jStr_HuXX_*/qHuM00() const { return _Nij[0][0]; }
  double/* jStr_HuXX_*/qHuI1() const { return _J1; }
  double/* jStr_HuXX_*/qHuI2() const { return _J2; }
  double/* jStr_HuXX_*/qHuI3() const { return _J3; }
  double/* jStr_HuXX_*/qHuI4() const { return _J4; }
  double/* jStr_HuXX_*/qHuI5() const { return _J5; }
  double/* jStr_HuXX_*/qHuI6() const { return _J6; }
  double/* jStr_HuXX_*/qHuI7() const { return _J7; }
  double/* jStr_HuXX_*/qHuI8() const { return _J8; }
  double/* jStr_HuXX_*/qHuH1() const { return sgn(_J1) * log(fabs(_J1)); }
  double/* jStr_HuXX_*/qHuH2() const { return sgn(_J2) * log(fabs(_J2)); }
  double/* jStr_HuXX_*/qHuH3() const { return sgn(_J3) * log(fabs(_J3)); }
  double/* jStr_HuXX_*/qHuH4() const { return sgn(_J4) * log(fabs(_J4)); }
  double/* jStr_HuXX_*/qHuH5() const { return sgn(_J5) * log(fabs(_J5)); }
  double/* jStr_HuXX_*/qHuH6() const { return sgn(_J6) * log(fabs(_J6)); }
  double/* jStr_HuXX_*/qHuH7() const { return sgn(_J7) * log(fabs(_J7)); }
  double/* jStr_HuXX_*/qHuH8() const { return sgn(_J8) * log(fabs(_J8)); }
  double/* jStr_HuXX_*/qHuTheta() const { return _ThetaHu2; }
  double/* jStr_HuXX_*/qHuEcc() const { return _eHu2; }

  double/* jStr_HuXX_*/Huqmu00() const { return _nu[0][0]; }
  double/* jStr_HuXX_*/Huqmu11() const { return _nu[1][1]; }
  double/* jStr_HuXX_*/Huqmu20() const { return _nu[2][0]; }
  double/* jStr_HuXX_*/Huqmu02() const { return _nu[0][2]; }
  double/* jStr_HuXX_*/Huqmu21() const { return _nu[2][1]; }
  double/* jStr_HuXX_*/Huqmu12() const { return _nu[1][2]; }
  double/* jStr_HuXX_*/Huqmu30() const { return _nu[3][0]; }
  double/* jStr_HuXX_*/Huqmu03() const { return _nu[0][3]; }

  double/* jStr_HuXX_*/Huqeta00() const { return _zeta[0][0]; }
  double/* jStr_HuXX_*/Huqeta11() const { return _zeta[1][1]; }
  double/* jStr_HuXX_*/Huqeta20() const { return _zeta[2][0]; }
  double/* jStr_HuXX_*/Huqeta02() const { return _zeta[0][2]; }
  double/* jStr_HuXX_*/Huqeta21() const { return _zeta[2][1]; }
  double/* jStr_HuXX_*/Huqeta12() const { return _zeta[1][2]; }
  double/* jStr_HuXX_*/Huqeta30() const { return _zeta[3][0]; }
  double/* jStr_HuXX_*/Huqeta03() const { return _zeta[0][3]; }

  // CMS-AN-2012/004
  // http://www.phys.uniroma1.it/DipWeb/dottorato/DOTT_FISICA/MENU/03DOTTORANDI/Seminari24/Semin_24_ott2011/Pandolfi2011.pdf
  // http://castello.web.cern.ch/castello/cms-note-higgs2l2j.pdf
  double/* jStr_sphXX_*/pft() const {
    return sqrt(Divide(_pT2sum, _pTsum * _pTsum));;
  }

  // https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=166288
  // pt-sum that suppressed forward tracks
  // http://arxiv.org/pdf/1005.4060v2
  double/* jStr_sphXX_*/beamThrust() const {
    return _sumB;
  }

  // http://elpub.bib.uni-wuppertal.de/servlets/DerivateServlet/Derivate-1480/dc1006.pdf
  double/* jStr_sphXX_*/circularity() const {
    return 2.* Divide(std::min(_sphericityLambda1, _sphericityLambda2), _sphericityLambda1 + _sphericityLambda2);
  }
  
  double/* jStr_sphXX_*/sphericity() const {
    return (3./2.) * (_sphericityLambda2 + _sphericityLambda3);
  }
  
  double/* jStr_sphXX_*/spherocity() const {
    return (3./2.) * (_spherocityLambda2 + _spherocityLambda3);
  }
  
  double/* jStr_sphXX_*/aplanarity() const {
    return (3./2.) * _sphericityLambda3;
  }
  
  double/* jStr_sphXX_*/aplanority() const {
    return (3./2.) * _spherocityLambda3;
  }
  
  double/* jStr_sphXX_*/Y() const {
    return (sqrt(3.)/2.) * (_sphericityLambda2 - _sphericityLambda3);
  }
  
  double/* jStr_sphXX_*/planarity() const {
    return (_sphericityLambda2 - _sphericityLambda3);
  }

  double/* jStr_sphXX_*/planority() const {
    return (_spherocityLambda2 - _spherocityLambda3);
  }

  double/* jStr_sphXX_*/Dshape() const {
    return 27.* _spherocityLambda1 * _spherocityLambda2 * _spherocityLambda3;
  }

  // http://www-conf.kek.jp/dis06/transparencies/WG4/hfs-tasevsky.ppt
  // http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node213.html
  double/* jStr_sphXX_*/Cshape() const {
    return 3.* ((_spherocityLambda1 * _spherocityLambda2) + (_spherocityLambda2 * _spherocityLambda3) + (_spherocityLambda3 * _spherocityLambda1));
  }

  // Check; C = 1 - 2nd Fox-Wolfram moment
  double/* jStr_sphXX_*/H2() const { 
    return 1. - Cshape();
  }
  // http://arxiv.org/pdf/1001.4082
  double/* jStr_sphXX_*/Fshape() const { 
    double e1 = _Mlin.eigenvalues().first;
    double e2 = _Mlin.eigenvalues().second;
    return Divide(e2, e1);
  }
  // http://arxiv.org/pdf/1001.4082
  double/* jStr_sphXX_*/G() const { 
    double f = Fshape();
    return Divide(4. * f, (1. + f) * (1. + f));
  }
  // http://arxiv.org/pdf/1001.4082 2D transverse sphericity
  double/* jStr_sphXX_*/ST2D() const { 
    double e1 = _Mlin.eigenvalues().first;
    double e2 = _Mlin.eigenvalues().second;
    return Divide(2. * e2, e1 + e2);
  }

  // http://arxiv.org/pdf/1205.3963v1.pdf
  // https://indico.cern.ch/event/100963/session/6/contribution/76/material/paper/1.pdf
  // http://arxiv-web3.library.cornell.edu/pdf/1110.2278v1.pdf
  // http://arxiv.org/pdf/1110.2278v1.pdf
  // Ortiz' tensor has an extra factor of 1/summed pT outside.
  double/* jStr_sphXX_*/sT2D() const {
    Hep2Matrix ortiz = Hep2Matrix(0.);
    ortiz[0][0] = _Mlin[0][0];
    ortiz[0][1] = _Mlin[0][1];
    ortiz[1][0] = _Mlin[1][0];
    ortiz[1][1] = _Mlin[1][1];
    ortiz /= _pTsum;
    double e1 = ortiz.eigenvalues().first;
    double e2 = ortiz.eigenvalues().second;
    return Divide(2. * e2, e1 + e2);
  }

  double/* jStr_sphXX_*/nConstituents() const {
    return _nConstituents;
  }

  double/* jStr_sphXX_*/meanpT() const {
    return Divide(_pTsum, _nConstituents);
  }

  double/* jStr_sphXX_*/maxpT() const {
    return _maxpT;
  }

  // http://arxiv.org/pdf/1001.4082 2D transverse spherOcity
  double/* jStr_sphXX_*/SoT2D() const {
  double s = NaN;
  if(_nConstituents < 2) s = NaN;
  if(_nConstituents == 2) s = 1.;
  double m = Divide(_minQtCrossN, _sumQt); DEBUG(_minQtCrossN); DEBUG(_sumQt);
  if(_nConstituents > 2) s = 0.25 * (CLHEP::pi * CLHEP::pi) * m * m; DEBUG(s);
  //if(s < 0. || s > 1.) return NaN;
  return s;
  }

  double/* jStr_sphXX_*/detMlin() const { 
    return _Mlin.determinant();
  }


    ///@{ Thrust scalar accessors
    /// The thrust scalar, \f$ T \f$, (maximum thrust).
    double thrust() const { return _thrusts->empty()? NaN: _thrusts->at(0); }
    /// The thrust major scalar, \f$ M \f$, (thrust along thrust major axis).
    double thrustMajor() const { return _thrusts->empty()? NaN: _thrusts->at(1); }
    /// The thrust minor scalar, \f$ m \f$, (thrust along thrust minor axis).
    double thrustMinor() const { return _thrusts->empty()? NaN: _thrusts->at(2); }
    /// The oblateness, \f$ O = M - m \f$ .
    double oblateness() const { return _thrusts->empty()? NaN: _thrusts->at(1) - _thrusts->at(2); }
    ///@}

    double acoplanarity() const { return _aco; }

private:
    CLHEP::HepLorentzVector& boost(CLHEP::HepLorentzVector i, const AliEmcalJet* parentJet = 0) { DEBUG(0);
      if(SPH_NOBOOST) return i;
    if(!parentJet) return i; DEBUG(parentJet);
    TLorentzVector *hlv= new TLorentzVector; DEBUG(hlv);
    parentJet->GetMom(*hlv);
    
    //TLorentzVector hlv;
    //parentJet->GetMom(hlv);
    CLHEP::HepLorentzVector theHLV = CLHEP::HepLorentzVector(hlv->X(), hlv->Y(), hlv->Z(), hlv->T());
    if(parentJet != 0 && (theHLV.restMass2() <= m2min || theHLV.findBoostToCM().mag2() >= b2max)) { // !ZMxpvTachyonic
      i.setX(NaN);
      i.setY(NaN);
      i.setZ(NaN);
      i.setT(NaN);
      parentJet = 0; // signal that we are to not boost, pass back NaN HLV instead
    }
    return !parentJet? i: i.boost(theHLV.findBoostToCM());
  }

  void compute(const AliEmcalJet* parentJet = 0, bool useConstituents = false) { DEBUG(_theJetCollection); DEBUG(parentJet); DEBUG(useConstituents);
    if((!_theJetCollection || _theJetCollection->GetEntriesFast() == 0) && useConstituents == false) {
      return;
    }
    double px, py, pz;
    
    std::vector<CLHEP::Hep3Vector> threeMomenta;
    threeMomenta.clear();

        if(useConstituents && !parentJet) { // loop over tracks in event
      _sumB = 0.; DEBUG(_fTracks); DEBUG(_theJetCollection); DEBUG(parentJet); DEBUG(useConstituents);

      Int_t ntracks = _fTracks->GetEntriesFast(); DEBUG(ntracks);
      // Iterate over its constituents...
      for(Int_t i = 0; i != ntracks; i++) {
        AliVTrack *track = static_cast<AliVTrack*>(_fTracks->At(i));
        //AliVTrack *track = static_cast<AliVTrack*>(jet->TrackAt(i, _fTracks));
        if(!track) continue;
        //TLorentzVector *hlv= new TLorentzVector;
        //AliVParticle *part = static_cast<AliVParticle*>(jet->TrackAt(i, fTracks));
        CLHEP::HepLorentzVector theHLV = CLHEP::HepLorentzVector(track->Px(), track->Py(), track->Pz(), track->E());
        px = boost(theHLV, parentJet).px();
        py = boost(theHLV, parentJet).py();
        pz = boost(theHLV, parentJet).pz();
        fillTensors(px, py, pz);
        fillMij(px, py, pz);
        _Mij /= _sumI;
        _Nij /= _sumJ;

        double pt = boost(theHLV, parentJet).perp();
        if(pt > _maxpT) _maxpT = pt;

        threeMomenta.push_back(boost(theHLV, parentJet).v());

        //CLHEP::Hep2Vector qt = CLHEP::Hep2Vector(theHLV.v()); // temp
        CLHEP::Hep3Vector* qt = new CLHEP::Hep3Vector( boost(theHLV, parentJet).v() ); DEBUG(parentJet); DEBUG(theHLV); DEBUG(qt);DEBUG(*qt);
        
        qt->setZ(0); DEBUG(*qt);
        _sumQt += qt->mag(); DEBUG(qt->mag()); DEBUG(_sumQt);
        _qtVect->push_back(qt); DEBUG(_qtVect->size());

      } DEBUG("tensors filled"); DEBUG(_sumQt);

      _nConstituents = fillVectors(); DEBUG(_nConstituents);
      calcThrust(threeMomenta); DEBUG(threeMomenta.size());
      _aco = calcAcoplanarity(threeMomenta); DEBUG(threeMomenta.size());
      // Int_t Nch = parentJet->Nch();
      // for(Int_t i 0; i != Nch; i++) {
      //   AliVTrack *track = static_cast<AliVTrack*>(parentJet->TrackAt(i, _fTracks));
      //   if(!track) continue;
      //   CLHEP::HepLorentzVector theHLV = CLHEP::HepLorentzVector(track->Px(), track->Py(), track->Pz(), track->E());
      //   px = boost(theHLV, parentJet).px();
      //   py = boost(theHLV, parentJet).py();
      //   pz = boost(theHLV, parentJet).pz();
      //   fillTensors(px, py, pz);
      } 
    else if(useConstituents && parentJet && parentJet->Nch()) {     DEBUG("I am here");// loop over constituents in jet boosted into own rest frame
      _sumB = 0.; DEBUG(_fTracks); DEBUG(_theJetCollection); DEBUG(parentJet); DEBUG(useConstituents);


      Int_t Nch = parentJet->Nch(); DEBUG(Nch);     DEBUG("I am here");
      for(Int_t i = 0; i != Nch; i++) {
        AliVTrack *track = static_cast<AliVTrack*>(parentJet->TrackAt(i, _fTracks));
        if(!track) continue;
        CLHEP::HepLorentzVector theHLV = CLHEP::HepLorentzVector(track->Px(), track->Py(), track->Pz(), track->E());
        px = boost(theHLV, parentJet).px();
        py = boost(theHLV, parentJet).py();
        pz = boost(theHLV, parentJet).pz();
        fillTensors(px, py, pz); DEBUG("tensors filled");
        double jetEta = parentJet->Eta();
        double jetPhi = parentJet->Phi();
        //double q = (track->Charge() > 0.? 2.: 1.); 
        double q = track->Charge(); 
        fillMij(px, py, pz, jetEta, jetPhi, q);
        _Mij /= _sumI;
        _Nij /= _sumJ;

        double pt = boost(theHLV, parentJet).perp();
        if(pt > _maxpT) _maxpT = pt;

        threeMomenta.push_back(boost(theHLV, parentJet).v());

        //CLHEP::Hep2Vector qt = CLHEP::Hep2Vector(theHLV.v()); // temp
        CLHEP::Hep3Vector* qt = new CLHEP::Hep3Vector( boost(theHLV, parentJet).v() ); DEBUG(parentJet); DEBUG(theHLV); DEBUG(qt);DEBUG(*qt);
        
        qt->setZ(0); DEBUG(*qt);
        _sumQt += qt->mag(); DEBUG(qt->mag()); DEBUG(_sumQt);
        _qtVect->push_back(qt); DEBUG(_qtVect->size());
       
      }     DEBUG("I am here");

      _nConstituents = fillVectors(); DEBUG(_nConstituents);
      calcThrust(threeMomenta); DEBUG(threeMomenta.size());
      _aco = calcAcoplanarity(threeMomenta); DEBUG(threeMomenta.size());
 //      JetConstituentIterator firstConstituent = JetConstituentIterator::first(parentJet);
 //      JetConstituentIterator lastConstituent = JetConstituentIterator::last(parentJet);
 //      _sumB = 0.;
 //      for(; firstConstituent != lastConstituent; ++firstConstituent) {
	// px = boost(firstConstituent.hlv(), parentJet).px();
	// py = boost(firstConstituent.hlv(), parentJet).py();
	// pz = boost(firstConstituent.hlv(), parentJet).pz();
	// fillTensors(px, py, pz);
 //      }
    }
    else { // loop over jets, or subjets if looking at jet boosted into own rest frame
    DEBUG("I am here");



      Int_t njets = _theJetCollection->GetEntriesFast(); DEBUG(njets);
      const TClonesArray* theJets = _theJetCollection; DEBUG(theJets);
      _sumB = 0.; DEBUG(_fTracks); DEBUG(_theJetCollection); DEBUG(parentJet); DEBUG(useConstituents);
      for(int i = 0; i != njets; i++) {

          const AliEmcalJet* aJet = dynamic_cast<AliEmcalJet*>(theJets->At(i)); //DEBUG(aJet);
          if(aJet == 0) continue;

          TLorentzVector* hlv = new TLorentzVector; DEBUG(hlv);
          aJet->GetMom(*hlv); DEBUG(hlv);
          CLHEP::HepLorentzVector theHLV = CLHEP::HepLorentzVector(hlv->X(), hlv->Y(), hlv->Z(), hlv->T());
          px = boost(theHLV, parentJet).px();
          py = boost(theHLV, parentJet).py();
          pz = boost(theHLV, parentJet).pz();
          fillTensors(px, py, pz);
          fillMij(px, py, pz);
          _Mij /= _sumI;
          _Nij /= _sumJ;

          double pt = boost(theHLV, parentJet).perp();
          if(pt > _maxpT) _maxpT = pt;

          threeMomenta.push_back(boost(theHLV, parentJet).v());

        // Bwaa-ha-ha! ROOT garbage collection make ka-boooooom!
        //delete aJet; 
        aJet = 0; // We're done with this jet; do everything onwards with HLVs
        
        //CLHEP::Hep2Vector qt = CLHEP::Hep2Vector(theHLV.v()); // temp
        CLHEP::Hep3Vector* qt = new CLHEP::Hep3Vector( boost(theHLV, parentJet).v() ); DEBUG(parentJet); DEBUG(theHLV); DEBUG(qt);DEBUG(*qt);
        
        qt->setZ(0); DEBUG(*qt);
        _sumQt += qt->mag(); DEBUG(qt->mag()); DEBUG(_sumQt);
        _qtVect->push_back(qt); DEBUG(_qtVect->size());
        
  //if((*firstJet) == 0) continue; DEBUG(0);
  //const Jet* aJet = (*firstJet)->clone(true, true); DEBUG(0); // Implicit new here
  //if(aJet == 0) continue; DEBUG(0);
  //fourMomenta->push_back(new CLHEP::HepLorentzVector(boost(aJet->hlv(), parentJet))); DEBUG(0);
  //delete aJet; // We're done with this jet; do everything onwards with HLVs
      }     DEBUG("I am here");

      _nConstituents = fillVectors(); DEBUG(_nConstituents);
      calcThrust(threeMomenta); DEBUG(threeMomenta.size());
      _aco = calcAcoplanarity(threeMomenta); DEBUG(threeMomenta.size());
 //      JetCollectionHelper::jetcollection_t theJets(_theJetCollection->begin(), _theJetCollection->end());
 //      JetCollectionHelper::jetcollection_t::iterator first = theJets.begin();
 //      JetCollectionHelper::jetcollection_t::iterator last  = theJets.end();
 //      _sumB = 0.;
 //      for(; first != last; ++first) {
	// if((*first) == 0) continue;
	// const Jet* aJet = (*first)->clone(true, true);
	// if(aJet == 0) continue;
	// px = boost(aJet->hlv(), parentJet).px();
	// py = boost(aJet->hlv(), parentJet).py();
	// pz = boost(aJet->hlv(), parentJet).pz();
	// fillTensors(px, py, pz);
	// delete aJet;
 //      }
    } DEBUG("tensors filled");
        DEBUG("I am here");
    // Matricies are symmetric:
    _sphericityTensor[1][0] = _sphericityTensor[0][1]; DEBUG(0);
    _sphericityTensor[2][0] = _sphericityTensor[0][2]; DEBUG(0);
    _sphericityTensor[2][1] = _sphericityTensor[1][2]; DEBUG(0);

    _spherocityTensor[1][0] = _spherocityTensor[0][1]; DEBUG(0);
    _spherocityTensor[2][0] = _spherocityTensor[0][2]; DEBUG(0);
    _spherocityTensor[2][1] = _spherocityTensor[1][2]; DEBUG(0);
    DEBUG(_pmag2Sum); DEBUG(_pmagSum); 
    _sphericityTensor /= _pmag2Sum? _pmag2Sum: NaN; DEBUG(0);
    _spherocityTensor /= _pmagSum? _pmagSum: NaN; DEBUG(0);
        DEBUG("I am here");
  }




 // Comparator
  struct Mod2Cmp {
    bool operator()(const CLHEP::Hep3Vector& a, const CLHEP::Hep3Vector& b) {
      return a.mag2() > b.mag2();
    }
    bool operator()(const CLHEP::Hep3Vector*& a, const CLHEP::Hep3Vector*& b) {
      return a->mag2() > b->mag2();
    }
  } mod2Cmp;

  // Do the general case thrust calculation
  void calcT(const std::vector<CLHEP::Hep3Vector>& momenta, double& t, CLHEP::Hep3Vector& taxis) {
    // This function implements the iterative algorithm as described in the
    // Pythia manual. We take eight (four) different starting vectors
    // constructed from the four (three) leading particles to make sure that
    // we don't find a local maximum.
    std::vector<CLHEP::Hep3Vector> p = momenta;
    //assert(p.size() >= 3);
    unsigned int n = 3;
    if (p.size() == 3) n = 3;
    std::vector<CLHEP::Hep3Vector> tvec;
    std::vector<double> tval;
    std::sort(p.begin(), p.end(), mod2Cmp);
    for(int i = 0 ; i < pow(2, n-1); ++i) {
      // Create an initial vector from the leading four jets
      CLHEP::Hep3Vector foo(0,0,0);
      int sign = i;
      for(unsigned int k = 0 ; k < n ; ++k) {
        (sign % 2) == 1 ? foo += p[k] : foo -= p[k];
        sign /= 2;
      }
      foo=foo.unit();

      // Iterate
      double diff=999.;
      while(diff>1e-5) {
        CLHEP::Hep3Vector foobar(0,0,0);
        for (unsigned int k=0 ; k<p.size() ; k++)
          foo.dot(p[k])>0 ? foobar+=p[k] : foobar-=p[k];
        diff=(foo-foobar.unit()).mag();
        foo=foobar.unit();
      }

      // Calculate the thrust value for the vector we found
      t=0.;
      for(unsigned int k=0 ; k<p.size() ; k++)
        t+=fabs(foo.dot(p[k]));

      // Store everything
      tval.push_back(t);
      tvec.push_back(foo);
    }

    // Pick the solution with the largest thrust
    t=0.;
    for (unsigned int i=0 ; i<tvec.size() ; i++)
      if (tval[i]>t){
        t=tval[i];
        taxis=tvec[i];
      }
  }



  // Do the full calculation
  void calcThrust(const std::vector<CLHEP::Hep3Vector>& fsmomenta) {
    // Make a vector of the three-momenta in the final state
    double momentumSum(0.0);
      for(int i = 0; i != fsmomenta.size(); i++) {
      momentumSum += (fsmomenta.at(i)).mag();
    }
    //MSG_DEBUG("Number of particles = " << fsmomenta.size());


    // Clear the caches
    _thrusts->clear();
    _thrustAxes->clear();


    // If there are fewer than 2 visible particles, we can't do much
    if (fsmomenta.size() < 2) {
      for (int i = 0; i < 3; ++i) {
        _thrusts->push_back(-1);
        _thrustAxes->push_back(new CLHEP::Hep3Vector(0,0,0));
      }
      return;
    }


    // Handle special case of thrust = 1 if there are only 2 particles
    if (fsmomenta.size() == 2) {
      CLHEP::Hep3Vector axis(0,0,0);
      _thrusts->push_back(1.0);
      _thrusts->push_back(0.0);
      _thrusts->push_back(0.0);
      axis = fsmomenta[0].unit();
      if (axis.z() < 0) axis = -axis;
      _thrustAxes->push_back(new CLHEP::Hep3Vector(axis));
      /// @todo Improve this --- special directions bad...
      /// (a,b,c) _|_ 1/(a^2+b^2) (b,-a,0) etc., but which combination minimises error?
      if (axis.z() < 0.75)
        _thrustAxes->push_back( new CLHEP::Hep3Vector((axis.cross(CLHEP::Hep3Vector(0,0,1))).unit()) );
      else
        _thrustAxes->push_back( new CLHEP::Hep3Vector((axis.cross(CLHEP::Hep3Vector(0,1,0))).unit()) );
      _thrustAxes->push_back( new CLHEP::Hep3Vector((*_thrustAxes)[0]->cross(*(*_thrustAxes)[1])) );
      return;
    }



    // Temporary variables for calcs
    CLHEP::Hep3Vector axis(0,0,0);
    double val = 0.;

    // Get thrust
    calcT(fsmomenta, val, axis);
    //MSG_DEBUG("Mom sum = " << momentumSum);
    _thrusts->push_back(val / momentumSum);
    // Make sure that thrust always points along the +ve z-axis.
    if (axis.z() < 0) axis = -axis;
    axis = axis.unit();
    //MSG_DEBUG("Axis = " << axis);
    _thrustAxes->push_back(new CLHEP::Hep3Vector(axis));

    // Get thrust major
    std::vector<CLHEP::Hep3Vector> threeMomenta;

      for(int i = 0; i != fsmomenta.size(); i++) {
      // Get the part of each 3-momentum which is perpendicular to the thrust axis
      const CLHEP::Hep3Vector vpar = fsmomenta.at(i).dot(axis.unit()) * axis.unit();
      threeMomenta.push_back(fsmomenta.at(i) - vpar);
    }
    calcT(threeMomenta, val, axis);
    _thrusts->push_back(val / momentumSum);
    if (axis.x() < 0) axis = -axis;
    axis = axis.unit();
    _thrustAxes->push_back(new CLHEP::Hep3Vector(axis));

    // Get thrust minor
    if ((*_thrustAxes)[0]->dot(*(*_thrustAxes)[1]) < 1e-10) {
      axis = (*_thrustAxes)[0]->cross(*(*_thrustAxes)[1]);
      _thrustAxes->push_back(new CLHEP::Hep3Vector(axis));
      val = 0.0;
      
        for(int i = 0; i != fsmomenta.size(); i++) {
        val += fabs(axis.dot(fsmomenta.at(i)));
      }
      _thrusts->push_back(val / momentumSum);
    } else {
      _thrusts->push_back(-1.0);
      _thrustAxes->push_back(new CLHEP::Hep3Vector(0,0,0));
    }

  }





int fillVectors() { DEBUG("begin fillVectors"); DEBUG(_instName);
  int N = _qtVect->size(); DEBUG(N);
  _qtCrossN->clear(); DEBUG(0);
  _qtDotN->clear(); DEBUG(0);

  for(int j = 0; j != N; j++) { DEBUG(j);
    double sumQtCrossN = 0.; DEBUG(0);
    double sumQtDotN = 0.; DEBUG(0);

    for(int i = 0; i != N; i++) {DEBUG(i);
      //if(i == j) continue; // forbid auto-correlations?
      sumQtCrossN += fabs( _qtVect->at(i)->cross( _qtVect->at(j)->unit() ).mag() ); DEBUG(sumQtCrossN);
      sumQtDotN += fabs( _qtVect->at(i)->dot( _qtVect->at(j)->unit() ) ); DEBUG(sumQtDotN);
    }
    _qtCrossN->push_back(sumQtCrossN); DEBUG(sumQtCrossN);
    _qtDotN->push_back(sumQtDotN); DEBUG(sumQtDotN);
  }
  std::sort(_qtCrossN->begin(), _qtCrossN->end()); DEBUG("sort1");
  std::sort(_qtDotN->begin(), _qtDotN->end()); DEBUG("sort2");
  
  if(_qtCrossN) {
    if(!_qtCrossN->empty()) {
      _minQtCrossN = _qtCrossN->at(0); DEBUG(_qtCrossN->size()); DEBUG(_minQtCrossN);
      _maxQtCrossN = _qtCrossN->at(N?N-1:0);
    }
  }

  if(_qtDotN) {
    if(!_qtDotN->empty()) {
      _minQtDotN = _qtDotN->at(0); DEBUG(_qtDotN->size()); DEBUG(_minQtDotN);
      _maxQtDotN = _qtDotN->at(N?N-1:0);
    }
  }

  // _maxQtCrossN = _qtCrossN->at(N?N-1:0);
  // _minQtDotN = _qtDotN->at(0);
  // _maxQtDotN = _qtDotN->at(N?N-1:0);

  // _minQtCrossN = *std::min_element(_qtCrossN->begin(), _qtCrossN->end()); DEBUG(_minQtCrossN);
  // _maxQtCrossN = *std::max_element(_qtCrossN->begin(), _qtCrossN->end());
  // _minQtDotN = *std::min_element(_qtDotN->begin(), _qtDotN->end());
  // _maxQtDotN = *std::max_element(_qtDotN->begin(), _qtDotN->end());
  DEBUG("end fillVectors"); DEBUG(_instName);
  return N;
  }


  void fillTensors(double px, double py, double pz) { DEBUG("fill tensors");
    _pmag2 = (px * px) + (py * py) + (pz * pz); DEBUG(_pmag2);
    _pmag = Sqrt(_pmag2); DEBUG(_pmag);
    _invpmag = Divide(1., _pmag); DEBUG(_invpmag);
    _pmagSum += _pmag; DEBUG(_pmagSum);
    _pmag2Sum += _pmag2; DEBUG(_pmag2Sum);

    _pT = sqrt((px * px) + (py * py)); DEBUG(_pT);
    _pT2 = _pT * _pT; DEBUG(_pT2);
    _pTsum += _pT; DEBUG(_pTsum);
    _pT2sum += _pT2; DEBUG(_pT2sum);
    
    _sphericityTensor[0][0] += px * px; DEBUG(_sphericityTensor[0][0]);
    _sphericityTensor[0][1] += px * py; DEBUG(_sphericityTensor[0][1]);
    _sphericityTensor[0][2] += px * pz; DEBUG(_sphericityTensor[0][2]);
    _sphericityTensor[1][1] += py * py; DEBUG(_sphericityTensor[1][1]);
    _sphericityTensor[1][2] += py * pz; DEBUG(_sphericityTensor[1][2]);
    _sphericityTensor[2][2] += pz * pz; DEBUG(_sphericityTensor[2][2]);
    
    _spherocityTensor[0][0] += px * px * _invpmag; DEBUG(_spherocityTensor[0][0]);
    _spherocityTensor[0][1] += px * py * _invpmag; DEBUG(_spherocityTensor[0][1]);
    _spherocityTensor[0][2] += px * pz * _invpmag; DEBUG(_spherocityTensor[0][2]);
    _spherocityTensor[1][1] += py * py * _invpmag; DEBUG(_spherocityTensor[1][1]);
    _spherocityTensor[1][2] += py * pz * _invpmag; DEBUG(_spherocityTensor[1][2]);
    _spherocityTensor[2][2] += pz * pz * _invpmag; DEBUG(_spherocityTensor[2][2]);

    double pt = Sqrt((px * px) + (py * py)); DEBUG(pt);
    _Mlin[0][0] += Divide(px * px, pt); DEBUG(_Mlin[0][0]);
    _Mlin[0][1] += Divide(px * py, pt); DEBUG(_Mlin[0][1]);
    _Mlin[1][0] += Divide(px * py, pt); DEBUG(_Mlin[1][0]);
    _Mlin[1][1] += Divide(py * py, pt); DEBUG(_Mlin[1][1]);

    double eta = fabs(pseudoRapidity(_pmag, pz)); DEBUG(eta);
    _sumB += pt * exp(-1.* pt * eta);
  }

  void fillMij(double px, double py, double pz, double jetEta = 0., double jetPhi = 0., double q = 0.) {
    double pt = Sqrt((px * px) + (py * py)); DEBUG(pt);
    double pmag2 = (px * px) + (py * py) + (pz * pz); DEBUG(pmag2);
    double pmag = Sqrt(pmag2); DEBUG(pmag);
    double eta = atanh(Divide(pz, pmag)); DEBUG(eta);// http://en.wikipedia.org/wiki/Pseudorapidity
    double phi = PhiCorr()(Arctan2(py, px));
    double x = eta - jetEta; double y = DeltaPhi()(phi, jetPhi);
    double Ixy = pt;
    _Mij[0][0] += intpow<0>(x) * intpow<0>(y) * Ixy;
    _Mij[0][1] += intpow<0>(x) * intpow<1>(y) * Ixy;
    _Mij[0][2] += intpow<0>(x) * intpow<2>(y) * Ixy;
    _Mij[0][3] += intpow<0>(x) * intpow<3>(y) * Ixy;
    _Mij[1][0] += intpow<1>(x) * intpow<0>(y) * Ixy;
    _Mij[1][1] += intpow<1>(x) * intpow<1>(y) * Ixy;
    _Mij[1][2] += intpow<1>(x) * intpow<2>(y) * Ixy;
    _Mij[1][3] += intpow<1>(x) * intpow<3>(y) * Ixy;
    _Mij[2][0] += intpow<2>(x) * intpow<0>(y) * Ixy;
    _Mij[2][1] += intpow<2>(x) * intpow<1>(y) * Ixy;
    _Mij[2][2] += intpow<2>(x) * intpow<2>(y) * Ixy;
    _Mij[2][3] += intpow<2>(x) * intpow<3>(y) * Ixy;
    _Mij[3][0] += intpow<3>(x) * intpow<0>(y) * Ixy;
    _Mij[3][1] += intpow<3>(x) * intpow<1>(y) * Ixy;
    _Mij[3][2] += intpow<3>(x) * intpow<2>(y) * Ixy;
    _Mij[3][3] += intpow<3>(x) * intpow<3>(y) * Ixy;
    _sumI += Ixy;

    double Jxy = q*pt;
    _Nij[0][0] += intpow<0>(x) * intpow<0>(y) * Jxy;
    _Nij[0][1] += intpow<0>(x) * intpow<1>(y) * Jxy;
    _Nij[0][2] += intpow<0>(x) * intpow<2>(y) * Jxy;
    _Nij[0][3] += intpow<0>(x) * intpow<3>(y) * Jxy;
    _Nij[1][0] += intpow<1>(x) * intpow<0>(y) * Jxy;
    _Nij[1][1] += intpow<1>(x) * intpow<1>(y) * Jxy;
    _Nij[1][2] += intpow<1>(x) * intpow<2>(y) * Jxy;
    _Nij[1][3] += intpow<1>(x) * intpow<3>(y) * Jxy;
    _Nij[2][0] += intpow<2>(x) * intpow<0>(y) * Jxy;
    _Nij[2][1] += intpow<2>(x) * intpow<1>(y) * Jxy;
    _Nij[2][2] += intpow<2>(x) * intpow<2>(y) * Jxy;
    _Nij[2][3] += intpow<2>(x) * intpow<3>(y) * Jxy;
    _Nij[3][0] += intpow<3>(x) * intpow<0>(y) * Jxy;
    _Nij[3][1] += intpow<3>(x) * intpow<1>(y) * Jxy;
    _Nij[3][2] += intpow<3>(x) * intpow<2>(y) * Jxy;
    _Nij[3][3] += intpow<3>(x) * intpow<3>(y) * Jxy;
    _sumJ += Jxy;
  }

void calcHu() {
  double xbar = _Mij[1][0]/_Mij[0][0];
  double ybar = _Mij[0][1]/_Mij[0][0];
  _mu[0][0] = _Mij[0][0];
  _mu[0][1] = 0;
  _mu[1][0] = 0;
  _mu[1][1] = _Mij[1][1] - (xbar * _Mij[0][1]);
  _mu[2][0] = _Mij[2][0] - (xbar * _Mij[1][0]);
  _mu[0][2] = _Mij[0][2] - (ybar * _Mij[0][1]);
  _mu[2][1] = _Mij[2][1] - (2. * xbar * _Mij[1][1]) - (ybar * _Mij[2][0]) + (2. * xbar * xbar * _Mij[0][1]);
  _mu[1][2] = _Mij[1][2] - (2. * ybar * _Mij[1][1]) - (xbar * _Mij[0][2]) + (2. * ybar * ybar * _Mij[1][0]);
  _mu[3][0] = _Mij[3][0] - (3. * xbar * _Mij[2][0]) + (2. * xbar * xbar * _Mij[1][0]);
  _mu[0][3] = _Mij[0][3] - (3. * ybar * _Mij[0][2]) + (2. * ybar * ybar * _Mij[0][1]);
  for(unsigned int i = 0; i != 4; i++) {
    for(unsigned int j = 0; j != 4; j++) {
      _eta[i][j] = Divide(_mu[i][j], pow(_mu[0][0], (1 + ((i + j) / 2))));
    }
  }

  CLHEP::HepMatrix muP = CLHEP::HepMatrix(3, 3, 0.);
  muP[2][0] = _mu[2][0] / _mu[0][0];
  muP[0][2] = _mu[0][2] / _mu[0][0];
  muP[1][1] = _mu[1][1] / _mu[0][0];

double lambdaPlus = (0.5 * (muP[2][0] - muP[0][2])) + (0.5 * sqrt( (4. * intpow<2>(muP[1][1])) + intpow<2>(muP[2][0] - muP[0][2])));
double lambdaMinus = (0.5 * (muP[2][0] - muP[0][2])) - (0.5 * sqrt( (4. * intpow<2>(muP[1][1])) + intpow<2>(muP[2][0] - muP[0][2])));
double lambda1 = std::max(fabs(lambdaPlus), fabs(lambdaMinus));
double lambda2 = std::min(fabs(lambdaPlus), fabs(lambdaMinus));


  _ThetaHu = 0.5 * atan2((2. * muP[1][1]), (muP[2][0] - muP[0][2]));
  _eHu = sqrt(1. - Divide(lambda2, lambda1));

  _I1 = _eta[2][0] + _eta[0][2];

  _I2 = intpow<2>(_eta[2][0] - _eta[0][2]) + (4. * intpow<2>(_eta[1][1]));

  _I3 = intpow<2>(_eta[3][0] - (3. * _eta[1][2])) + intpow<2>((3. * _eta[2][1]) - _eta[0][3]);
  
  _I4 = intpow<2>(_eta[3][0] + _eta[1][2]) + intpow<2>(_eta[2][1] + _eta[0][3]);

  _I5 = ( (_eta[3][0] - (3. * _eta[1][2])) * (_eta[3][0] + _eta[1][2]) * ( intpow<2>(_eta[3][0] + _eta[1][2]) - (3. * intpow<2>(_eta[2][1] + _eta[0][3])) ) ) + 
  ( ( (3. * _eta[2][1]) - _eta[0][3]) * (_eta[2][1] + _eta[0][3]) * ( (3. * intpow<2>(_eta[3][0] + _eta[1][2])) - intpow<2>(_eta[2][1] + _eta[0][3]) ) );
  
  _I6 = ( (_eta[2][0] - _eta[0][2]) * ( intpow<2>(_eta[3][0] + _eta[1][2]) - intpow<2>(_eta[2][1] + _eta[0][3]) ) ) +
  (4. * _eta[1][1] * (_eta[3][0] + _eta[1][2]) * (_eta[2][1] + _eta[0][3]));

  _I7 = ( ( (3. * _eta[2][1]) - _eta[0][3]) * (_eta[3][0] + _eta[1][2]) * ( intpow<2>(_eta[3][0] + _eta[1][2]) - (3. * intpow<2>(_eta[2][1] + _eta[0][3])) ) ) - 
  ( (_eta[3][0] - (3. * _eta[1][2])) * (_eta[2][1] + _eta[0][3]) * ( (3. * intpow<2>(_eta[3][0] + _eta[1][2])) - intpow<2>(_eta[2][1] + _eta[0][3]) ) );
  
  _I8 = (_eta[1][1] * ( intpow<2>(_eta[3][0] + _eta[1][2]) - intpow<2>(_eta[0][3] + _eta[2][1]) ) ) -
  ( (_eta[2][0] - _eta[0][2]) * (_eta[3][0] + _eta[1][2]) * (_eta[0][3] + _eta[2][1]) );

  }

  void calcHu2() {
  double xtilde = _Nij[1][0]/_Nij[0][0];
  double ytilde = _Nij[0][1]/_Nij[0][0];
  _nu[0][0] = _Nij[0][0];
  _nu[0][1] = 0;
  _nu[1][0] = 0;
  _nu[1][1] = _Nij[1][1] - (xtilde * _Nij[0][1]);
  _nu[2][0] = _Nij[2][0] - (xtilde * _Nij[1][0]);
  _nu[0][2] = _Nij[0][2] - (ytilde * _Nij[0][1]);
  _nu[2][1] = _Nij[2][1] - (2. * xtilde * _Nij[1][1]) - (ytilde * _Nij[2][0]) + (2. * xtilde * xtilde * _Nij[0][1]);
  _nu[1][2] = _Nij[1][2] - (2. * ytilde * _Nij[1][1]) - (xtilde * _Nij[0][2]) + (2. * ytilde * ytilde * _Nij[1][0]);
  _nu[3][0] = _Nij[3][0] - (3. * xtilde * _Nij[2][0]) + (2. * xtilde * xtilde * _Nij[1][0]);
  _nu[0][3] = _Nij[0][3] - (3. * ytilde * _Nij[0][2]) + (2. * ytilde * ytilde * _Nij[0][1]);
  for(unsigned int i = 0; i != 4; i++) {
    for(unsigned int j = 0; j != 4; j++) {
      _zeta[i][j] = Divide(_nu[i][j], pow(_nu[0][0], (1 + ((i + j) / 2))));
    }
  }

  CLHEP::HepMatrix nuP = CLHEP::HepMatrix(3, 3, 0.);
  nuP[2][0] = _nu[2][0] / _nu[0][0];
  nuP[0][2] = _nu[0][2] / _nu[0][0];
  nuP[1][1] = _nu[1][1] / _nu[0][0];

double lambdaPlus = (0.5 * (nuP[2][0] - nuP[0][2])) + (0.5 * sqrt( (4. * intpow<2>(nuP[1][1])) + intpow<2>(nuP[2][0] - nuP[0][2])));
double lambdaMinus = (0.5 * (nuP[2][0] - nuP[0][2])) - (0.5 * sqrt( (4. * intpow<2>(nuP[1][1])) + intpow<2>(nuP[2][0] - nuP[0][2])));
double lambda1 = std::max(lambdaPlus, lambdaMinus);
double lambda2 = std::min(lambdaPlus, lambdaMinus);


  _ThetaHu2 = 0.5 * atan2((2. * nuP[1][1]), (nuP[2][0] - nuP[0][2]));
  _eHu2 = sqrt(1. - Divide(lambda2, lambda1));

  _J1 = _zeta[2][0] + _zeta[0][2];

  _J2 = intpow<2>(_zeta[2][0] - _zeta[0][2]) + (4. * intpow<2>(_zeta[1][1]));

  _J3 = intpow<2>(_zeta[3][0] - (3. * _zeta[1][2])) + intpow<2>((3. * _zeta[2][1]) - _zeta[0][3]);
  
  _J4 = intpow<2>(_zeta[3][0] + _zeta[1][2]) + intpow<2>(_zeta[2][1] + _zeta[0][3]);

  _J5 = ( (_zeta[3][0] - (3. * _zeta[1][2])) * (_zeta[3][0] + _zeta[1][2]) * ( intpow<2>(_zeta[3][0] + _zeta[1][2]) - (3. * intpow<2>(_zeta[2][1] + _zeta[0][3])) ) ) + 
  ( ( (3. * _zeta[2][1]) - _zeta[0][3]) * (_zeta[2][1] + _zeta[0][3]) * ( (3. * intpow<2>(_zeta[3][0] + _zeta[1][2])) - intpow<2>(_zeta[2][1] + _zeta[0][3]) ) );
  
  _J6 = ( (_zeta[2][0] - _zeta[0][2]) * ( intpow<2>(_zeta[3][0] + _zeta[1][2]) - intpow<2>(_zeta[2][1] + _zeta[0][3]) ) ) +
  (4. * _zeta[1][1] * (_zeta[3][0] + _zeta[1][2]) * (_zeta[2][1] + _zeta[0][3]));

  _J7 = ( ( (3. * _zeta[2][1]) - _zeta[0][3]) * (_zeta[3][0] + _zeta[1][2]) * ( intpow<2>(_zeta[3][0] + _zeta[1][2]) - (3. * intpow<2>(_zeta[2][1] + _zeta[0][3])) ) ) - 
  ( (_zeta[3][0] - (3. * _zeta[1][2])) * (_zeta[2][1] + _zeta[0][3]) * ( (3. * intpow<2>(_zeta[3][0] + _zeta[1][2])) - intpow<2>(_zeta[2][1] + _zeta[0][3]) ) );
  
  _J8 = (_zeta[1][1] * ( intpow<2>(_zeta[3][0] + _zeta[1][2]) - intpow<2>(_zeta[0][3] + _zeta[2][1]) ) ) -
  ( (_zeta[2][0] - _zeta[0][2]) * (_zeta[3][0] + _zeta[1][2]) * (_zeta[0][3] + _zeta[2][1]) );

  }

  const CLHEP::HepMatrix* sphericityTensor() const {
    return _theJetCollection? &_sphericityTensor: 0;
  }

  const CLHEP::HepMatrix* spherocityTensor() const {
    return _theJetCollection? &_spherocityTensor: 0;
  }

    const CLHEP::HepMatrix* sphericityTensor(int /* dummy */) const {
    return _fTracks? &_sphericityTensor: 0;
  }

  const CLHEP::HepMatrix* spherocityTensor(int /* dummy */) const {
    return _fTracks? &_spherocityTensor: 0;
  }

  bool normalisedEigenvalues(const CLHEP::HepMatrix* M, double* root) const { DEBUG(0); DEBUG(M); 
    EigenSolverSym3x3 solve(*M); DEBUG(0);
    solve(root); DEBUG(0);
    double sum = NaN;DEBUG(0);
    if(!isnan(root[0]) && !isnan(root[1]) && !isnan(root[2])) sum = root[0] + root[1] + root[2];
    double norm = Divide(1., sum); DEBUG(0);
    root[0] = Multiply(root[0], norm); DEBUG(0);
    root[1] = Multiply(root[1], norm); DEBUG(0);
    root[2] = Multiply(root[2], norm); DEBUG(0);
    return sum != 0.;
  }

double calcAcoplanarity(const std::vector<CLHEP::Hep3Vector>& fsmomenta) {
  if(fsmomenta.size() <= 1) return NaN; // too few tracks to do anything
  int nperm = fsmomenta.size(); DEBUG(nperm);
  CLHEP::Hep3Vector akovec = CLHEP::Hep3Vector(0.,0.,0.);
  // loop over permutations, find acoplanarity axis
  double akopl = 99999999999.; // or arbitrary very big number, eg, 99999.
  for(int i = 0; i != nperm; i++) {
    for(int j = 0; j != nperm; j++) {
      if(i == j) continue;
      DEBUG(fsmomenta.at(i));
      DEBUG(fsmomenta.at(j));
  // trial acoplanarity plane
      CLHEP::Hep3Vector pcros = fsmomenta.at(i).cross(fsmomenta.at(j)); DEBUG(pcros);
      double ptot = pcros.mag(); DEBUG(ptot);
      if(ptot == 0.) continue;
      double aktry = 0.;
      for(int k = 0; k != nperm; k++) { DEBUG(fsmomenta.at(k));
        aktry += pcros.dot(fsmomenta.at(k));
      } // loop over k
      aktry /= ptot; DEBUG(aktry);
      if(aktry < akopl) {
        akovec = pcros / ptot; DEBUG(akovec);
        akopl = aktry; DEBUG(akopl);
      } // minimisation test 
      //DEBUG(aktry);
  }// loop over j
  } // loop over i
  double psum = 0.;
  akopl = 0.; DEBUG(akovec);
  for(int i = 0; i != nperm; i++) {
    akopl += fabs(akovec.dot(fsmomenta.at(i)));
    psum += fsmomenta.at(i).mag();
  } // loop over i
  double q = (akopl / psum); DEBUG(psum);
  akopl = 4. * q * q; DEBUG(akopl);
  return akopl;
} // end of function

  const TClonesArray* _theJetCollection;
  std::vector<CLHEP::Hep3Vector*>* _qtVect;
  std::vector<CLHEP::Hep3Vector*>* _3momenta;
  std::vector<double>* _thrusts;
  std::vector<CLHEP::Hep3Vector*>* _thrustAxes;
  std::vector<double>* _qtDotN;
  std::vector<double>* _qtCrossN;
  CLHEP::HepMatrix _sphericityTensor;
  CLHEP::HepMatrix _spherocityTensor;
  CLHEP::HepMatrix _Mij;
  CLHEP::HepMatrix _Nij;
  CLHEP::HepMatrix _mu;
  CLHEP::HepMatrix _nu;
  CLHEP::HepMatrix _eta;
  CLHEP::HepMatrix _zeta;
  double _ThetaHu, _eHu, _I1, _I2, _I3, _I4, _I5, _I6, _I7, _I8;
  double _ThetaHu2, _eHu2, _J1, _J2, _J3, _J4, _J5, _J6, _J7, _J8;
  Hep2Matrix _Mlin;
  double _aco;
  double _pmag;
  double _pmag2;
  double _pmagSum;
  double _pmag2Sum;
  double _pT;
  double _pT2;
  double _pTsum;
  double _pT2sum;
  double _sumI;
  double _sumJ;
  double _sumB;
  double _sumQt;
  double _maxpT;
  int _nConstituents; // tracks, jets, subjets, cemins, bemins, memims, etc...
  double _invpmag;
  double _sphericityLambda1, _sphericityLambda2, _sphericityLambda3;
  double _spherocityLambda1, _spherocityLambda2, _spherocityLambda3;
  double _minQtCrossN, _maxQtCrossN, _minQtDotN, _maxQtDotN;
  TClonesArray* _fTracks;
  char* _instName;
};

// Ancient fortran code for calculating acoplanarity:
// From: http://www.fe.infn.it/~mandreot/Didattica/VitaMuoni/cern/links/pro/px114.car
// (OPAL PX library, circa 1989)
// *.*********************************************************
// *. ------
// *. PXAKO4
// *. ------
// *. Routine to calculate the acoplanarity distribution.
// *. The algorithm is an approximate one based on a Jade
// *. program by S.Bethke and E.Elsen.
// *. The algorithm has the following characteristics:
// *.   (1) The missing momentum is added as an extra particle
// *.       in order to "complete" momentum conservation.
// *.   (2) The acoplanarity plane is assumed to be defined by
// *.       the momentum vectors of two of the particles in the
// *.       event. The plane so-defined which results in the
// *.       smallest acoplanarity value is taken as the "true"
// *.       acoplanarity plane.
// *.   (3) Should the number of particles in the event exceed
// *.       the cutoff NTPER (see explanation of the argument list
// *.       below) only the NTPER particles having the largest
// *.       momentum are used to define the acoplanarity plane.
// *.       The missing momentum vector used to complete momentum
// *.       conservation is calculated using only these NTPER particles.
// *.   (4) The final acoplanarity value is calculated using all
// *.       particles, whether or not they were used in the search
// *.       for the acoplanarity plane.  This acoplanarity calculation
// *.       also includes a fictitious particle introduced to complete
// *.       momentum conservation, should one be necessary.
// *. Usage     :
// *.
// *.      INTEGER  ITKDM,MXTRK,NTPER
// *.      PARAMETER  (ITKDM=4.or.more,MXTRK=1.or.more)
// *.      PARAMETER  (NTPER=80.or.so)
// *.      INTEGER  NTRAK,IERR
// *.      REAL  PTRAK (ITKDM,MXTRK),AKOVEC (3.or.more)
// *.      REAL  AKOPL
// *.
// *.      NTRAK = 1.to.MXTRK
// *.      CALL PXAKO4 (NTRAK,ITKDM,MXTRK,PTRAK,
// *.     +            NTPER,AKOPL,AKOVEC,IERR)
// *.
// *. INPUT     : NTRAK   Total number of particles
// *. INPUT     : ITKDM   First dimension of PTRAK
// *. INPUT     : MXTRK   Maximum number particles allowed in PTRAK
// *. INPUT     : PTRAK   Particle 4-momentum (Px,Py,Pz,E)
// *. INPUT     : NTPER   Maximum number of particles to use for the
// *.                     Acoplanarity calculation (for speed purposes)
// *. OUTPUT    : AKOPL   The acoplanarity value
// *. OUTPUT    : AKOVEC  the normalized acoplanarity axis vector
// *. OUTPUT    : IERR    = 0 if all is OK
// *.
// *. CALLS     : PXADDV,PXMAGV,PXCOPV,PXSORV,PXCRO3,PXDOTV
// *. CALLED    : By User
// *.
// *. AUTHOR    :  M.Weber/J.W.Gary
// *. CREATED   :  05-Mar-89
// *. LAST MOD  :  05-Mar-89
// *.
// *. Modification Log.
// *.
// *.*********************************************************
// +SEQ,DECLARE.
//       INTEGER LMXTK
//       PARAMETER  (LMXTK=250)
//       INTEGER  KIX (LMXTK)
//       INTEGER  NTRAK,ITKDM,MXTRK,IERR,NTPER,NPERM,IX,IP1,IP2,IP3,
//      +         INX1,INX2,IOFF,IKCNT
//       REAL  PTRAK (ITKDM,*),AKOVEC (*),PCROS (3),PMISS (4),
//      +      PMAG (LMXTK)
//       REAL  AKOPL,PTOT,AKTRY,DOTP,PSUM
//       SAVE  IKCNT
//       LOGICAL  SORT
//       DATA  IKCNT / 0 /

//       IERR = 0
//       IF (NTRAK.LE.1.OR.NTRAK.GE.MXTRK) THEN
//           WRITE (6,FMT='('' PXAKO4: Error, NTRAK+1 ='',I6,
//      +           '' greater than MXTRK ='',I6)') NTRAK+1,MXTRK
//           IERR = -1
//           GO TO 990
//       END IF
//       IF (NTRAK.GE.LMXTK) THEN
//           WRITE (6,FMT='('' PXAKO4: Error, NTRAK+1 ='',I6,
//      +           '' greater than LMXTK ='',I6)') NTRAK+1,LMXTK
//           IERR = -1
//           GO TO 990
//       END IF
// *  sort momenta if number of tracks too large
// *  ---- ------- -- ------ -- ------ --- -----
//       NPERM = NTRAK
//       IF (NPERM.GT.NTPER) THEN
//           DO 110 IP1 = 1,NTRAK
//               CALL PXMAGV (3,PTRAK (1,IP1),PMAG (IP1))
//  110      CONTINUE
//           CALL PXSORV (NTRAK,PMAG,KIX,'I')
//           NPERM = NTPER
//           SORT = .TRUE.
//           IKCNT = IKCNT + 1
//           IF (IKCNT.LE.5) WRITE (6,FMT='(
//      +       '' PXKOP4: NTRAK,NPERM ='',2I10)') NTRAK,NPERM
//       ELSE
//           DO 115 IP1 = 1,NTRAK
//               KIX (IP1) = IP1
//  115      CONTINUE
//           SORT = .FALSE.
//       END IF
// *  missing momentum for particles to be used to find plane
// *  ------- -------- --- --------- -- -- ---- -- ---- -----
//       CALL PXZERV (4,PMISS)
//       DO 120 IP1 = 1,NPERM
// C >>> This loop for Bethke/Elsen (almost)      DO 120 IP1 = 1,NTRAK
//           INX1 = KIX ((NTRAK+1)-IP1)
//           CALL PXADDV (4,PTRAK (1,INX1),PMISS,PMISS)
//  120  CONTINUE
//       CALL PXMAGV (3,PMISS,PTOT)
//       IOFF = 1
//       IF (PTOT.GT.(1.E-4)) THEN
//           NPERM = NPERM + 1
//           IOFF = 2
//           DO 125 IX = 1,3
//               PTRAK (IX,NTRAK+1) = - PMISS (IX)
//  125      CONTINUE
//           KIX (NTRAK+1) = NTRAK + 1
//       END IF
// *  loop over permutations, find acoplanarity axis
// *  ---- ---- ------------  ---- ------------ ----
//       AKOPL = 999999.
//       DO 170 IP1 = 1,NPERM-1
//           DO 160 IP2 = IP1+1,NPERM
// *           trial acoplanarity plane
// *           ----- ------------ -----
//               INX1 = KIX ((NTRAK+IOFF)-IP1)
//               INX2 = KIX ((NTRAK+IOFF)-IP2)
//               CALL PXCRO3 (PTRAK (1,INX1),PTRAK (1,INX2),PCROS)
//               CALL PXMAGV (3,PCROS,PTOT)
//               IF (PTOT.EQ.0) THEN
//                   WRITE (6,FMT='('' PXAKO4: PTOT ='',E12.4)') PTOT
//                   GO TO 160
//               END IF
// *           trial acoplanarity value
// *           ----- ------------ -----
//               AKTRY = 0.
//               DO 140 IP3 = 1,NPERM
//                   INX1 = KIX ((NTRAK+IOFF)-IP3)
//                   CALL PXDOTV (3,PCROS,PTRAK (1,INX1),DOTP)
//                   AKTRY = AKTRY + ABS (DOTP)
//  140          CONTINUE
//               AKTRY = AKTRY / PTOT
// *           store new acoplanarity value and axis
// *           ----- --- ------------ ----- --- ----
//               IF (AKTRY.LT.AKOPL) THEN
//                   DO 150 IX = 1,3
//                       AKOVEC (IX) = PCROS (IX) / PTOT
//  150              CONTINUE
//                   AKOPL = AKTRY
//               END IF
//  160      CONTINUE
//  170  CONTINUE
// *  missing momentum for all particles
// *  ------- -------- --- --- ---------
//       IF (SORT) CALL PXMIS4 (NTRAK,ITKDM,PTRAK,PTRAK (1,NTRAK+1))
//       PSUM = 0.
//       AKOPL = 0.
//       DO 180 IP1 = 1,(NTRAK+IOFF)-1
//           CALL PXDOTV (3,AKOVEC,PTRAK (1,IP1),DOTP)
//           AKOPL = AKOPL + ABS (DOTP)
//           CALL PXMAGV (3,PTRAK (1,IP1),PTOT)
//           PSUM = PSUM + PTOT
//  180  CONTINUE
//       AKOPL = 4. * (AKOPL / PSUM)**2

//  990  RETURN
//       END


// Older implementation
// arXiv:1010.3698v1 [hep-ph]
// class SphericitySpherocity {
// public:
//   SphericitySpherocity(const JetCollection* theJetCollection): 
//   _theJetCollection(theJetCollection), _pmag(0), _pmag2(0), _pmagSum(0), _pmag2Sum(0),
//   _sphericityLambda1(NaN), _sphericityLambda2(NaN), _sphericityLambda3(NaN),
//   _spherocityLambda1(NaN), _spherocityLambda2(NaN), _spherocityLambda3(NaN) {
//     _sphericityTensor = CLHEP::HepMatrix(3, 3, 0.);
//     _spherocityTensor = CLHEP::HepMatrix(3, 3, 0.);
    
//     if(_theJetCollection) {
//       compute();
      
//       double lambda[3] = {NaN, NaN, NaN};
//       if(normalisedEigenvalues(sphericityTensor(), lambda)) {
// 	_sphericityLambda1 = lambda[0];
// 	_sphericityLambda2 = lambda[1];
// 	_sphericityLambda3 = lambda[2];
//       }
      
//       if(normalisedEigenvalues(spherocityTensor(), lambda)) {
// 	_spherocityLambda1 = lambda[0];
// 	_spherocityLambda2 = lambda[1];
// 	_spherocityLambda3 = lambda[2];
//       }
//     }
//   }

//   // http://elpub.bib.uni-wuppertal.de/servlets/DerivateServlet/Derivate-1480/dc1006.pdf
//   double circularity() const {
//     return 2.* Divide(std::min(_sphericityLambda1, _sphericityLambda2), _sphericityLambda1 + _sphericityLambda2);
//   }

//   double sphericity() const {
//     return (3./2.) * (_sphericityLambda2 + _sphericityLambda3);
//   }

//   double spherocity() const {
//     return (3./2.) * (_spherocityLambda2 + _spherocityLambda3);
//   }

//   double aplanarity() const {
//     return (3./2.) * _sphericityLambda3;
//   }

//   double aplanority() const {
//     return (3./2.) * _spherocityLambda3;
//   }

//   double Y() const {
//     return (sqrt(3.)/2.) * (_sphericityLambda2 - _sphericityLambda3);
//   }

//   double DShape() const {
//     return 27.* _spherocityLambda1 * _spherocityLambda2 * _spherocityLambda3;
//   }
  
//   double Cshape() const {
//     return 3.* ((_spherocityLambda1 * _spherocityLambda2) + (_spherocityLambda2 * _spherocityLambda3) + (_spherocityLambda3 * _spherocityLambda1));
//   }

// private:
//   void compute() {
//     if(!_theJetCollection || _theJetCollection->size() == 0) {
//       return;
//     }
//     JetCollectionHelper::jetcollection_t theJets(_theJetCollection->begin(), _theJetCollection->end());
//     JetCollectionHelper::jetcollection_t::iterator first = theJets.begin();
//     JetCollectionHelper::jetcollection_t::iterator last  = theJets.end();
    
//     double px, py, pz;

//     for(; first != last; ++first) {
//       if((*first) == 0) continue;
//       const Jet* aJet = (*first)->clone(true, true);
//       if(aJet == 0) continue;

//       px = aJet->px();
//       py = aJet->py();
//       pz = aJet->pz();

//       _pmag2 = (px * px) + (py * py) + (pz * pz);
//       _pmag = Sqrt(_pmag2);
//       _invpmag = Divide(1., _pmag);
//       _pmagSum += _pmag;
//       _pmag2Sum += _pmag2;

//       _sphericityTensor[0][0] += px * px;
//       _sphericityTensor[0][1] += px * py;
//       _sphericityTensor[0][2] += px * pz;
//       _sphericityTensor[1][1] += py * py;
//       _sphericityTensor[1][2] += py * pz;
//       _sphericityTensor[2][2] += pz * pz;

//       _spherocityTensor[0][0] += px * px * _invpmag;
//       _spherocityTensor[0][1] += px * py * _invpmag;
//       _spherocityTensor[0][2] += px * pz * _invpmag;
//       _spherocityTensor[1][1] += py * py * _invpmag;
//       _spherocityTensor[1][2] += py * pz * _invpmag;
//       _spherocityTensor[2][2] += pz * pz * _invpmag;
//       delete aJet;
//     }
//     // Matricies are symmetric:
//     _sphericityTensor[1][0] = _sphericityTensor[0][1];
//     _sphericityTensor[2][0] = _sphericityTensor[0][2];
//     _sphericityTensor[2][1] = _sphericityTensor[1][2];

//     _spherocityTensor[1][0] = _spherocityTensor[0][1];
//     _spherocityTensor[2][0] = _spherocityTensor[0][2];
//     _spherocityTensor[2][1] = _spherocityTensor[1][2];

//     _sphericityTensor /= _pmag2Sum;
//     _spherocityTensor /= _pmagSum;
//   }

//   const CLHEP::HepMatrix* sphericityTensor() const {
//     return _theJetCollection? &_sphericityTensor: 0;
//   }

//   const CLHEP::HepMatrix* spherocityTensor() const {
//     return _theJetCollection? &_spherocityTensor: 0;
//   }

//   bool normalisedEigenvalues(const CLHEP::HepMatrix* M, double* root) const {
//     EigenSolverSym3x3 solve(*M);
//     solve(root);
//     double sum = NaN;
//     if(!isnan(root[0]) && !isnan(root[1]) && !isnan(root[2])) sum = root[0] + root[1] + root[2];
//     double norm = Divide(1., sum);
//     root[0] = Multiply(root[0], norm);
//     root[1] = Multiply(root[1], norm);
//     root[2] = Multiply(root[2], norm);
//     return sum != 0.;
//   }

//   const JetCollection* _theJetCollection;
//   CLHEP::HepMatrix _sphericityTensor;
//   CLHEP::HepMatrix _spherocityTensor;
//   double _pmag;
//   double _pmag2;
//   double _pmagSum;
//   double _pmag2Sum;
//   double _invpmag;
//   double _sphericityLambda1, _sphericityLambda2, _sphericityLambda3;
//   double _spherocityLambda1, _spherocityLambda2, _spherocityLambda3;
// };

#endif // SPHERICITYSPHEROCITY_H


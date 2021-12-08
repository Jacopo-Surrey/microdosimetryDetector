//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Authors: S. Guatelli and J. Davis, CMRP, UOW
// 
//

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "AnalysisManager.hh"
#include "G4Accumulable.hh"
#include <map>

#include "RunMessenger.hh"
#include "DetectorMessenger.hh"

class G4AccumulableManager;

class RunMessenger;

class RunAction : public G4UserRunAction
{
public:
    RunAction(AnalysisManager* analysis, DetectorMessenger* detector);

   ~RunAction();

//public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

	void IncreaseHitCount();
	// To be called once per event when at least a hit is recorded
	// inside the detector.
	// Increments the total number of hits recorded by one,
	// then check if they have reached hitsRequired, and
	// if so stops the simulation.
	
	void setHitsRequired(G4int hits) { hitsRequired = hits; }

private:
	// !!! Change the following depending on the desired statistics !!!
	G4int hitsRequired;

	G4AccumulableManager* accumulableManager;
	
	AnalysisManager* analysisMan;
	DetectorMessenger* detectorMess;

	G4Accumulable<G4int> iHits;

	RunMessenger* messenger;

};
#endif






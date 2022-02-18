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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//

#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "SteppingAction.hh"
#include "AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4StepPoint.hh"
#include "G4RunManager.hh"

#include "DetectorConstruction.hh"

SteppingAction::SteppingAction(AnalysisManager* pAnalysis)
{ 
	analysis = pAnalysis;
	fSecondary = 0;
	
	kinScorerName = "kinScorer_phys";
	
	//activeVolumeName = "SV_phys_1";
	//G4cout << "SteppingAction: outputting sensitive volume: " << activeVolumeName << G4endl;
}

SteppingAction::~SteppingAction()
{ 
}

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
	G4StepPoint* preStepPoint = aStep -> GetPreStepPoint();
	//G4StepPoint* postStepPoint = aStep -> GetPostStepPoint();
	
	G4String preStepVolumeName = preStepPoint -> GetPhysicalVolume() -> GetName();
	//G4String postStepVolumeName = postStepPoint -> GetPhysicalVolume() -> GetName();
	
	// if it's the first (!) step inside the desired (very thin!) volume
	if( (preStepPoint -> GetStepStatus() == fGeomBoundary) && (preStepVolumeName == kinScorerName) )
	{
		// only consider ions, protons, and neutrons for path length
		G4String particleName = aStep -> GetTrack() -> GetParticleDefinition() -> GetParticleName();
		G4int particleZ = aStep -> GetTrack() -> GetParticleDefinition() -> GetAtomicNumber();
		
		if ( (particleZ > 0) || (particleName == "neutron"))
		{	
			G4double kin = preStepPoint -> GetKineticEnergy();
			
			G4int eventID = G4RunManager::GetRunManager() ->
								GetCurrentEvent() -> GetEventID();
			
			analysis -> ScoreKineticEnergy(kin, particleZ, eventID);
		}
	}
	
	/*if( (postStepVolumeName == kinScorerName) && (preStepVolumeName != kinScorerName)  )
	{
		G4double kin = postStepPoint -> GetKineticEnergy();
		
		G4cout << "SCORED KINETIC ENERGY: " << kin << G4endl;
	}*/
	
	/*
  G4SteppingManager*  steppingManager = fpSteppingManager;
  G4Track* theTrack = aStep -> GetTrack();

  // check if it is alive
  if(theTrack-> GetTrackStatus() == fAlive) { return; }

  // Retrieve the secondary particles
    fSecondary = steppingManager -> GetfSecondary();

   for(size_t lp1=0;lp1<(*fSecondary).size(); lp1++)
     { 
       // Retrieve the info about the generation of secondary particles 
       G4String volumeName = (*fSecondary)[lp1] -> GetVolume() -> GetName(); // volume where secondary was generated 
       G4String secondaryParticleName =  (*fSecondary)[lp1]->GetDefinition() -> GetParticleName();  // name of the secondary
       G4double secondaryParticleKineticEnergy =  (*fSecondary)[lp1] -> GetKineticEnergy(); // kinetic energy
      // G4String process = (*fSecondary)[lp1]-> GetCreatorProcess()-> GetProcessName();   // process creating it
       G4double charge = (*fSecondary)[lp1] -> GetDynamicParticle() -> GetDefinition() -> GetPDGCharge();
       G4int AA = (*fSecondary)[lp1] -> GetDynamicParticle() -> GetDefinition() -> GetBaryonNumber();
          
     if (volumeName == activeVolumeName)
	 {
	   if ((secondaryParticleName == "proton") ||
               (secondaryParticleName == "neutron")||
               (secondaryParticleName == "alpha") ||
               (secondaryParticleName == "deuton") || 
               (secondaryParticleName == "triton") || 
               (secondaryParticleName == "He3") || 
	       (secondaryParticleName =="GenericIon"))
               analysis -> FillSecondaries(AA, charge, secondaryParticleKineticEnergy); 	
          }
   }*/
}


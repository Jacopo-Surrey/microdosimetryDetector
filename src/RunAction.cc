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
// Authors: S. Guatelli and J. Davis, CMPR, UOW
// 

#include "RunAction.hh"
#include "AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4RunManager.hh"

RunAction::RunAction(AnalysisManager* analysis)
: G4UserRunAction(),
  iHits(0)
{ 
	analysisMan = analysis;

	accumulableManager = G4AccumulableManager::Instance();	//see B1
	accumulableManager -> RegisterAccumulable(iHits);
	
	// Minimum value for acceptable statistics -- default value
	hitsRequired = 10000;
	
	messenger = new RunMessenger(this);
}

RunAction::~RunAction()
{
}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int run_number = aRun->GetRunID();
  G4cout << "### Run " << run_number << " start." << G4endl;

  accumulableManager -> Reset();

  // Create ROOT file, histograms and ntuple
  analysisMan -> book();
 
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "Number of events = " << aRun->GetNumberOfEvent() << G4endl;

// Close the output ROOT file with the results
   analysisMan -> finish(); 

}

void RunAction::IncreaseHitCount()
{
	iHits += 1;
	
	// Make the accumulable thread-safe, in MT
	// Might be expensive, as it has to be called after
	// every hit in the SV
	//accumulableManager->Merge();
	
	if( iHits.GetValue() >= hitsRequired )
	{
		G4RunManager::GetRunManager() -> AbortRun(true);
		
		G4cout << iHits.GetValue() << " reached, stopping run..." << G4endl;
	}
}


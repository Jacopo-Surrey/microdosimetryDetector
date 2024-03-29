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
// Authors: Susanna Guatelli and Francesco Romano
// susanna@uow.edu.au, francesco.romano@ct.infn.it

// Modified by Jacopo Magini: j.magini@surrey.ac.uk

#ifndef DetectorConstruction_H 
#define DetectorConstruction_H 1
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "AnalysisManager.hh"
#include "DetectorMessenger.hh"

class G4VPhysicalVolume;
class DetectorMessenger;
class G4LogicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction(AnalysisManager* analysis, DetectorMessenger* messenger);
    ~DetectorConstruction();

	// virtual methods to override
    G4VPhysicalVolume* Construct();
    void ConstructSDandField();
	
	G4bool IsDetectorInsidePhantom() { return usingPhantom; }
	
private:
	AnalysisManager* analysis;
	DetectorMessenger* messenger;
	
	G4String detectorType;
	
	G4double detectorPositionDepth;
	
	G4double detectorSizeWidth;
	G4double detectorSizeThickness;
	G4double secondStageSizeDim;
	G4double secondStageThickness;
	
	G4int pixelKinScoring;
	
	G4bool usingPhantom;
	
	G4bool multiSV;
	G4double multiSVbreadth;
	G4double SVspacing;
	G4int nOfSV; // per row
	
	G4double phantomSize; // width
	G4double highPRegionBufferSize; // size of the water region around the detector edge
		// where the precision (cut) is still high
	G4double requiredWidth;
	
	//volumes to keep track of -- WRITE BETTER
	G4VPhysicalVolume* physical_world;
	G4LogicalVolume* logical_motherVolumeForDetector;
	
	// Methods called by Construct() depending on the chosen setup
	void ConstructWithWaterPhantom();
	void ConstructWithoutWaterPhantom(); // causes segfault when applying physics
	
	void ConstructDiamondDetector();	// commented out
	void ConstructMicroDiamondDetector();
	void ConstructWPMicroDiamondDetector();	//review later -- highP region size to be changed accordingly -- move SV to middle of highP
	void ConstructTelescopeDetector();	// to do: move SV to middle of highP
	void ConstructSiliconDetector();	// commented out
	void ConstructWaterPixelDetector();
	
	G4bool checkOverlap = false;
};
#endif

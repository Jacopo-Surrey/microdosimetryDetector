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
//
// Created by Jacopo Magini: j.magini@surrey.ac.uk

#include "DetectorMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

#include "G4UIcmdWithoutParameter.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"

DetectorMessenger::DetectorMessenger(AnalysisManager* analysis_manager)
{
    changeTheGeometryDir = new G4UIdirectory("/geometrySetup/");
    changeTheGeometryDir -> SetGuidance("Geometry setup");
    
    changeTheDetectorCmd = new G4UIcmdWithAString("/geometrySetup/selectDetector", this);
    changeTheDetectorCmd -> SetGuidance("Select the type of detector you wish to use");
    changeTheDetectorCmd -> SetParameterName("Material", false);
	changeTheDetectorCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);
	
	changeDetectorPositionDir = new G4UIdirectory("/geometrySetup/detectorPosition/");
	changeDetectorPositionDir -> SetGuidance("Modify the placement of the detector");
	
	changeDetectorPositionDepthCmd = new G4UIcmdWithADoubleAndUnit("/geometrySetup/detectorPosition/setDepth", this);
	changeDetectorPositionDepthCmd -> SetGuidance("Set the detector depth inside the water phantom");
	changeDetectorPositionDepthCmd -> SetParameterName("Depth", false);
	changeDetectorPositionDepthCmd -> SetRange("Depth >= 10. && Depth <= 290.");	//the phantom is 300.*mm deep
	changeDetectorPositionDepthCmd -> SetUnitCategory("Length");
	changeDetectorPositionDepthCmd -> SetDefaultUnit("mm");
	changeDetectorPositionDepthCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);
	
	changeDetectorDimensionDir = new G4UIdirectory("/geometrySetup/detectorDimension/");
	changeDetectorDimensionDir -> SetGuidance("Modify the dimensions of the detector");
	
	changeDetectorSizeWidthCmd = new G4UIcmdWithADoubleAndUnit("/geometrySetup/detectorDimension/setWidth", this);
	changeDetectorSizeWidthCmd -> SetGuidance("Set the width of the detector");
	changeDetectorSizeWidthCmd -> SetParameterName("Width", false);
	changeDetectorSizeWidthCmd -> SetRange("Width >= 0.1 && Width <= 1000.");	// DOES IT MAKE SENSE???
	changeDetectorSizeWidthCmd -> SetUnitCategory("Length");
	changeDetectorSizeWidthCmd -> SetDefaultUnit("um");
	changeDetectorSizeWidthCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);
	
	changeDetectorSizeThicknessCmd = new G4UIcmdWithADoubleAndUnit("/geometrySetup/detectorDimension/setThickness", this);
	changeDetectorSizeThicknessCmd -> SetGuidance("Set the thickness of the detector");
	changeDetectorSizeThicknessCmd -> SetParameterName("Thickness", false);
	changeDetectorSizeThicknessCmd -> SetRange("Thickness >= 0.1 && Thickness <= 1000.");	// DOES IT MAKE SENSE???
	changeDetectorSizeThicknessCmd -> SetUnitCategory("Length");
	changeDetectorSizeThicknessCmd -> SetDefaultUnit("um");
	changeDetectorSizeThicknessCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);
	
	enableWaterPhantomCmd = new G4UIcmdWithABool("/geometrySetup/enableWaterPhantom", this);
	enableWaterPhantomCmd -> SetGuidance("If true, the detector is placed inside a water phantom");
	enableWaterPhantomCmd -> SetParameterName("Phantom", false);
	enableWaterPhantomCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);
	
	useMultipleSVCmd = new G4UIcmdWithABool("/geometrySetup/useMultipleSV", this);
	useMultipleSVCmd -> SetGuidance("Set to true for multiple SV, false for a single one");
	useMultipleSVCmd -> SetParameterName("MultipleSV", true);
	useMultipleSVCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);
	
	applyChangesToGeometryCmd = new G4UIcmdWithoutParameter("/geometrySetup/applyChanges",this);
    applyChangesToGeometryCmd -> SetGuidance("Apply selected changes to the geometry");
	applyChangesToGeometryCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);
	
	analysis = analysis_manager;
	
	// default values
	detectorType = "MicroDiamond";
	detectorDepth = 10*mm;
	detectorWidth = 100.*um;
	detectorThickness = 8.*um;
	usingPhantom = true;	// FIX ME: currently the program crashes at various stages if this is set to false
	multiSV = false;
}

DetectorMessenger::~DetectorMessenger()
{
    delete changeTheDetectorCmd;
	delete changeDetectorPositionDepthCmd;
	delete changeDetectorSizeWidthCmd;
	delete changeDetectorSizeThicknessCmd;
	delete enableWaterPhantomCmd;
	delete useMultipleSVCmd;
	
	delete applyChangesToGeometryCmd;
	
	delete changeDetectorPositionDir;
	delete changeDetectorDimensionDir;
	delete changeTheGeometryDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String commandContent)
{

	if( command == changeTheDetectorCmd )
	{
		if( commandContent == "Diamond" || commandContent == "MicroDiamond" || commandContent == "Silicon" )
		{
			detectorType = commandContent;
			
			G4cout << "Detector type changed to " << commandContent << G4endl;
			G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
			
			/*
			DetectorConstruction* newDetector = new DetectorConstruction(analysis, commandContent);	//PUT THESE LINES IN NEW COMMAND
			G4RunManager *runManager = G4RunManager::GetRunManager();
			runManager -> SetUserInitialization(newDetector);
			runManager -> GeometryHasBeenModified();
			*/
		}
		
		else
		{
			G4cout <<"Unknown detector type: " << commandContent << ". Geometry not changed" << G4endl;
		}
	}
	
	else if( command == changeDetectorPositionDepthCmd )
	{
		detectorDepth = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
		
		G4cout << "Detector depth in water changed to " << commandContent << G4endl;
		if( usingPhantom == false )	G4cout << "However the water phantom is not enabled. Enable it with '/geometrySetup/enableWaterPhantom true' or this value will be ignored" << G4endl;
		G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
	}
	
	else if( command == changeDetectorSizeWidthCmd )
	{
		detectorWidth = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
		
		G4cout << "Detector width changed to " << commandContent << G4endl;
		G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
	}
	
	else if( command == changeDetectorSizeThicknessCmd )
	{
		detectorThickness = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(commandContent);
		
		G4cout << "Detector thickness changed to " << commandContent << G4endl;
		G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
	}
	
	else if( command == enableWaterPhantomCmd )
	{
		usingPhantom = G4UIcmdWithABool::GetNewBoolValue(commandContent);
		
		if( usingPhantom == true )	G4cout << "Water phantom enabled" << G4endl;
		else if( usingPhantom == false )	G4cout << "Water phantom disabled" << G4endl;
		else	G4cout << "Error: " << commandContent << "is not a bool" << G4endl;
		
		G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
	}
	
	else if( command == useMultipleSVCmd )
	{
		multiSV = G4UIcmdWithABool::GetNewBoolValue(commandContent);
		
		if( multiSV == true )	G4cout << "Now using multiple SV" << G4endl;
		else if( multiSV == false )	G4cout << "Now using a single SV" << G4endl;
		else	G4cout << "Error: " << commandContent << "is not a bool" << G4endl;
		
		G4cout << "Run /geometrySetup/applyChanges to apply" << G4endl;
	}
	
	else if( command == applyChangesToGeometryCmd )
	{
		G4RunManager* runManager = G4RunManager::GetRunManager();
		
		// if geometryHasChanged == true
		DetectorConstruction* newDetector = new DetectorConstruction(analysis, this);
		runManager -> SetUserInitialization(newDetector);
		runManager -> GeometryHasBeenModified();
		
		// (nested?) if regionsHaveChanged == true
		G4VUserPhysicsList* newPhysics = new PhysicsList(this);
		runManager -> SetUserInitialization(newPhysics);
		runManager -> PhysicsHasBeenModified();
		
		G4cout << "Changes to geometry have been applied" << G4endl;		
	}
}

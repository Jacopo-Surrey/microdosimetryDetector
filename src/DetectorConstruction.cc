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

#include "DetectorConstruction.hh"
#include "globals.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
//#include "G4SubtractionSolid.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "SensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"

// set the active SV to output
// being static, they can be accessed from the SteppinAction
#ifndef USING_SILICON 
	G4int DetectorConstruction::sensitiveVolumeToOutput = 1;	// 0=50x50, 1=300x300, 2=100x100, 3=200x200
	G4double DetectorConstruction::detector_thickness = 8.*micrometer; // set detector thickness
	G4double DetectorConstruction::dd = 24.*millimeter; // detecotr-source distance. NB change GPS source settings as well!!
#else
	G4int DetectorConstruction::sensitiveVolumeToOutput = 12;	// the one in the very middle 
#endif

DetectorConstruction::DetectorConstruction(AnalysisManager* analysis_manager)
{
analysis = analysis_manager;
}

DetectorConstruction::~DetectorConstruction(){

}

#ifndef USING_SILICON
//#ifndef USING_SILICON	// by default, build the diamond detector
G4VPhysicalVolume* DetectorConstruction::Construct()
{

	//Define each individual element
	//Define Nitrogen
	G4double A = 14.01 * g/mole;
	G4double Z = 7;
	G4Element* elN = new G4Element ("Nitrogen", "N", Z, A);

	//Define Oxygen
	A = 16.0 * g/mole;
	Z = 8;
	G4Element* elO = new G4Element ("Oxygen", "O", Z, A);

	//Define Hydrogen 
	A = 1.01 * g/mole;
	Z = 1;
	G4Element* elH = new G4Element ("Hydrogen", "H", Z, A);
	
	//Define Boron
	A = 10.8 * g/mole;
	Z = 5;
	G4Element* elB = new G4Element ("Boron", "B", Z, A);

	//Define Carbon
	A = 12.01 * g/mole;
	Z = 6;
	G4Element* elC = new G4Element ("Carbon", "C", Z, A);

	//Define Air   
	G4Material* Air = new G4Material("Air", 1.29*mg/cm3, 2);
	Air -> AddElement(elN, 70*perCent);
	Air -> AddElement(elO, 30*perCent);

	//Define diamond
	A = 12.01 * g/mole;
	Z = 6;
	G4Material* diamond = new G4Material("diamond", Z, A, 3.515*g/cm3);
			
	//Define p-type diamond (boron doped diamond)
	G4Material* p_diamond = new G4Material("p_diamond", 3.514*g/cm3, 2);
	// Boron concentration used is 1e20 cm-3, considering the diamond density and a Boron atomic weight of 10.811u
	p_diamond -> AddElement(elC, 99.94887*perCent);
	p_diamond -> AddElement(elB, 0.05113*perCent);

	//Define chromium contact
	G4Material* chromium = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cr");

	//Define Aluminium contacts (AlContact)
	//A = 26.981 * g/mole;
	//Z = 13;
	//G4Material* AlContact = new G4Material("AlContact", Z, A, 2.7 *g/cm3);

	//Define Gold contact (AuContact)
	//A = 196.97 * g/mole;
	//Z = 79;
	//G4Material* AuContact = new G4Material("AuContact", Z, A, 19.3 *g/cm3);

	//define water
	G4Material* water = new G4Material("water", 1*g/cm3, 2);
	water -> AddElement(elH, 2);
	water -> AddElement(elO, 1);
	
	//Define Vacuum
	G4double vacuumDensity = 1.e-25 *g/cm3;
	G4double pressure = 3.e-18*pascal;
	G4double temperature = 2.73*kelvin;
	G4Material* vacuum = new G4Material("Galactic", Z=1., A=1.01*g/mole,
			         vacuumDensity,kStateGas,temperature,pressure);

	//Define volumes
	// World volume
	// if box
	G4double worldx = 6.*mm /2.;  //half length!!!!
	G4double worldy = 6.*mm /2.;
	// if cylinder
	G4double worldrz = 1.5*mm;  
	G4double worldrZ = 3.*mm;
	// the height is the same for both shapes
	G4double worldz = (dd+2.) /2.; //half length!!!!


	// World volume, containing all geometry
	G4Box* world = new G4Box("world_vol", worldx, worldy, worldz);
        //G4Cons* world = new G4Cons("world_vol", 0.*mm, worldrz, 0.*mm, worldrZ, worldz, 0*deg, 360*deg);

	G4LogicalVolume* logical_world = new G4LogicalVolume(world, vacuum, "world_log", 0,0,0);

	//set the logical world volume invisible
	logical_world -> SetVisAttributes(G4VisAttributes::GetInvisible());

	G4VPhysicalVolume* physical_world = new G4PVPlacement(0,
								G4ThreeVector(),
								logical_world, 
								"world_phys",
								0, 
								false, 
								0);


/*	// Whole box geometry
	// I ignore the Aluminium walls. The inside of the chamber is filled with air
	G4double chamVol_x = 250*millimeter /2.;
	G4double chamVol_y = 150*millimeter /2.;
	G4double chamVol_z = 116*millimeter /2.;

	G4Box* chamVol_box = new G4Box("chamVol_box", chamVol_x, chamVol_y, chamVol_z);
 
	G4LogicalVolume* logical_chamVol = new G4LogicalVolume(chamVol_box, Air, "chamVol_log",0,0,0);
 
	new G4PVPlacement(0, G4ThreeVector(0,0,0), logical_chamVol,"chamVol_phys",
					logical_world, 
					false, 0, true);
 
	//logical_chamVol -> SetVisAttributes(G4VisAttributes(G4Colour(255,255,255))); //white
	logical_chamVol -> SetVisAttributes(G4VisAttributes::GetInvisible());
 
	// 4 later: smaller air volume where I can lower the cuts with G4Region
	G4double highPVol_x = 8*millimeter /2.; 
	G4double highPVol_y = 8*millimeter /2.;
	G4double highPVol_z = 6*millimeter /2.;

	G4Box* highPVol_box = new G4Box("highPVol_box", highPVol_x, highPVol_y, highPVol_z);
 
	G4LogicalVolume* logical_highPVol = new G4LogicalVolume(highPVol_box, Air, "highPVol_log",0,0,0);
 
	G4double sourceDetectorAxis_x = chamVol_x - 70*millimeter;
 
	G4ThreeVector highP_position = G4ThreeVector( sourceDetectorAxis_x , 0 , - chamVol_z + 20*millimeter );
 
	new G4PVPlacement(0, highP_position, logical_highPVol,"highPVol_phys",
				logical_chamVol, 
				false, 0, true);
 
	logical_highPVol -> SetVisAttributes(G4VisAttributes(G4Colour(255,255,255))); //white
 
	// Sensitive volumes
	G4double SVspacing_ee = 200.*micrometer; //edge-edge distance

	G4double SVthickness = detector_thickness /2.; // hald the detector thickness defined at the beginning of this document
 
	// 4 sensitive volumes
	G4double SVside[4] = { 50.*micrometer /2., 300.*micrometer /2., 100.*micrometer /2., 200.*micrometer /2.}; // half side!!
 
	G4Box* SV_box[4];
	G4LogicalVolume* logical_SV[4];
 
	std::ostringstream name;
 
	//prepare its colour
	G4VisAttributes SVcolour(G4Colour(0.5, 0.5, 0.5));
	SVcolour.SetForceSolid(true);
 
	G4ThreeVector SVposition ={ -(1.5*SVspacing_ee + 2*SVside[2] +  2*SVside[1]), 0, SVthickness};	//position of the first edge (left) 
	for( int i=0; i<4; i++)
	{
		name << "SV_box_" << i;
		SV_box[i] = new G4Box(name.str(), SVside[i], SVside[i], SVthickness);
		name.str(""); //clears the string 
		
		name << "SV_log_" << i;
		logical_SV[i] = new G4LogicalVolume(SV_box[i], diamond, name.str(), 0,0,0);
		name.str("");
		
		name << "SV_phys_" << i;
		SVposition[0] += SVside[i];	// update the position to the centre 
		new G4PVPlacement(0, SVposition, logical_SV[i], name.str(),
					logical_highPVol,
					false, 0, true);
		name.str("");
		SVposition[0] += SVside[i] + SVspacing_ee; // update x-position to the edge of the next SV
		
		logical_SV[i] -> SetVisAttributes(SVcolour);
	}
*/
	// simplified geometry
	// distance between source and microdosimeter
	G4double chamVol_rz = 1.*millimeter;
	G4double chamVol_rZ = 2.6*millimeter;
	G4double chamVol_z = (dd+1.) /2.; //half length!!!!

        G4Cons* chamVol_cons = new G4Cons("chamVol_cons", 0.*mm, chamVol_rz, 0.*mm, chamVol_rZ, chamVol_z, 0*deg, 360*deg);
 
	G4LogicalVolume* logical_chamVol = new G4LogicalVolume(chamVol_cons, Air, "chamVol_log",0,0,0);
	//G4LogicalVolume* logical_chamVol = new G4LogicalVolume(chamVol_cons, vacuum, "chamVol_log",0,0,0);
 
	new G4PVPlacement(0, G4ThreeVector(0,0,0), logical_chamVol,"chamVol_phys",
					logical_world, 
					false, 0, true);
 
	logical_chamVol -> SetVisAttributes(G4VisAttributes(G4Colour(255,255,255))); //white
	//logical_chamVol -> SetVisAttributes(G4VisAttributes::GetInvisible());
 
	// 4 later: smaller air volume where I can lower the cuts with G4Region
	G4double highPVol_x = 1.5*millimeter /2.; 
	G4double highPVol_y = 0.7*millimeter /2.;
	G4double highPVol_z = 0.7*millimeter /2.; // has to be at least 0.31 as it contains p-type and substrate on the back-half (301 um)

	G4Box* highPVol_box = new G4Box("highPVol_box", highPVol_x, highPVol_y, highPVol_z);
 
	G4LogicalVolume* logical_highPVol = new G4LogicalVolume(highPVol_box, Air, "highPVol_log",0,0,0);
	//G4LogicalVolume* logical_highPVol = new G4LogicalVolume(highPVol_box, vacuum, "highPVol_log",0,0,0);
 
	G4double detectorPosition_z = -dd/2.;
 
	G4ThreeVector highP_position = G4ThreeVector( 0. , 0 , detectorPosition_z );
 
	new G4PVPlacement(0, highP_position, logical_highPVol,"highPVol_phys",
				logical_chamVol, 
				false, 0, true);
 
	logical_highPVol -> SetVisAttributes(G4VisAttributes(G4Colour(255,255,255))); //white
 
	// Sensitive volumes
	G4double SVspacing_ee = 200.*micrometer; //edge-edge distance

	G4double SVthickness = detector_thickness /2.; // half the detector thickness
 
	// 4 sensitive volumes
	G4double SVside[4] = { 50.*micrometer /2., 300.*micrometer /2., 100.*micrometer /2., 200.*micrometer /2.}; // half side!!
 
	G4Box* SV_box[4];
	G4LogicalVolume* logical_SV[4];
 
	std::ostringstream name;
 
	//prepare colours
	G4VisAttributes SVcolour(G4Colour(0.5, 0.5, 0.5));
	SVcolour.SetForceSolid(true);
	
	G4VisAttributes fe_colour(G4Colour::Brown());
	fe_colour.SetForceSolid(true);
	G4VisAttributes pDcolour(G4Colour::Blue());
	pDcolour.SetForceSolid(true);

	// chromium front-electrode
	G4Box* fe_box[4];
	G4LogicalVolume* logical_fe[4];
	G4double fet = 50.*nanometer /2.; // front-electrode half-thickness

	// p-type diamond
	G4Box* pD_box[4];
	G4LogicalVolume* logical_pD[4];
	G4double pDt = 1.*micrometer /2.; // p-type diamond back-electrode half-thickness
 
	G4ThreeVector SVposition ={ -(1.5*SVspacing_ee + 2*SVside[2] +  2*SVside[1]), 0, SVthickness};	//position of the first edge (left)
	G4ThreeVector fe_position = SVposition;
	fe_position[2] = 2*SVthickness + fet; // updating z-coordinate considering they are halves thicknesses
	G4ThreeVector pDposition = SVposition;
	pDposition[2] = -pDt;

	for( int i=0; i<4; i++)
	{
		name << "SV_box_" << i;
		SV_box[i] = new G4Box(name.str(), SVside[i], SVside[i], SVthickness);
		name.str(""); //clears the string 
		
		name << "SV_log_" << i;
		logical_SV[i] = new G4LogicalVolume(SV_box[i], diamond, name.str(), 0,0,0);
		name.str("");
		
		name << "SV_phys_" << i;
		SVposition[0] += SVside[i];	// update the position to the centre 
		new G4PVPlacement(0, SVposition, logical_SV[i], name.str(),
					logical_highPVol,
					false, 0, true);
		name.str("");
		
		logical_SV[i] -> SetVisAttributes(SVcolour);

		// chromium front-electrode
		name << "fe_box_"<< i;
		fe_box[i] = new G4Box(name.str(), SVside[i], SVside[i], fet);
		name.str("");
		name << "fe_log_"<< i;
		logical_fe[i]= new G4LogicalVolume(fe_box[i], chromium, name.str(), 0,0,0);
		name.str("");
		name << "fe_phys_" << i;
		fe_position[0]=SVposition[0];
		new G4PVPlacement(0, fe_position, logical_fe[i], name.str(),
					logical_highPVol,
					false, 0, true);
		name.str("");

		logical_fe[i] -> SetVisAttributes(fe_colour);

		// p-type diamond back-electrode
		name << "pD_box_"<< i;
		pD_box[i] = new G4Box(name.str(), SVside[i], SVside[i], pDt);
		name.str("");
		name << "pD_log_"<< i;
		logical_pD[i]= new G4LogicalVolume(pD_box[i], p_diamond, name.str(), 0,0,0);
		name.str("");
		name << "pD_phys_" << i;
		pDposition[0]=SVposition[0];
		new G4PVPlacement(0, pDposition, logical_pD[i], name.str(),
					logical_highPVol,
					false, 0, true);
		name.str("");

		logical_pD[i] -> SetVisAttributes(pDcolour);


		SVposition[0] += SVside[i] + SVspacing_ee; // update x-position to the edge of the next SV
		
	}



	// HPHT diamondsubstrate, chosen to simulate only one big substrate for simplicity as it does not affect the simulation results. Actually I should not be a common substrate and its dimension (but the thickness that is correct) are differt
	G4double subsx = 1.5*millimeter /2.; // substrate half-side x (taken as the highPvol)
	G4double subsy = 0.7*millimeter /2.; // substrate half-side y (taken as the highPvol)
	G4double subt = 300.*micrometer /2.; // substrate half-thicknesss
	G4Box* sub_box = new G4Box("sub_box", subsx, subsy, subt);
	G4LogicalVolume* logical_sub = new G4LogicalVolume(sub_box, diamond, "sub_log", 0,0,0);
	new G4PVPlacement(0, {0,0, -(2*pDt+subt)}, logical_sub, "sub_phys",
				logical_highPVol,
				false, 0, true);

	logical_sub -> SetVisAttributes(SVcolour);
	
	// high precision region
	G4Region* highPRegion = new G4Region("highPRegion");
	highPRegion->AddRootLogicalVolume(logical_highPVol);
        
	return physical_world; 
}

#else	// if the flag is on, build the reference silicon instead
G4VPhysicalVolume* DetectorConstruction::Construct()
{
	//load NIST database
	G4NistManager* nist = G4NistManager::Instance();
	nist->SetVerbose(1);

	//MEMO: remove individual materials
	//Define each individual element
	//Define Nitrogen
	G4double A = 14.01 * g/mole;
	G4double Z = 7;
	G4Element* elN = new G4Element ("Nitrogen", "N", Z, A);

	//Define Oxygen
	A = 16.0 * g/mole;
	Z = 8;
	G4Element* elO = new G4Element ("Oxygen", "O", Z, A);

	//Define Hydrogen 
	A = 1.01 * g/mole;
	Z = 1;
	G4Element* elH = new G4Element ("Hydrogen", "H", Z, A);

	//Define Boron
	A = 10.8 * g/mole;
	Z = 5;
	G4Element* elB = new G4Element ("Boron", "B", Z, A);

	//Define Carbon
	A = 12.01 * g/mole;
	Z = 6;
	G4Element* elC = new G4Element ("Carbon", "C", Z, A);

	//Define Silicon
	A = 28.085 * g/mole;
	Z = 14;
	G4Element* elSi = new G4Element ("Silicon", "Si", Z, A);
	
	//Define Air   
	G4Material* Air = new G4Material("Air", 1.29*mg/cm3, 2);
	Air -> AddElement(elN, 70*perCent);
	Air -> AddElement(elO, 30*perCent);
 
	//Define PMMA (C502H8)
	// NIST reference 
	G4Material* PMMA = new G4Material("PMMA", 1.19*g/cm3, 3);
	PMMA -> AddElement(elC, 5);
	PMMA -> AddElement(elO, 2);
	PMMA -> AddElement(elH, 8);

	//define water
	G4Material* water = new G4Material("water", 1*g/cm3, 2);
	water -> AddElement(elH, 2);
	water -> AddElement(elO, 1);

	//define materials
	G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");
	G4Material* SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	
	
	//Define volumes
	// World volume
	G4double worldx = 1*m /2.;  //half length!!!!
	G4double worldy = 0.5*m /2.;
	G4double worldz = 0.5*m /2.;
	
	// World volume, containing all geometry
	G4Box* world = new G4Box("world_box", worldx, worldy, worldz);

	G4LogicalVolume* logical_world = new G4LogicalVolume(world, Air, "world_log", 0,0,0);

	//set the logical world volume invisible
	logical_world -> SetVisAttributes(G4VisAttributes::GetInvisible());

	G4VPhysicalVolume* physical_world = new G4PVPlacement(0,
								G4ThreeVector(),
								logical_world, 
								"world_phys",
								0, 
								false, 
								0);
	
	//water phantom
	G4double phantom_x = 300.*mm /2.;
	G4double phantom_y = 300.*mm /2.;
	G4double phantom_z = 300.*mm /2.;
	
	G4Box* phantom_box = new G4Box("phantom_box", phantom_x, phantom_y, phantom_z);
	
	G4LogicalVolume* logical_phantom = new G4LogicalVolume(phantom_box, water, "phantom_log", 0,0,0);
	
	G4ThreeVector phantom_position = G4ThreeVector( phantom_x, 0, 0 );	//the phantom starts at x=0
	
	 new G4PVPlacement(0, phantom_position, logical_phantom,"phantom_phys",
				logical_world, 
				false, 0, true);
	
	logical_phantom -> SetVisAttributes(G4VisAttributes(G4Colour(0., 0.2, 0.6)));
	
	//smaller volume where I can lower the cuts with G4Region
	G4double highPVol_x = 6.*mm /2.; 
	G4double highPVol_y = 11.*mm /2.;
	G4double highPVol_z = 11.*mm /2.;

	G4Box* highPVol_box = new G4Box("highPVol_box", highPVol_x, highPVol_y, highPVol_z);
 
	G4LogicalVolume* logical_highPVol = new G4LogicalVolume(highPVol_box, water, "highPVol_log",0,0,0);
	
	G4double detectorDepth = 10.*mm;	// CHANGE ME
	
	G4double detectorCentre_x = -phantom_x + detectorDepth + highPVol_x;
	
	G4ThreeVector highP_position = G4ThreeVector( detectorCentre_x , 0 , 0 );
 
	new G4PVPlacement(0, highP_position, logical_highPVol,"highPVol_phys",
				logical_phantom, 
				false, 0, true);
 
	logical_highPVol -> SetVisAttributes(G4VisAttributes(G4Colour(0., 0., 1.)));
	
	// I need to set the size of the SV now, because some other parameters depend on it
	G4double SV_thick = 25.*um /2.; 
	
	// PMMA
	G4double PMMA_x = SV_thick;
	G4double PMMA_y = 5.*mm /2.;
	G4double PMMA_z = 5.*mm /2.;
	
	G4Box* PMMA_box = new G4Box("PMMA_box", PMMA_x, PMMA_y, PMMA_z);
	
	G4LogicalVolume* logical_PMMA = new G4LogicalVolume(PMMA_box, PMMA, "PMMA_log", 0,0,0);
	
	new G4PVPlacement(0, G4ThreeVector(), logical_PMMA, "PMMA_phys",
					logical_highPVol,
					false, 0, true);
	
	logical_PMMA -> SetVisAttributes(G4VisAttributes(G4Colour(0., 1., 0.)));
	
	// sensitive volumes 
	G4double SV_radius = SV_thick *2;	// this is the full length, not the half lenght!
	G4double SVspacing_ee = 10.*um;	// distance between edges of two SV
	
	G4Tubs* SV_cyl = new G4Tubs("SV_cyl", 0., SV_radius, SV_thick, 0., 2*M_PI*rad);
	
	G4LogicalVolume* logical_SV[25];	// 5x5 volumes (arbitrary number)
	
	G4RotationMatrix* cylRot = new G4RotationMatrix;
	cylRot->rotateY(M_PI/2.*rad);
	
	G4ThreeVector SVposition;
 
	std::ostringstream name;
	
	// prepare its colour
	G4VisAttributes SVcolour(G4Colour(0.5, 0.5, 0.5));
	SVcolour.SetForceSolid(true);
	
	G4int index;
	for(int i=-2; i<3; i++)
	{
		for(int j=-2; j<3; j++)
		{
			index = i+2 + (j+2)*5;
			name << "SV_log_" << index;
			logical_SV[index] = new G4LogicalVolume(SV_cyl, silicon, name.str(), 0,0,0);
			name.str("");
		
			name << "SV_phys_" << index;
			SVposition = { 0., i*SVspacing_ee + i*2.*SV_radius, j*SVspacing_ee + j*2.*SV_radius };
			new G4PVPlacement(cylRot, SVposition, logical_SV[index], name.str(),
						logical_PMMA,
						false, 0, true);
			name.str("");
		
			logical_SV[index] -> SetVisAttributes(SVcolour);
		}
	}
	
	// Si02 layer
	G4double oxyde_x = 1.*um /2.;
	G4double oxyde_y = PMMA_y;
	G4double oxyde_z = PMMA_z;
	
	G4Box* oxyde_box = new G4Box("oxyde_box", oxyde_x, oxyde_y, oxyde_z);
	
	G4LogicalVolume* logical_oxyde = new G4LogicalVolume(oxyde_box, SiO2, "oxyde_log", 0,0,0);
	
	G4ThreeVector oxyde_position = G4ThreeVector( PMMA_x + oxyde_x, 0, 0 );
	new G4PVPlacement(0, oxyde_position, logical_oxyde, "oxyde_phys",
					logical_highPVol,
					false, 0, true);
	
	logical_oxyde -> SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 0.6)));
	
	// high precision region
	G4Region* highPRegion = new G4Region("highPRegion");
	highPRegion->AddRootLogicalVolume(logical_highPVol);
	
	return physical_world;
}
#endif

void DetectorConstruction::ConstructSDandField()
{
	std::ostringstream outName; outName << "SV_log_" << sensitiveVolumeToOutput;

	SensitiveDetector* SD = new SensitiveDetector("SD", "DetectorHitsCollection", analysis);
	G4SDManager::GetSDMpointer()->AddNewDetector(SD);
	SetSensitiveDetector(outName.str(), SD);	
}

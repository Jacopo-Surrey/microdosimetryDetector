#include "PhysicsList.hh"

#include <G4EmLivermorePhysics.hh>
#include <G4DecayPhysics.hh> 
#include <G4ProductionCutsTable.hh>
#include <G4SystemOfUnits.hh>

#include <G4HadronElasticPhysicsHP.hh>
#include <G4HadronPhysicsQGSP_BIC_HP.hh>
#include <G4IonBinaryCascadePhysics.hh>

#include <G4RegionStore.hh>

#include "G4StepLimiterPhysics.hh"

PhysicsList::PhysicsList(DetectorMessenger* detectorMessenger)
{
	//Low energy EM physics 
	RegisterPhysics(new G4EmLivermorePhysics());
	
	//Decay physics
	RegisterPhysics(new G4DecayPhysics());
	
	//Hadronic physics
	RegisterPhysics(new G4HadronElasticPhysicsHP());
	RegisterPhysics(new G4HadronPhysicsQGSP_BIC_HP());
	RegisterPhysics(new G4IonBinaryCascadePhysics());

	//Add step limiter (as defined in DetectorConstruction)
	RegisterPhysics(new G4StepLimiterPhysics());
	
	// check if I'm using a phantom, and therefore if I need
	// a higher precision in the inner region
	usingRegion = detectorMessenger->GetUsingPhantomBool();
}


void PhysicsList::SetCuts()
{
	if( usingRegion == true )	defaultCutValue = 100.*micrometer;
	else if( usingRegion == false )	defaultCutValue = 0.1*micrometer;
	
	// default cut values
	fCutForGamma     = defaultCutValue;
	fCutForElectron  = defaultCutValue;
	fCutForPositron  = defaultCutValue;
	fCutForProton    = defaultCutValue;
	
	// set cuts everywhere
	SetCutValue(fCutForGamma, "gamma");
	SetCutValue(fCutForElectron, "e-");
	SetCutValue(fCutForPositron, "e+");
	SetCutValue(fCutForProton, "proton");
	
	if( usingRegion == true )
	{
		// set different cuts in the high precision region
		G4String regName = "highPRegion";
		G4Region* region = G4RegionStore::GetInstance()->GetRegion(regName);
		
		G4ProductionCuts* cuts = new G4ProductionCuts;
		G4double highPCut = 0.1*um;
		
		cuts->SetProductionCut(highPCut, G4ProductionCuts::GetIndex("gamma"));
		cuts->SetProductionCut(highPCut, G4ProductionCuts::GetIndex("e-"));
		cuts->SetProductionCut(highPCut, G4ProductionCuts::GetIndex("e+"));
		cuts->SetProductionCut(highPCut, G4ProductionCuts::GetIndex("proton"));
		
		region->SetProductionCuts(cuts);
	}
	
	// Update the production cuts table energy range
	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250.*eV,10.*GeV);  
}

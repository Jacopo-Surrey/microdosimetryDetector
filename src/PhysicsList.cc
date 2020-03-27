#include "PhysicsList.hh"

#include <G4EmLivermorePhysics.hh>
#include <G4DecayPhysics.hh> 
#include <G4ProductionCutsTable.hh>
#include <G4SystemOfUnits.hh>

#include <G4HadronElasticPhysicsHP.hh>
#include <G4HadronPhysicsQGSP_BIC_HP.hh>
#include <G4IonBinaryCascadePhysics.hh>

PhysicsList::PhysicsList()
{
  //Low energy EM physics 
  RegisterPhysics(new G4EmLivermorePhysics());
  
  //Decay physics
  RegisterPhysics(new G4DecayPhysics());

  //Hadronic physics
  RegisterPhysics(new G4HadronElasticPhysicsHP());
  RegisterPhysics(new G4HadronPhysicsQGSP_BIC_HP());
  RegisterPhysics(new G4IonBinaryCascadePhysics());
}


void PhysicsList::SetCuts()
{
  // The method SetCuts() is mandatory in the interface. Here, one just use 
  // the default SetCuts() provided by the base class.
  G4VUserPhysicsList::SetCuts();
  
  // Update the production cuts table energy range
  // G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV,100.*GeV);  
}

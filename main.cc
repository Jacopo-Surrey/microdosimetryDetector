//Based on advanced/Radioprotection

#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4MTRunManager.hh"
#include "AnalysisManager.hh"
#include "ActionInitialization.hh"
#include "AnalysisManager.hh"
#include "G4UIExecutive.hh"


int main(int argc, char** argv)
{

#ifdef G4MULTITHREADED
  G4MTRunManager* pRunManager = new G4MTRunManager;
  pRunManager->SetNumberOfThreads(4); // Is equal to 2 by default
#else  
  G4RunManager* pRunManager = new G4RunManager;
#endif

  AnalysisManager* analysis = new AnalysisManager();
 
  DetectorConstruction* detector = new DetectorConstruction(analysis);  
    
  pRunManager -> SetUserInitialization(detector);

  G4VUserPhysicsList* physics = new PhysicsList();
  
  pRunManager -> SetUserInitialization(physics); 

   // User action initialization  

  ActionInitialization* actions = new ActionInitialization(analysis);
  pRunManager->SetUserInitialization(actions);

  G4VisManager* visManager = new G4VisExecutive();
  visManager->Initialize();
 
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  if(argc == 1){

    /*this sets up the user interface to run in interactive mode */
    G4UIExecutive* ui = new G4UIExecutive(argc, argv); 
    G4cout << " UI session starts ..." << G4endl;
    UImanager -> ApplyCommand("/control/execute vis.mac");
    //now, we run in interactive mode so tell the UI manager to read the vis.mac macro file and 
    //UI->ApplyCommand("/control/execute vis.mac"); 
    ui -> SessionStart();
    delete ui;

  } else { 
    //otherwise we run in batch mode
    G4String command = "/control/execute ";//create first part of command
    G4String fileName = argv[1];//second part is the file name that was typed at the command line 
    UImanager->ApplyCommand(command+fileName);//join the two and pass to the UI manager for interpretation
  }

  delete visManager;
  delete analysis; 
  delete pRunManager;
  
  return 0;
}

//Simple example to access the qtree in DIANA root files 
//To run do:
//.x example1_QChain.C
//
 

{
  //the tree will always be named qtree 
  QChain aChain("qtree"); 
 
  //QChain uses list files to conveniently chain qtrees from individual files in the list  
  aChain.Add("/nfs/cuore1/data/CUORE0/OfficialProcessed_v02.30/ds2130/Unblinded_202410_C.list"); 
  
  //Use TTree Draw to select some data 
  Int_t nSelected = aChain.Draw("Baseline:BaselineRMS:Time", "Channel == 1", "GOFF") ; //used GOFF to avoid drawing this 3D histogram !! 

  //Another way to access the number of selected rows
  Int_t nSelected2 = aChain.GetSelectedRows();
  
  if(nSelected >0){
    
    Double_t* blArray = aChain.GetV1();//these are arrays of dimension nSelected 
    Double_t* rmsArray = aChain.GetV2();
    Double_t* timeArray = aChain.GetV3();
    //You can do up to GetV4();

    //calculate the rms of the baseline value 
    
    Double_t blMean =0.;
    Double_t blRMS =0.;
    //Get the mean of the baselines
    for(Int_t i=0;i<nSelected;i++){
      blMean += blArray[i];
    }
    blMean/=nSelected;
    
    //Get the RMS of the baselines
    for(Int_t i=0;i<nSelected;i++){
      blRMS += pow(blArray[i]-blMean,2);
    }
    blRMS/=nSelected;
    blRMS = sqrt(blRMS);
    

    //Set up a split canvas
    TCanvas* myCan = new TCanvas("myCan","myCan", 500,800);
    myCan->Divide(1,2);//this splits vertically
    myCan->cd(1);//The top part
    
    //Now let's graph the baseline vs Time and draw a line a mean +/- 1 RMS
    TGraph* blGraph = new TGraph(nSelected,timeArray,blArray);
    blGraph->GetHistogram()->SetXTitle("Time (s)");
    blGraph->GetHistogram()->SetYTitle("Baseline (mV)");
    blGraph->GetHistogram()->SetTitle("Baseline vs Time with +/- 1 RMS band");
    blGraph->Draw("AP");

    TLine* li = new TLine();
    li->SetLineColor(kRed);
    li->DrawLine(timeArray[0],(blMean-blRMS),timeArray[nSelected-1],(blMean-blRMS));//line at blMean-1 RMS
    li->DrawLine(timeArray[0],(blMean+blRMS),timeArray[nSelected-1],(blMean+blRMS));//line at blMean + 1 RMS 
   

    //Next plot BaselineRMS vs Time and compare to calculated RMS   
    myCan->cd(2);//The bottom part
    
    //Now let's graph the baseline vs Time and draw a line a mean +/- 1 RMS
    TGraph* blRMSGraph = new TGraph(nSelected,timeArray,rmsArray);
    blRMSGraph->GetHistogram()->SetXTitle("Time (s)");
    blRMSGraph->GetHistogram()->SetYTitle("Baseline RMS (mV)");
    blRMSGraph->GetHistogram()->SetTitle("BaselineRMS  vs Time with line at RMS of Baselines");
    blRMSGraph->Draw("AP");
    //Use the TLine from before 
    li->SetLineColor(kBlue);
    li->DrawLine(timeArray[0], blRMS,timeArray[nSelected-1],blRMS);//line at RMS of the baselines

    //

    //remember to cleanup , for now commented out so the canvas persists
    //delete blGraph;
    //delete blRMSGraph;
    //delete li
    //delete myCan;
  }
}

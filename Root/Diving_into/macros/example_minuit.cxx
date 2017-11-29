/* Example of an extended log likelihood fit 
      control part of this macro is general   */

// -------- begin of user code -------
// global variables for this macro
TF1 *nPDF;        // probability density function for the fit
TNtuple *inpdata; //n-tuple to hold input data
// Info for initialisation of MINUIT
int NFitPar=3; // specify number of fit parameters
//------------------------------------------------------------------------------
int initialize_fit(TMinuit* minfit){ // initialisation of FIT 
 // Define a probability density function, normalized to N !
 //       exponential in range [0,5.] plus off-set
    nPDF=new TF1("eplusconstPDF","[2]*((1.-[1])*(exp(-x/[0])-exp(-5./[0]))/[0]+[1]/(5.))",0.,5.);
 // input data come from a file and are stored in an NTuple
    inpdata=new TNtuple("InputData","InputData","x");
    cout << "\nNtuple contains " << inpdata->ReadFile("expob.dat") 
         << " entries.\n\n";
    minfit->DefineParameter(0,    // Param index
                            "tau", // Param name
                            1,     // Param initial value
                            0.1,   // Param initial error, 0 for fix marameter
                            0,     // Param lower limit
                            0);    // Param upper limit
    minfit->DefineParameter(1,"offset",0.5, 0.1, 0, 1);
    minfit->DefineParameter(2,"norm",  150, 10,  0, 0);
    return 0;}
//------------------------------------------------------------------------------
//The function to be minimized, called by MINUIT, must have this form.
void the_function(Int_t &npar,                 // Optional
                  Double_t* derivatives_array, // optional
                  Double_t& function_val,      // the function value
                  Double_t* par,               // the array of parameters
                  Int_t internal_flag){        // internal flag
 // calculate extended negative log likelihood 
    function_val=0.;
    // pass on parameters to PDF
    nPDF->SetParameters(par[0],par[1],par[2]);
    // calculate -log L, i.e. loop over ntuple
    float *ntrow;
    for (int i=0; i < inpdata->GetEntries(); ++i){
        inpdata->GetEntry(i); ntrow=inpdata->GetArgs();
        function_val -= log(nPDF->Eval(ntrow[0]));}
    // add a Poission-term to take into account normalisation 
    function_val -=inpdata->GetEntries()*log(par[2])-par[2];
    function_val *=2.; // mult. by 2, as usual in ROOT, i.e. Dchi2=D(-logL) 
}
//------------------------------------------------------------------------------
void end_of_fit(TMinuit* minfit){
   // compare data with fit at the end
   TCanvas *cfit = new TCanvas("cfit","results",10,10,400,400); 
   cfit->cd();
   inpdata->Draw("x"); TH1F ht(*htemp); // access to histogram
   ht.SetLineWidth(2); ht.SetLineColor(kBlue);
   // PDF must be scaled to take into account bin width
   ht.Eval(nPDF); ht.Scale(ht.GetBinWidth(1));
   ht.SetName("Data");ht.SetTitle("Fit to data;x;N_{Events}");
   ht.DrawClone("C SAME");}
// -------- end of user code -------
// ***************************************************************************
// -------- start of general code  -------
// Function to access info on fit (and print it)
void printFit(TMinuit *minfit) {
  using namespace TMath;
  using namespace std;
    char line[200];
    Double_t vline[25];Double_t eline[25];
    cout << "\n\n\n";
    // ----------------------------------------------------------------------
    cout << "Fitted parameters: " << endl;
    cout << "       NO.     NAME       VALUE    ERROR " << endl;
    for (int n = 0; n < minfit->fNu; n++) {         
        sprintf(line, "     %4d %9s %12.5g %8.3g",n+1,
         (const char*)minfit->fCpnam[n],minfit->fU[n],minfit->fWerr[n]);
    cout << line << endl; }
    // ----------------------------------------------------------------------
    cout << "  Correlation Matrix: " << endl;
    cout << "  NO.    GLOBAL";
    for (Int_t id = 1; id <= minfit->fNu; ++id) 
      cout<< "       " <<id; cout<< endl;
        for (int i = 1; i <= minfit->fNu; ++i) {
            int ix  = minfit->fNexofi[i-1];
            int ndi = i*(i + 1) / 2;
            for (Int_t j = 1; j <= minfit->fNu; ++j){
                int m = Max(i,j); int n = Min(i,j);
                int ndex=m*(m-1)/2+n; int ndj=j*(j+1)/2;
                vline[j-1] = minfit->fVhmat[ndex-1]/
        sqrt(fabs(minfit->fVhmat[ndi-1]*minfit->fVhmat[ndj-1]));}
            sprintf(line, "  %2d   %8.3g ",ix,minfit->fGlobcc[i-1]);
            cout << line;
            for (Int_t it = 1; it <= minfit->fNu; ++it) {
                sprintf(line, "  %6.3f",vline[it-1]); cout << line; }
            cout << endl; }
    // ----------------------------------------------------------------------
    cout << "  Covariance Matrix: " << endl;
    double dxdi, dxdj;
    for (int i = 1; i <= minfit->fNu; ++i) {
        int ix  = minfit->fNexofi[i-1];
        int ndi = i*(i + 1) / 2;
        minfit->mndxdi(minfit->fX[i-1], i-1, dxdi);
        for (Int_t j = 1; j <= minfit->fNu; ++j) {
          minfit->mndxdi(minfit->fX[j-1], j-1, dxdj);
          int m=TMath::Max(i,j);int n=TMath::Min(i,j);
          int ndex=m*(m-1)/2+n;int ndj=j*(j+1)/2;
          eline[j-1] = dxdi*minfit->fVhmat[ndex-1]*dxdj*minfit->fUp; }
        for (Int_t it = 1; it <= minfit->fNu; ++it) {
         sprintf(line, " %10.3e",eline[it-1]); cout << line; }
         cout << endl; }
}
//------------------------------------------------------------------------------
void formatGraph(TGraph*g,int col,int msize,int lwidth){
     g->SetLineColor(col);    g->SetMarkerColor(col);
     g->SetMarkerSize(msize); g->SetLineWidth(lwidth);}

void plotContours(TMinuit* minfit, int p1, int p2) {
 // Get confidence contours of parameters
    int ic=0; 
    minfit->SetPrintLevel(0); // not print all countour points
    minfit->mncomd("Set ERR 4",ic); // Set the the contour level
    TGraph* cont_2sigma = (TGraph*) minfit->Contour(50);// contour w. 50 points
    minfit->mncomd("Set ERR 1",ic); // Set the the contour level
    TGraph* cont_1sigma = (TGraph*) minfit->Contour(50);// contour w. 50 points

    // The minimum of the graph and its 1 sigma error
    TGraphErrors min_g(1); min_g.SetMarkerStyle(22);
    min_g.SetPoint(0,minfit->fU[0],minfit->fU[1]);
    min_g.SetPointError(0,minfit->fWerr[0],minfit->fWerr[1]);

    // Maquillage of the Graphs
    formatGraph(cont_1sigma,kRed,p2,p1);
    formatGraph(cont_2sigma,kGreen,p2,p1);
    cont_2sigma->SetTitle("Contours;#tau;off-set");

    TCanvas *cresult = new TCanvas("cresult","",10,410,400,400);
    cresult->cd();
    cont_2sigma->DrawClone("APC");cont_1sigma->DrawClone("SamePC");
    min_g.DrawClone("PSame");}
//------------------------------------------------------------------------------
// main program for MINUIT fit
int example_minuit(){
    TMinuit* myminuit=new TMinuit(NFitPar); //initialize global pointer 
    if(initialize_fit(myminuit)!=0) return -1;
   // Standard control of a fit with MINUIT
    int ic=0;  // integer for condition code
    myminuit->SetFCN(the_function);
    myminuit->mncomd("MIN", // Start minimization (SIMPLEX first, then MIGRAD)
                    ic);    // 0 if command executed normally
    myminuit->mncomd("MINOS",ic); // Call MINOS for asymmetric errors
    myminuit->mncomd("HESSE",ic); // Call HESSE for correct error matrix
    end_of_fit(myminuit);   // Call user-defined fit summary
    printFit(myminuit);     // retrieve output from minuit
    plotContours(myminuit,2,3); //contour lines of fit parameters (2 and 3)
    return 0; }

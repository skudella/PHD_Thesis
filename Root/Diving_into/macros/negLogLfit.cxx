/* negative log likelihood fit (unbinned) */

TF1 *PDF; TNtuple *inpdata;

// --------  interface to the minimizer: function of parameters
void fFCN(Int_t &npar, double *gin=0, double &n2lL, double *par, int iflag) {
 float *ntrow;
 // calculate negative log likelihood 
 n2lL=0.;
// set parameters of PDF
 PDF->SetParameters(par[0],par[1]);
// calculate -log L, i.e. loop over ntuple
 for (int i=0; i < inpdata->GetEntries(); ++i){
    inpdata->GetEntry(i); ntrow=inpdata->GetArgs();
    n2lL -= log(PDF->Eval(ntrow[0])); 
 }
 n2lL *= 2.; //multiply by two (as common elsewhere in ROOT)
}

// ----------- main program, fit control
void negLogLfit(){

// define a probability density function, normalized to 1!
PDF=new TF1("eplusconstPDF","(1.-[1])*(exp(-x/[0])-exp(-5./[0]))/[0]+[1]/(5.)",0.,5.);
//                         exponential in range [0,5.] plus off-set 

// input data come from a file and are stored in an NTuple
inpdata=new TNtuple("InputData","InputData","x");
// read data from file and store in ntuple
ifstream inp; double x;
inp.open("expob.dat");
 while(!(inp >> x)==0){inpdata->Fill(x);}
inp.close();

// create fitter instance and initialize
 TVirtualFitter::SetDefaultFitter("Minuit");
 TVirtualFitter *fit=TVirtualFitter::Fitter(0,2);
 fit->SetFCN(fFCN); //assign function to be minimized
// set initial values of parameters
//                 #   name    val err  low up
 fit->SetParameter(0, "tau",   1., 0.1,  0, 0);
 fit->SetParameter(1, "off",   0.5, 0.1,  0, 0);

// run the fit
 double arglist[2]={5000,0.01};   // {max. number of function calls, tolerance}
  fit->ExecuteCommand("MINIMIZE", arglist, 2);      // performs SIMPLEX + MIGRAD algorithms
  fit->ExecuteCommand("MINOS", arglist, 0);         // MINOS error evaluation 

// retrieve output 
   int nvpar,nparx; double amin,edm,errdef;
   if(fit->GetStats(amin,edm,errdef,nvpar,nparx)==3){
     cout<<endl<<"*==* Fit converged:"
         << " nlL="<<amin<<" edm="<<edm<<" nvpar="<<nvpar<<" nparx="<<nparx<<endl;
   fit->PrintResults(4,amin); }
// get covariance Matrix an print it 
    TMatrixD *covMatrix = new TMatrixD(nparx,nparx,fit->GetCovarianceMatrix());
    covMatrix->Print();

//plot data, fit result, and parameter contours
   TCanvas *c = new TCanvas("c","contours",10,10,400,600);
   c->Divide(1,2);
   c->cd(1);
   inpdata->Draw("x"); 
   TH1F *ht  =(TH1F *) htemp->Clone(); // access to histogram 
   ht->SetLineWidth(2);
   ht->SetLineColor(kBlue);
   // PDF must be scaled to take account of # of Entries and bin width
   ht->Eval(PDF); ht->Scale( inpdata->GetEntries() * ht->GetBinWidth(1) );
   ht->Draw("C SAME");

// plot contours
   c->cd(2);
   //Get contour for parameter 0 versus parameter 1  for ERRDEF=2 
   gMinuit->SetErrorDef(4); //note 4 and not 2!
   TGraph *gr2 = (TGraph*)gMinuit->Contour(40,0,1);
   gr2->SetTitle("1#sigma and 2#sigma contours ;tau;off-set");
   gr2->SetFillColor(42);
   gr2->Draw("alf");
   //Get contour for parameter 0 versus parameter 2 for ERRDEF=1  
   gMinuit->SetErrorDef(1);
   TGraph *gr1 = (TGraph*)gMinuit->Contour(40,0,1);
   gr1->SetFillColor(38);
   gr1->Draw("lf");

//clean up
   delete inpdata; delete PDF;
}

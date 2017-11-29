{
//
//Menu for Tutorial "Diving into ROOT"
//

TControlBar bar("vertical");

bar.SetNameTitle("Menu Diving into ROOT", "ROOT Tutorial");

bar.AddButton("define and plot function",".x slits.cxx","plot intensity of light falling onto a double slit ");
bar.AddButton("TGrapErrors",".x macro1.cxx","display data and save as image");
bar.AddButton("display data",".x macro2.cxx","read data from file and display graph");
bar.AddButton("polar graph",".x macro3.cxx","display polar graph");
bar.AddButton("2-dim graph",".x macro4.cxx","display 2-dimensional data");
bar.AddButton("fill hisgogram",".x macro5.cxx","fill histogram with poisson numbers");
bar.AddButton("Histogram Operations",".x macro6.cxx","divide and add histograms");
bar.AddButton("2d Histograms",".x macro7.cxx","2-dimensional histograms");
bar.AddButton("write histogram",".x write_to_file.cxx","write histogram to root file");
bar.AddButton("read histogram",".x read_from_file.cxx","read histogram from root file");
bar.AddButton("save ntuple",".x write_ntuple_to_file.cxx","save n-tuple to file");
bar.AddButton("read ntuple",".x read_ntuple_from_file.cxx","read ntuple from file");
bar.AddButton("function fit",".x macro8.cxx","fit functions to data");
bar.AddButton("Toy MC",".x macro9.cxx","toy Monte Carlo: compare chi2 and logL fits");
bar.AddButton("log L fit",".x example_minuit.cxx","example of log likelihood fit based on MINUIT");
bar.AddButton("EXIT",".q","EXIT");

bar.Show();
gROOT->SaveContext();
}


/*
	Read several previously produced N-Tuples and print on screen its content

	you can easily create some files with the following statement:
	for i in 0 1 2 3 4 5; do root -l -x -b -q "write_ntuple_to_file.cxx(\"conductivity_experiment_${i}.root\", 100)"; done
*/

void read_ntuple_with_chain(){
	// initiate a TChain with the name of the TTree to be processed
	TChain in_chain("cond_data");
	in_chain.Add("conductivity_experiment*.root"); // add files, wildcards work

	// define variables and assign them to the corresponding branches
	float pot, cur, temp, pres;
	my_tuple->SetBranchAddress("Potential", &pot);
	my_tuple->SetBranchAddress("Current", &cur);
	my_tuple->SetBranchAddress("Temperature", &temp);
	my_tuple->SetBranchAddress("Pressure", &pres);

	cout << "Potential\tCurrent\tTemperature\tPressure\n";
	for (size_t irow=0; irow<in_chain.GetEntries(); ++irow){
		in_chain.GetEntry(irow); // loads all variables that have been connected to branches
		cout << pot << "\t" << cur << "\t" << temp << "\t" << pres << endl;
	}
}

//tissue.cpp
//=========================

//=========================
//Include Dependencies

//=========================
// Public Member Functions for Tissue.cpp

Tissue::Tissue(string filename) {

	ifstream ifs(filename.c_str());

	if(!ifs) {
		cout << filename << " is not available" << endl;
		return;
	}

	stringstream ss;
	string line;
	string temp;
	char trash;
	int rank;
	double height, width;
	Coord corner;
	double x, y;
	Cell* curr;

	while (getline(ifs,line)) {
		ss.str(line)

		getline(ss,temp,':');

		if (temp == "CellRank") {
			ss >> rank;
		}
		else if (temp == "Corner") {
			ss >> x >> trash >> y;
			Coord loc(x,y);
			corner = loc;
		}
		else if (temp == "Height") {
			ss >> height;
		}
		else if (temp == "Width") {
			ss >> width;
		}
		else if (temp == "End_Cell") {
			//create new cell with collected data and push onto vector 
			curr = new Cell(rank, corner, height, width);
			cells.push_back(curr);
		}

		ss.clear();
	}

	ifs.close();
}

void Tissue::calc_New_Forces() {

	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->calc_New_Forces();
	}

	return;
}

void Tissue::update_Cell_Locations() {
	
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->update_Node_Locations();
	}

	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->update_Wall_Angles();
	}

	return;
}

void print_Data_Output(ofstream & ofs) {

}

void print_VTK_File(ofstream& ofs) {

}

void grow_Cells(const int Ti) {
	
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->add_Wall_Node(Ti);
		cells.at(i)->add_Cyt_Node(Ti);
	}

	return;
}

//=========================
//End of tissue.cpp

#include "IO.hpp"
#include <iostream>
#include <fstream>

using namespace std;

IO::IO(char *input, char *output) {
	this->output = output;
	readInputfile(input);
}
IO::~IO() {

}

void IO::readInputfile(char *filename) {
	string line;

	// Default values if no input file available
	para.xLength = 1.0;
	para.yLength = 1.0;
	para.iMax = 128;
	para.jMax = 128;
	para.tEnd = 16.5;
	para.deltaX = 0.1;
	para.deltaY = 0.1;
	para.deltaT = 0;
	para.tau = 0.5;
	para.deltaVec = 1.4;
	para.iterMax = 100;
	para.eps = 0.001;
	para.omg = 1.7;
	para.alpha = 0.9;
	para.re = 1000;
	para.gx = 0;
	para.gy = 0;
	para.ui = 0;
	para.vi = 0;
	para.pi = 0;

	ifstream infile(filename);

	while (getline(infile, line)) {
		// Read values from inputvals.txt
		int pos_equ = line.find("=");
		string before_equ = line.substr(0, pos_equ);
		string after_equ = line.substr(pos_equ + 1, line.length());
		if (!before_equ.compare("xLength"))
			para.xLength = atof(after_equ.c_str());
		if (!before_equ.compare("yLength"))
			para.yLength = atof(after_equ.c_str());
		if (!before_equ.compare("iMax"))
			para.iMax = atof(after_equ.c_str());
		if (!before_equ.compare("jMax"))
			para.jMax = atof(after_equ.c_str());
		if (!before_equ.compare("tEnd"))
			para.tEnd = atof(after_equ.c_str());
		if (!before_equ.compare("tau"))
			para.tau = atof(after_equ.c_str());
		if (!before_equ.compare("deltaT"))
			para.deltaT = atof(after_equ.c_str());
		if (!before_equ.compare("deltaX"))
			para.deltaX = atof(after_equ.c_str());
		if (!before_equ.compare("deltaY"))
			para.deltaY = atof(after_equ.c_str());
		if (!before_equ.compare("deltaVec"))
			para.deltaVec = atof(after_equ.c_str());
		if (!before_equ.compare("iterMax"))
			para.iterMax = atof(after_equ.c_str());
		if (!before_equ.compare("eps"))
			para.eps = atof(after_equ.c_str());
		if (!before_equ.compare("omg"))
			para.omg = atof(after_equ.c_str());
		if (!before_equ.compare("alpha"))
			para.alpha = atof(after_equ.c_str());
		if (!before_equ.compare("re"))
			para.re = atof(after_equ.c_str());
		if (!before_equ.compare("gx"))
			para.gx = atof(after_equ.c_str());
		if (!before_equ.compare("gy"))
			para.gy = atof(after_equ.c_str());
		if (!before_equ.compare("ui"))
			para.ui = atof(after_equ.c_str());
		if (!before_equ.compare("vi"))
			para.vi = atof(after_equ.c_str());
		if (!before_equ.compare("pi"))
			para.pi = atof(after_equ.c_str());
	}

}

#define Element(field,ic) ((field)[(ic)[0]][(ic)[1]])

RealType IO::interpolateVelocityU(RealType x, RealType y, GridFunctionType & u,
		const PointType & delta) {

	RealType deltaX = delta[0];
	RealType deltaY = delta[1];

	MultiIndexType index;

	// Computation of u(x,y)
	index[0] = ((int) (x / deltaX)) + 1;
	index[1] = ((int) ((y + (deltaY / 2)) / deltaY)) + 1;

	// The coordinates of the cell corners

	RealType x1 = (index[0] - 1) * deltaX;
	RealType x2 = index[0] * deltaX;
	RealType y1 = ((index[1] - 1) - 0.5) * deltaY;
	RealType y2 = (index[1] - 0.5) * deltaY;

	MultiIndexType offset;

	offset[0] = index[0] - 1;
	offset[1] = index[1] - 1;

	RealType u1 = Element(u, offset); // datafields->u->getField ()[i - 1][j - 1];

	offset[0] = index[0];
	offset[1] = index[1] - 1;

	RealType u2 = Element(u, offset); //datafields->u->getField ()[i][j - 1];

	offset[0] = index[0] - 1;
	offset[1] = index[1];

	RealType u3 = Element(u, offset); //datafields->u->getField ()[i - 1][j];
	RealType u4 = Element(u, index);

	RealType uInterploated = (1.0 / (deltaX * deltaY))
			* ((x2 - x) * (y2 - y) * u1 + (x - x1) * (y2 - y) * u2
					+ (x2 - x) * (y - y1) * u3 + (x - x1) * (y - y1) * u4);

	return uInterploated;
}

RealType IO::interpolateVelocityV(RealType x, RealType y, GridFunctionType & v,
		const PointType & delta) {
	RealType deltaX = delta[0];
	RealType deltaY = delta[1];

	// Computation of v(x,y)
	MultiIndexType index;
	index[0] = ((int) ((x + (deltaX / 2)) / deltaX)) + 1;
	index[1] = ((int) (y / deltaY)) + 1;

	// The coordinates of the cell corners

	RealType x1 = ((index[0] - 1) - 0.5) * deltaX;
	RealType x2 = (index[0] - 0.5) * deltaX;
	RealType y1 = (index[1] - 1) * deltaY;
	RealType y2 = index[1] * deltaY;

	MultiIndexType offset;

	offset[0] = index[0] - 1;
	offset[1] = index[1] - 1;

	RealType v1 = Element(v, offset); //datafields->v->getField ()[i - 1][j - 1];

	offset[0] = index[0];
	offset[1] = index[1] - 1;

	RealType v2 = Element(v, offset); //datafields->v->getField ()[i][j - 1];

	offset[0] = index[0] - 1;
	offset[1] = index[1];

	RealType v3 = Element(v, offset); //datafields->v->getField ()[i - 1][j];

	RealType v4 = Element(v, index); //datafields->v->getField ()[i][j];

	RealType vInterpolated = (1.0 / (deltaX * deltaY))
			* ((x2 - x) * (y2 - y) * v1 + (x - x1) * (y2 - y) * v2
					+ (x2 - x) * (y - y1) * v3 + (x - x1) * (y - y1) * v4);
	return vInterpolated;
}

void IO::writeVTKFile(const 	MultiIndexType & griddimension,
								GridFunctionType u,
								GridFunctionType v,
								GridFunctionType p,
								const PointType & delta,
								int step) {
	RealType deltaX = delta[0];
	RealType deltaY = delta[1];

	IndexType iMax = griddimension[0] - 1;
	IndexType jMax = griddimension[1] - 1;

	char numstr[21];
	sprintf(numstr, "%d", step);
	std::string filename;
	filename.append("./");
	filename.append("field_");
	filename.append(numstr);
	filename.append(".vts");

	std::filebuf fb;
	fb.open(const_cast<char *>(filename.c_str()), std::ios::out);
	std::ostream os(&fb);

	os << "<?xml version=\"1.0\"?>" << std::endl
			<< "<VTKFile type=\"StructuredGrid\">" << std::endl
			<< "<StructuredGrid WholeExtent=\"" << "0" << " " << (iMax - 1)
			<< " " << "0" << " " << (jMax - 1) << " " << "0" << " " << "0"
			<< " " << "\" GhostLevel=\"" << "1" << "\">" << std::endl
			<< "<Piece Extent=\"" << "0" << " " << (iMax - 1) << " " << "0"
			<< " " << (jMax - 1) << " " << "0" << " " << "0" << " " << "\">"
			<< std::endl << "<Points>" << std::endl
			<< "<DataArray type=\"Float64\" format=\"ascii\" NumberOfComponents=\"3\"> "
			<< std::endl;
	for (int i = 0; i < iMax; ++i) {
		for (int j = 0; j < jMax; ++j) {
			os << std::scientific << i * deltaX << " " << j * deltaY << " "
					<< 0.0 << std::endl;
		}
	}
	os << "</DataArray>" << std::endl << "</Points>" << std::endl
			<< "<PointData Vectors=\"field\"  Scalars=\"P\">" << std::endl
			<< "<DataArray Name=\"field\" NumberOfComponents=\"3\" type=\"Float64\" >"
			<< std::endl;
	for (int i = 0; i < iMax; ++i) {
		RealType x = i * deltaX;

		for (int j = 0; j < jMax; ++j) {
			RealType y = j * deltaY;

			os << std::scientific << interpolateVelocityU(x, y, u, delta) << " "
					<< interpolateVelocityV(x, y, v, delta) << " " << 0.
					<< std::endl;
		}

	}
	os << "</DataArray>" << std::endl
			<< "<DataArray type=\"Float64\" Name=\"P\" format=\"ascii\">"
			<< std::endl;
	for (int i = 0; i <= iMax; ++i) {
		for (int j = 0; j <= jMax; ++j) {
			os << std::scientific << p[i][j] << " ";

		}
		os << std::endl;

	}

	os << "</DataArray>" << std::endl << "</PointData>" << std::endl
			<< "</Piece>" << std::endl << "</StructuredGrid>" << std::endl
			<< "</VTKFile>" << std::endl;
	fb.close();
}

void IO::writeVTKMasterfile(const 	MultiIndexType & griddimension,
									GridFunctionType u,
									GridFunctionType v,
									GridFunctionType p,
									const PointType & delta,
									int step) {

	IndexType iMax = griddimension[0];
	IndexType jMax = griddimension[1];

	// 2x1-Gitter
	IndexType localgriddimension[2];
	localgriddimension[0] = (griddimension[0])/2 + 1 ;
	localgriddimension[1] = griddimension[1];

	// als Argumente uebergeben ?!
	IndexType stencilwidth = 1;
	IndexType processorgridcoords[2];
	// processorgridcoords??
	processorgridcoords[0] = 0;
	processorgridcoords[1] = 1;

	char numstr[21];
	sprintf(numstr, "%d", step);
	std::string filename;
	filename.append("./");
	filename.append("field_");
	filename.append(numstr);
	filename.append(".pvtr");

	IndexType x1 = -1; //processorgridcoords[0] * (localgriddimension[0]) - stencilwidth;
	IndexType x2 = localgriddimension[0]; //(processorgridcoords[0] + 1) * (localgriddimension[0]) + stencilwidth - 1;
	IndexType x3 = -1;//processorgridcoords[1] * (localgriddimension[1]) - stencilwidth;
	IndexType x4 = localgriddimension[1]; //(processorgridcoords[1] + 1) * (localgriddimension[1]) + stencilwidth - 1;

	IndexType x5 = localgriddimension[0]; //processorgridcoords[0] * (localgriddimension[0]) - stencilwidth;
	IndexType x6 = localgriddimension[0]*2-2; //(processorgridcoords[0] + 1) * (localgriddimension[0]) + stencilwidth - 1;
	IndexType x7 = -1;//processorgridcoords[1] * (localgriddimension[1]) - stencilwidth;
	IndexType x8 = localgriddimension[1]; //(processorgridcoords[1] + 1) * (localgriddimension[1]) + stencilwidth - 1;


	std::filebuf fb;
	fb.open(const_cast<char *>(filename.c_str()), std::ios::out);
	std::ostream os(&fb);

	os << "<?xml version=\"1.0\"?>" << std::endl
	<< "<VTKFile type=\"PRectilinearGrid\">" << std::endl
	<< "<PRectilinearGrid WholeExtent=\"" << "0" << " " << (iMax - 1)
	<< " " << "0" << " " << (jMax - 1) << " " << "0" << " " << "0"
	<< " " << "\" GhostLevel=\"" << "1" << "\">" << std::endl
	<< "<PCoordinates>" << std::endl
	<< "<PDataArray type=\"Float64\">" << std::endl
	<< "<PDataArray type=\"Float64\">" << std::endl
	<< "<PDataArray type=\"Float64\">" << std::endl
	<< "</PCoordinates>" << std::endl
	<< "<Piece Extent=\"" << x1 << " " << x2 << " " << x3 << " " << x4 << " " << "0" << " " << "0"
	<< " " << "\" Source=\"" << "solution_processor_0.vtr" << "\">" << std::endl
	<< "<Piece Extent=\"" << x5 << " " << x6 << " " << x7 << " " << x8 << " " << "0" << " " << "0"
	<< " " << "\" Source=\"" << "solution_processor_1.vtr" << "\">" << std::endl
	<< "<PPointData>" << std::endl
	<< "<PDataArray type=\"Float64\" Name=\"p\">" << std::endl
	<< "</PPointData>" << std::endl
	<< "<PRectilinearGrid>" << std::endl
	<< "</VTKFile>" << std::endl;
	fb.close();

}

void IO::writeVTKSlavefile(const 	MultiIndexType & griddimension,
									GridFunctionType u,
									GridFunctionType v,
									GridFunctionType p,
									const PointType & delta,
									int world_rank,
									int d) {

	IndexType iMax = griddimension[0];
		IndexType jMax = griddimension[1];

		// 2x1-Gitter
		IndexType localgriddimension[2];
		localgriddimension[0] = (griddimension[0])/2 + 1 ;
		localgriddimension[1] = griddimension[1];

		// als Argumente uebergeben ?!
		IndexType stencilwidth = 1;
		IndexType processorgridcoords[2];
		// processorgridcoords??
		processorgridcoords[0] = 0;
		processorgridcoords[1] = 1;

	char my_rank[2];
	sprintf(my_rank, "%d", world_rank);
	std::string filename;
	filename.append("./");
	filename.append("solution_processor_");
	filename.append(my_rank);
	filename.append(".vtr");

	IndexType x1 = -1; //processorgridcoords[0] * (localgriddimension[0]) - stencilwidth;
	IndexType x2 = localgriddimension[0]; //(processorgridcoords[0] + 1) * (localgriddimension[0]) + stencilwidth - 1;
	IndexType x3 = -1;//processorgridcoords[1] * (localgriddimension[1]) - stencilwidth;
	IndexType x4 = localgriddimension[1]; //(processorgridcoords[1] + 1) * (localgriddimension[1]) + stencilwidth - 1;

	IndexType x5 = localgriddimension[0]; //processorgridcoords[0] * (localgriddimension[0]) - stencilwidth;
	IndexType x6 = localgriddimension[0]*2-2; //(processorgridcoords[0] + 1) * (localgriddimension[0]) + stencilwidth - 1;
	IndexType x7 = -1;//processorgridcoords[1] * (localgriddimension[1]) - stencilwidth;
	IndexType x8 = localgriddimension[1]; //(processorgridcoords[1] + 1) * (localgriddimension[1]) + stencilwidth - 1;

	std::filebuf fb;
	fb.open(const_cast<char *>(filename.c_str()), std::ios::out);
	std::ostream os(&fb);

	os << "<?xml version=\"1.0\"?>" << std::endl
	<< "<VTKFile type=\"PRectilinearGrid>" << std::endl
	<< "<PRectilinearGrid WholeExtent=\"" << "0" << " " << (iMax - 1)
	<< " " << "0" << " " << (jMax - 1) << " " << "0" << " " << "0"
	<< " " << "\" GhostLevel=\"" <<  "1" << "\">" << std::endl;
	if(world_rank==0){
	 	os << "<Piece Extent=\"" << x1 << " " << x2 << " " << x3 << " " << x4 << " " << "0" << " " << "0"  << "\">" << std::endl;
	}
	if(world_rank==1){
	 	os << "<Piece Extent=\"" << x5 << " " << x6 << " " << x7 << " " << x8 << " " << "0" << " " << "0"  << "\">" << std::endl;
	}
	os << "<DataArray type=\"Float64\" format=\"ascii\"> " << std::endl;
		for (int i = 0; i <= iMax-1; ++i) {
			os << std::scientific << i * delta[0] << " ";
		}
	os << "</DataArray>" << std::endl
	<< "<DataArray type=\"Float64\" format=\"ascii\"> " << std::endl;
		for (int i = 0; i <= jMax-1; ++i) {
			os << std::scientific << i * delta[1] << " ";
		}
	os << "</DataArray>" << std::endl
	<< "<DataArray type=\"Float64\" format=\"ascii\"> " << "0" << " " << "0" << " " << "</DataArray>" << std::endl
	<< "</Coordinates>" << std::endl
	<< "<PointData>" << std::endl
	<< "<PDataArray type=\"Float64\" Name=\"p\" format=\"ascii\">" << std::endl;
	for (int i = 0; i <= iMax; ++i){
		for (int j = 0; j <= jMax; ++j){
			os << std::scientific << p[i][j] << " ";
		}
		os << std::endl;
	}
	// wie/wo u und v ins file schreiben??
	os << "</DataArray>" << std::endl
	<< "</PointData>" << std::endl
	<< "</Piece>" << std::endl
	<< "</PRectilinearGrid>" << std::endl
	<< "</VTKFile>" << std::endl;
	fb.close();

}

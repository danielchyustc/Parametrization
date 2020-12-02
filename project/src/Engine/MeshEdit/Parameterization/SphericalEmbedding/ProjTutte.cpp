#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/ProjTutte.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool ProjTutte::Run() {
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::Minsurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	InitParaMap();

	UpdateTriMesh();
	return true;
}

void ProjTutte::InitParaMap()
{
	/*
	// find flat point
	double min_curv_met = 0.0;
	size_t i_flat = 0;
	for (size_t i = 0; i < nV; i++)
	{
		double K = GaussCurvature(i);
		double H = MeanCurvature(i);
		double curv_met = 2 * H * H - K;
		if (curv_met < min_curv_met)
		{
			min_curv_met = curv_met;
			i_flat = i;
		}
	}
	auto vi_flat = heMesh->Vertices().at(i_flat);

	cout << "Projected Tutte: Flattest point found." << endl;

	// compute area
//	double area_i = VoronoiArea(i_flat);
	double total_area = 0;
	for (size_t i = 0; i < nV; i++)
	{
		total_area += VoronoiArea(i);
	}
	double R = sqrt(total_area / (4 * area_i));

	cout << "Projected Tutte: Area computed." << endl;
	*/

	double area_i = 0;
	size_t i_flat = 0;
	double total_area = 0;
	for (size_t i = 0; i < nV; i++)
	{
		double area = VoronoiArea(i);
		total_area += area;
		if (area > area_i)
		{
			area_i = area;
			i_flat = i;
		}
	}
	auto vi_flat = heMesh->Vertices().at(i_flat);
	double R = sqrt(total_area / (4 * area_i));

	cout << "Projected Tutte: Area computed." << endl;

	SparseMatrix<double, RowMajor> A(nV, nV);
	MatrixXd B(nV, 2);
	SparseMatrix<double, RowMajor> zero(1, nV);
	// construct matrix A, B
	if (!cot_mat_available_)
	{
		BuildCotanWeightMatrix();
	}
	A = CotW;
	zero.setZero();
	A.row(i_flat) = zero;
	A.insert(i_flat, i_flat) = 1;
	B.setZero();

	cout << "Projected Tutte: Sparse matrix constructed" << endl;

	// set boundary to given value
	size_t d = vi_flat->Degree();
	size_t index = 0;
	for (auto vi : vi_flat->AdjVertices())
	{
		size_t i = heMesh->Index(vi);
		A.row(i) = zero;
		A.coeffRef(i, i) = 1;
		B(i, 0) = R * cos(index * 2 * 3.1415926 / d);
		B(i, 1) = R * sin(index * 2 * 3.1415926 / d);
		index++;
	}
	A.makeCompressed();

	cout << "Projected Tutte: Boundary values assigned." << endl;

	// solve sparse linear equations
	SparseLU<SparseMatrix<double>> solver;
	solver.compute(A);
	MatrixXd sol = solver.solve(B);

	cout << "Projected Tutte: Sparse linear system solved." << endl;

	// update paramap
	para_map_.resize(3, nV);
	for (size_t i = 0; i < nV; i++)
	{
		if (i == i_flat)
		{
			para_map_.col(i) = Vector3d(0, 0, 1);
		}
		else
		{
			double x = sol(i, 0);
			double y = sol(i, 1);
			double r = sqrt(x * x + y * y);
			//float theta = 4 * (1 - 1 / R) * atan(pow(r/R, 1));
			//para_map_.col(i) = Vector3f(sin(theta) * x / r, sin(theta) * y / r, -cos(theta));
			para_map_.col(i) = Vector3d(2 * x, 2 * y, r * r - 1) / (r * r + 1);
		}
	}

	cout << "Projected Tutte: Parameterization map updated." << endl;

//	MoveBarycenterToOrigin();
	init_status_ = true;
	finish_status_ = true;
}
#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/TutteEmbedding.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;


bool TutteEmbedding::Run() {
	// TODO
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::Minsurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	UpdateParaMap();
	UpdateTriMesh();
	return true;
}

void TutteEmbedding::UpdateParaMap()
{
	SparseMatrix<double> A(nV, nV);
	MatrixXd B(nV, 2);
	// construct matrix A, B
	for (size_t i = 0; i < nV; i++)
	{
		auto vi = heMesh->Vertices().at(i);
		A.insert(i, i) = 0;
		B.row(i) = Vector3d::Zero();
		for (auto vj : vi->AdjVertices())
		{
			double weight;
			switch (weight_type_)
			{
			case kEqual:
				weight = 1;
				break;
			case kCotangent:
				weight = CotanWeight(vi, vj);
				break;
			default:
				break;
			}
			A.coeffRef(i, i) += weight;
			A.insert(i, heMesh->Index(vj)) = -weight;
		}
	}

	// set boundary to given value
	size_t nB = heMesh->Boundaries()[0].size();
	for (size_t k = 0; k < nB; k++)
	{
		auto vi = heMesh->Boundaries()[0][k]->Origin();
		size_t i = heMesh->Index(vi);
		for (auto vj : vi->AdjVertices())
		{
			A.coeffRef(i, heMesh->Index(vj)) = 0;
		}
		A.coeffRef(i, i) = 1;
		double ratio;
		switch (bound_shape_)
		{
		case kCircle:
			B(i, 0) = 0.5 + 0.5 * cos(k * 2 * 3.1415926 / nB);
			B(i, 1) = 0.5 + 0.5 * sin(k * 2 * 3.1415926 / nB);
			break;
		case kSquare:
			ratio = (double)k / nB;
			if (ratio < 0.25)
			{
				B(i, 0) = ratio * 4;
				B(i, 1) = 0;
			}
			else if (ratio >= 0.25 && ratio < 0.5)
			{
				B(i, 0) = 1;
				B(i, 1) = (ratio - 0.25) * 4;
			}
			else if (ratio >= 0.5 && ratio < 0.75)
			{
				B(i, 0) = 1 - (ratio - 0.5) * 4;
				B(i, 1) = 1;
			}
			else
			{
				B(i, 0) = 0;
				B(i, 1) = 1 - (ratio - 0.75) * 4;
			}
			break;
		default:
			break;
		}
	}

	// solve sparse linear equations
	SparseLU<SparseMatrix<double>> solver;
	solver.compute(A);
	para_map_ = solver.solve(B);

	cout << "Success::Parametrize::Parametrize:" << endl
		<< "\t" << "parametrization successfully constructed" << endl;
}

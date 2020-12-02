#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/ASAP.h>

#include <Engine/Primitive/TriMesh.h>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SVD> 

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool ASAP::Run() {
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::Minsurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	InitFlatTri();
	UpdateParaMap();
	UpdateTriMesh();
	finish_status_ = true;
	return true;
}

void ASAP::UpdateTriMesh()
{
	// half-edge structure -> triangle mesh
	vector<pointf3> positions;
	vector<unsigned> indice;
	vector<normalf> normals = triMesh->GetNormals();
	vector<pointf2> texcoords;
	positions.reserve(nV);
	indice.reserve(3 * nT);
	texcoords.reserve(nV);

	double min_x = para_map_(0, 0);
	double min_y = para_map_(0, 1);
	double max_x = para_map_(0, 0);
	double max_y = para_map_(1, 0);

	for (size_t i = 0; i < nV; i++)
	{
		double x = para_map_(i, 0);
		double y = para_map_(i, 1);
		if (x < min_x)
			min_x = x;
		else if (x > max_x)
			max_x = x;

		if (y < min_y)
			min_y = y;
		else if (y > max_y)
			max_y = y;
	}

	for (size_t i = 0; i < nV; i++)
	{
		switch (view_mode_)
		{
		case kPlane:
			positions.push_back({ para_map_(i, 0), para_map_(i, 1), 0 });
			break;
		case kMesh:
			positions.push_back(heMesh->Vertices().at(i)->pos.cast_to<pointf3>());
			break;
		default:
			break;
		}
		double x = (para_map_(i, 0) - min_x) / (max_x - min_x);
		double y = (para_map_(i, 1) - min_y) / (max_y - min_y);
		texcoords.push_back({ x, y });
	}
	for (auto f : heMesh->Polygons())
	{
		for (auto v : f->BoundaryVertice())
		{
			indice.push_back(static_cast<unsigned>(heMesh->Index(v)));
		}
	}

	triMesh->Init(indice, positions, normals, texcoords);
}

void ASAP::InitFlatTri()
{
	flat_tri_ = (Matrix<double, 3, 2>*)malloc(nT * sizeof(Matrix<double, 3, 2>));
	if (flat_tri_ == nullptr)
	{
		cout << "No enough memory!" << endl;
	}

	for (size_t i = 0; i < nT; i++)
	{
		auto ti = heMesh->Polygons().at(i);
		auto vi0 = ti->HalfEdge()->Origin();
		auto vi1 = ti->HalfEdge()->End();
		auto vi2 = ti->HalfEdge()->Next()->End();
		vecf3 ei1 = vi1->pos - vi0->pos;
		vecf3 ei2 = vi2->pos - vi0->pos;
		double cos_theta = ei1.cos_theta(ei2);
		flat_tri_[i] = Matrix<double, 3, 2>();
		flat_tri_[i].row(0) = RowVector2d::Zero();
		flat_tri_[i].row(1) = RowVector2d(ei1.norm(), 0);
		flat_tri_[i].row(2) = RowVector2d(ei2.norm() * cos_theta, ei2.norm() * sqrt(1 - pow(cos_theta, 2)));
	}
}

void ASAP::UpdateParaMap()
{
	SparseMatrix<double> A(2 * nV + 2 * nT, 2 * nV + 2 * nT);
	VectorXd B(2 * nV + 2 * nT);

	int nB = heMesh->Boundaries()[0].size();
	auto v1 = heMesh->Boundaries()[0][0]->Origin();
	auto v2 = heMesh->Boundaries()[0][nB / 2]->Origin();
	int i1 = heMesh->Index(v1);
	int i2 = heMesh->Index(v2);

	A.insert(i1, i1) = 1;
	A.insert(nV + i1, nV + i1) = 1;
	A.insert(i2, i2) = 1;
	A.insert(nV + i2, nV + i2) = 1;
	B(i1) = 0;
	B(nV + i1) = 0;
	B(i2) = 1;
	B(nV + i2) = 1;

	// init A and B
	for (size_t i = 0; i < nV; i++)
	{
		if (i != i1 && i != i2)
		{
			auto vi = heMesh->Vertices().at(i);
			A.insert(i, i) = 0;
			A.insert(nV + i, nV + i) = 0;
			B(i) = 0;
			B(nV + i) = 0;
			for (auto vj : vi->AdjVertices())
			{
				int j = heMesh->Index(vj);
				A.insert(i, j) = 0;
				A.insert(nV + i, nV + j) = 0;
			}
		}
	}

	for (size_t t = 0; t < nT; t++)
	{
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
		A.insert(2 * nV + t, 2 * nV + t) = 0;
		A.insert(2 * nV + nT + t, 2 * nV + nT + t) = 0;
		B(2 * nV + t) = 0;
		B(2 * nV + nT + t) = 0;
		for (size_t k = 0; k < 3; k++, edge = edge->Next())
		{
			auto vi = edge->Origin();
			int i = heMesh->Index(vi);
			if (i != i1 && 1 != i2)
			{
				A.insert(i, 2 * nV + t) = 0;
				A.insert(i, 2 * nV + nT + t) = 0;
				A.insert(nV + i, 2 * nV + t) = 0;
				A.insert(nV + i, 2 * nV + nT + t) = 0;
			}
			A.insert(2 * nV + t, i) = 0;
			A.insert(2 * nV + t, nV + i) = 0;
			A.insert(2 * nV + nT + t, i) = 0;
			A.insert(2 * nV + nT + t, nV + i) = 0;
		}
	}

	for (size_t t = 0; t < nT; t++)
	{
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
		for (size_t k = 0; k < 3; k++, edge = edge->Next())
		{
			auto vi = edge->Origin();
			auto vj = edge->End();
			int i = heMesh->Index(vi);
			int j = heMesh->Index(vj);
			RowVector2d xi = flat_tri_[t].row(k);
			RowVector2d xj = flat_tri_[t].row((k + 1) % 3);
			double cot_ij = Cotan(vi, vj);

			if (i != i1 && i != i2)
			{
				A.coeffRef(i, i) += cot_ij;
				A.coeffRef(i, j) -= cot_ij;
				A.coeffRef(i, 2 * nV + t) -= cot_ij * (xi(0) - xj(0));
				A.coeffRef(i, 2 * nV + nT + t) -= cot_ij * (xi(1) - xj(1));

				A.coeffRef(nV + i, nV + i) += cot_ij;
				A.coeffRef(nV + i, nV + j) -= cot_ij;
				A.coeffRef(nV + i, 2 * nV + t) -= cot_ij * (xi(1) - xj(1));
				A.coeffRef(nV + i, 2 * nV + nT + t) += cot_ij * (xi(0) - xj(0));
			}

			if (j != i1 && j != i2)
			{
				A.coeffRef(j, i) -= cot_ij;
				A.coeffRef(j, j) += cot_ij;
				A.coeffRef(j, 2 * nV + t) -= cot_ij * (xj(0) - xi(0));
				A.coeffRef(j, 2 * nV + nT + t) -= cot_ij * (xj(1) - xi(1));

				A.coeffRef(nV + j, nV + i) -= cot_ij;
				A.coeffRef(nV + j, nV + j) += cot_ij;
				A.coeffRef(nV + j, 2 * nV + t) -= cot_ij * (xj(1) - xi(1));
				A.coeffRef(nV + j, 2 * nV + nT + t) += cot_ij * (xj(0) - xi(0));
			}

			A.coeffRef(2 * nV + t, i) += cot_ij * (xi(0) - xj(0));
			A.coeffRef(2 * nV + t, j) += cot_ij * (xj(0) - xi(0));
			A.coeffRef(2 * nV + t, nV + i) += cot_ij * (xi(1) - xj(1));
			A.coeffRef(2 * nV + t, nV + j) += cot_ij * (xj(1) - xi(1));
			A.coeffRef(2 * nV + t, 2 * nV + t) -= cot_ij * (xi - xj).squaredNorm();

			A.coeffRef(2 * nV + nT + t, i) += cot_ij * (xi(1) - xj(1));
			A.coeffRef(2 * nV + nT + t, j) += cot_ij * (xj(1) - xi(1));
			A.coeffRef(2 * nV + nT + t, nV + i) -= cot_ij * (xi(0) - xj(0));
			A.coeffRef(2 * nV + nT + t, nV + j) -= cot_ij * (xj(0) - xi(0));
			A.coeffRef(2 * nV + nT + t, 2 * nV + nT + t) -= cot_ij * (xi - xj).squaredNorm();
		}
	}

	SparseLU<SparseMatrix<double>> solver;
	solver.compute(A);
	VectorXd solution = solver.solve(B);

	para_map_.resize(nV, 2);
	for (size_t i = 0; i < nV; i++)
	{
		para_map_(i, 0) = solution(i);
		para_map_(i, 1) = solution(nV + i);
	}
}
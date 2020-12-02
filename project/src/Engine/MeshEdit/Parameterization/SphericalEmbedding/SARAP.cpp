#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/SARAP.h>
#include <Engine/Primitive/TriMesh.h>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SVD> 

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool SARAP::Run() {
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::Minsurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}
	cout << "1" << endl;
	InitRadius();
	cout << "2: " << radius << endl;
	InitParaMap();
	para_map_ *= radius;
	energy_ = Energy();
	cout << "3" << endl;
	InitTetra();
	cout << "4" << endl;
	InitParaSolver();
	cout << "5" << endl;
	iter_count_ = 0;

	UpdateTriMesh();
	return true;
}


bool SARAP::Iterate()
{
	iter_count_++;
	cout << iter_count_ << "th iteration" << endl;
//	UpdateTetra();
	UpdateParaMap();
//	UpdateRadius();
	UpdateTriMesh();
	return true;
}

bool SARAP::IterateTillDone()
{
	while (!finish_status_)
	{
//		UpdateTetra();
		UpdateParaMap();
//		UpdateRadius();
	}
	UpdateTriMesh();
	return true;
}

void SARAP::InitRadius()
{
	double area = 0;
	for (size_t t = 0; t < nT; t++)
	{
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
		auto v0 = edge->Origin();
		auto v1 = edge->End();
		auto v2 = edge->Next()->End();
		vecf3 e20 = v0->pos - v2->pos;
		vecf3 e21 = v1->pos - v2->pos;
		area += 0.5 * e20.cross(e21).norm();
	}
	radius = sqrt(area / (4 * 3.1415926));
}

void SARAP::InitTetra()
{
	tetra_array_ = (Matrix3d*)malloc(nT * sizeof(Matrix3d));
	if (tetra_array_ == nullptr)
	{
		cout << "No enough memory!" << endl;
	}

	for (size_t t = 0; t < nT; t++)
	{
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
		auto v0 = edge->Origin();
		auto v1 = edge->End();
		auto v2 = edge->Next()->End();
		vecf3 e01 = v1->pos - v0->pos;
		vecf3 e02 = v2->pos - v0->pos;
		vecf3 e10 = v0->pos - v1->pos;
		vecf3 e12 = v2->pos - v1->pos;
		vecf3 e20 = v0->pos - v2->pos;
		vecf3 e21 = v1->pos - v2->pos;
		double w0 = 2 * e01.cos_theta(e02) * e01.sin_theta(e02);
		double w1 = 2 * e10.cos_theta(e12) * e10.sin_theta(e12);
		double w2 = 2 * e20.cos_theta(e21) * e20.sin_theta(e21);
		vecf3 circum_center = (w0 * v0->pos + w1 * v1->pos + w2 * v2->pos) / (w0 + w1 + w2);
		double circum_radius = (circum_center - v0->pos).norm();
		double height = sqrt(pow(radius, 2) - pow(circum_radius, 2));
		tetra_array_[t] = Matrix3d();
		tetra_array_[t].col(0) = Vector3d(circum_radius, 0, height);
		tetra_array_[t].col(1) = Vector3d(circum_radius * (2 * pow(e10.cos_theta(e12), 2) - 1), circum_radius * w1, height);
		tetra_array_[t].col(2) = Vector3d(circum_radius * (2 * pow(e20.cos_theta(e21), 2) - 1), -circum_radius * w2, height);
	//	cout << tetra_array_[t] << endl;
	}
}

void SARAP::InitParaSolver()
{
	SparseMatrix<double> A(nV, nV);
//	A.insert(0, 0) = 1;
	cout << "Initializing sparse matrix..." << endl;
	// add cotangent weight
	for (size_t t = 0; t < nT; t++)
	{
		Matrix3d xt = tetra_array_[t];
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
		auto vt0 = edge->Origin();
		auto vt1 = edge->End();
		auto vt2 = edge->Next()->End();
		int i0 = heMesh->Index(vt0);
		int i1 = heMesh->Index(vt1);
		int i2 = heMesh->Index(vt2);
		double w01 = Weight(xt, 0, 1);
		double w12 = Weight(xt, 1, 2);
		double w02 = Weight(xt, 0, 2);
	//	cout << "w01: " << w01 << "	w02: " << w02 << "	w12: " << w12 << endl;
		double w0 = Weight(xt, 0);
		double w1 = Weight(xt, 1);
		double w2 = Weight(xt, 2);
	//	cout << "w0: " << w0 << "	w1: " << w1 << "	w2: " << w2 << endl;
	//	if (i0 != 0)
		{
			A.coeffRef(i0, i0) += w0 + w01 + w02;
			A.coeffRef(i0, i1) -= w01;
			A.coeffRef(i0, i2) -= w02;
		}
	//	if (i1 != 0)
		{
			A.coeffRef(i1, i1) += w1 + w01 + w12;
			A.coeffRef(i1, i0) -= w01;
			A.coeffRef(i1, i2) -= w12;
		}
	//	if (i2 != 0)
		{
			A.coeffRef(i2, i2) += w2 + w02 + w12;
			A.coeffRef(i2, i0) -= w02;
			A.coeffRef(i2, i1) -= w12;
		}
	}
	cout << "Prefactoring sparse matrix..." << endl;
	para_solver_.compute(A);
	cout << "Solver initialized." << endl;
}


void SARAP::UpdateTetra()
{
	for (size_t t = 0; t < nT; t++)
	{
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
		auto v0 = edge->Origin();
		auto v1 = edge->End();
		auto v2 = edge->Next()->End();
		vecf3 e01 = v1->pos - v0->pos;
		vecf3 e02 = v2->pos - v0->pos;
		vecf3 e10 = v0->pos - v1->pos;
		vecf3 e12 = v2->pos - v1->pos;
		vecf3 e20 = v0->pos - v2->pos;
		vecf3 e21 = v1->pos - v2->pos;
		double w0 = 2 * e01.cos_theta(e02) * e01.sin_theta(e02);
		double w1 = 2 * e10.cos_theta(e12) * e10.sin_theta(e12);
		double w2 = 2 * e20.cos_theta(e21) * e20.sin_theta(e21);
		vecf3 circum_center = (w0 * v0->pos + w1 * v1->pos + w2 * v2->pos) / (w0 + w1 + w2);
		double circum_radius = (circum_center - v0->pos).norm();
		double height = sqrt(pow(radius, 2) - pow(circum_radius, 2));
		tetra_array_[t] = Matrix3d();
		tetra_array_[t].col(0) = Vector3d(circum_radius, 0, height);
		tetra_array_[t].col(1) = Vector3d(circum_radius * (2 * pow(e10.cos_theta(e12), 2) - 1), circum_radius * w1, height);
		tetra_array_[t].col(2) = Vector3d(circum_radius * (2 * pow(e20.cos_theta(e21), 2) - 1), -circum_radius * w2, height);
	}
}

void SARAP::UpdateParaMap()
{
	MatrixXd B(nV, 3);
	B.setZero();

	for (size_t t = 0; t < nT; t++)
	{
		Matrix3d xt = tetra_array_[t];
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
		auto vt0 = edge->Origin();
		auto vt1 = edge->End();
		auto vt2 = edge->Next()->End();
		int i0 = heMesh->Index(vt0);
		int i1 = heMesh->Index(vt1);
		int i2 = heMesh->Index(vt2);
		double w01 = Weight(xt, 0, 1);
		double w12 = Weight(xt, 1, 2);
		double w02 = Weight(xt, 0, 2);
		double w0 = Weight(xt, 0);
		double w1 = Weight(xt, 1);
		double w2 = Weight(xt, 2);
		Vector3d xt0 = xt.col(0);
		Vector3d xt1 = xt.col(1);
		Vector3d xt2 = xt.col(2);
		Vector3d ut0 = para_map_.col(i0);
		Vector3d ut1 = para_map_.col(i1);
		Vector3d ut2 = para_map_.col(i2);
		Matrix3d S = w0 * xt0 * ut0.transpose()
			 + w1 * xt1 * ut1.transpose()
			 + w2 * xt2 * ut2.transpose()
			 + w01 * (xt1 - xt0) * (ut1 - ut0).transpose()
			 + w02 * (xt0 - xt2) * (ut0 - ut2).transpose()
			 + w12 * (xt2 - xt1) * (ut2 - ut1).transpose();
		JacobiSVD<Matrix3d> svd(S, ComputeFullU | ComputeFullV);
		Matrix3d U = svd.matrixU();
		Matrix3d V = svd.matrixV();
		Matrix3d L = V * U.transpose();
		B.row(i0) += (L * (w01 * (xt0 - xt1) + w02 * (xt0 - xt2) + w0 * xt0)).transpose();
		B.row(i1) += (L * (w01 * (xt1 - xt0) + w12 * (xt1 - xt2) + w1 * xt1)).transpose();
		B.row(i2) += (L * (w02 * (xt2 - xt0) + w12 * (xt2 - xt1) + w2 * xt2)).transpose();
	}

	MatrixXd solution = para_solver_.solve(B);
	for (size_t i = 0; i < nV; i++)
	{
		para_map_.col(i) = radius * solution.row(i).normalized().transpose();
	}
	double energy = Energy();
	if (fabs(energy - energy_) < 1e-5)
	{
		finish_status_ = true;
		cout << "Approaching completed." << endl;
	}
	energy_ = energy;
}

void SARAP::UpdateRadius()
{
	double area = 0;
	for (size_t t = 0; t < nT; t++)
	{
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
		auto v0 = edge->Origin();
		auto v1 = edge->End();
		auto v2 = edge->Next()->End();
		int i0 = heMesh->Index(v0);
		int i1 = heMesh->Index(v1);
		int i2 = heMesh->Index(v2);
		Vector3d e20 = para_map_.col(i0) - para_map_.col(i2);
		Vector3d e21 = para_map_.col(i1) - para_map_.col(i2);
		area += 0.5 * e20.cross(e21).norm();
	}
	radius = sqrt(area / (4 * 3.1415926));
}

double SARAP::Weight(Matrix3d xt, int i, int j)
{
	int k = 0;
	while (k == i || k == j)
		k++;
	Vector3d vi = xt.col(i);
	Vector3d vj = xt.col(j);
	Vector3d vk = xt.col(k);
//	Vector3f eki = vi - vk;
//	Vector3f ekj = vj - vk;
//	return eki.dot(ekj) / eki.cross(ekj).norm();
	Vector3d nik = vi.cross(vk).normalized();
	Vector3d njk = vj.cross(vk).normalized();
	return radius * nik.dot(njk) / nik.cross(njk).norm();
}

double SARAP::Weight(Matrix3d xt, int i)
{
	int j = 0, k = 0;
	while (j == i)
		j++;
	while (k == i || k == j)
		k++;
	Vector3d vi = xt.col(i);
	Vector3d vj = xt.col(j);
	Vector3d vk = xt.col(k);
	Vector3d nijk = (vj - vi).cross(vk - vi).normalized();
	Vector3d njk = vj.cross(vk).normalized();
	return (vj - vk).norm() * nijk.dot(njk) / nijk.cross(njk).norm();
//	return 0.1;
}

double SARAP::Energy()
{
	VectorXd mesh_area_ratio(nV);
	VectorXd para_area_ratio(nV);

	double mesh_total_area = 0;
	for (size_t i = 0; i < nV; i++)
	{
		double mesh_area = VoronoiArea(i);
		mesh_area_ratio(i) = mesh_area;
		mesh_total_area += mesh_area;
	}
	mesh_area_ratio /= mesh_total_area;

	double para_total_area = 0;
	for (size_t i = 0; i < nV; i++)
	{
		double para_area = 0;
		auto vi = heMesh->Vertices().at(i);
		for (auto vj : vi->AdjVertices())
		{
			size_t j = heMesh->Index(vj);
			size_t k = heMesh->Index(vi->HalfEdgeTo(vj)->Next()->End());
			vecd3 eij = { para_map_(0, j) - para_map_(0, i), para_map_(1, j) - para_map_(1, i), para_map_(2, j) - para_map_(2, i) };
			vecd3 eik = { para_map_(0, k) - para_map_(0, i), para_map_(1, k) - para_map_(1, i), para_map_(2, k) - para_map_(2, i) };
			para_area += eij.cross(eik).norm() / 6;
		}
		para_area_ratio(i) = para_area;
		para_total_area += para_area;
	}
	para_area_ratio /= para_total_area;

	return (para_area_ratio - mesh_area_ratio).norm();
}
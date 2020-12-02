#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/CdVA.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool CdVA::Run() {
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::Minsurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	InitParaMap();
	energy_ = Energy();
	InitSolver();
	UpdateTriMesh();
	return true;
}

bool CdVA::Iterate()
{
	for(int i = 0; i < 100; i++)
		UpdateParaMap();
	UpdateTriMesh();
	return true;
}

bool CdVA::IterateTillDone()
{
	while (!finish_status_)
		UpdateParaMap();
	UpdateTriMesh();
	return true;
}

void CdVA::InitSolver()
{
	L.resize(nV, nV);
	double factor = 2;
	cout << "Initializing sparse matrix..." << endl;
	if (!cot_mat_available_)
	{
		BuildCotanWeightMatrix();
	}
	L.setIdentity();
	L += factor * CotW;
	L.makeCompressed();
	cout << "Prefactoring sparse matrix..." << endl;
	solver.compute(L);
	cout << "Solver initialized." << endl;
}

void CdVA::UpdateParaMap()
{
	if (!finish_status_)
	{
		iter_count_++;
		cout << iter_count_ << "th iteration" << endl;

		VectorXd para_mesh_ratio(nV);
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
				para_area += eij.cross(eik).norm();
			}
			para_mesh_ratio(i) = para_area / VoronoiArea(i);
			if (para_mesh_ratio(i) > 5) para_mesh_ratio(i) = 5;
			if (para_mesh_ratio(i) < 0.2) para_mesh_ratio(i) = 0.2;
		}

		MatrixXd d = para_map_.transpose();
		for (size_t i = 0; i < nV; i++)
		{
			d.row(i) *= para_mesh_ratio(i);
		}

		d = L * d;
		for (size_t i = 0; i < nV; i++)
		{
			Vector3d si = para_map_.col(i);
			d.row(i) = d.row(i) * si * si.transpose();
		}
		MatrixXd x = solver.solve(d);
		for (size_t i = 0; i < nV; i++)
		{
			para_map_.col(i) = x.row(i).transpose() / para_mesh_ratio(i);
		//	para_map_.col(i).normalize();
		}
		Vector3d center = Barycenter(para_map_);
		for (size_t i = 0; i < nV; i++)
		{
			para_map_.col(i) = (para_map_.col(i) - center).normalized();
		}
		double energy = Energy();
		if (fabs(energy - energy_) < 1e-5)
		{
			finish_status_ = true;
			cout << "Approaching completed." << endl;
		}
		energy_ = energy;
	}
}

double CdVA::Energy()
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
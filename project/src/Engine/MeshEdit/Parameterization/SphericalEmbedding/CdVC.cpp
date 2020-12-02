#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/CdVC.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool CdVC::Run() {
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::Minsurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	InitParaMap();
	energy_ = Energy(para_map_);
	InitSolver();
	UpdateTriMesh();
	return true;
}

bool CdVC::Iterate()
{
	for(int i = 0; i < 100; i++)
		UpdateParaMap();
	UpdateTriMesh();
	return true;
}

bool CdVC::IterateTillDone()
{
	while (!finish_status_)
		UpdateParaMap();
	UpdateTriMesh();
	return true;
}

void CdVC::InitSolver()
{
	L.resize(nV, nV);
	double factor = (double)nV / 100;
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

void CdVC::UpdateParaMap()
{
	if (!finish_status_)
	{
		iter_count_++;
		cout << iter_count_ << "th iteration" << endl;
		MatrixXd d = L * para_map_.transpose();
		for (size_t i = 0; i < nV; i++)
		{
			Vector3d si = para_map_.col(i);
			d.row(i) = d.row(i) * si * si.transpose();
		}
		MatrixXd x = solver.solve(d);
		for (size_t i = 0; i < nV; i++)
		{
			para_map_.col(i) = x.row(i).transpose();
		//	para_map_.col(i).normalize();
		}
		Vector3d center = Barycenter(para_map_);
		for (size_t i = 0; i < nV; i++)
		{
			para_map_.col(i) = (para_map_.col(i) - center).normalized();
		}
		double energy = Energy(para_map_);
		if (fabs(energy - energy_) < 1e-3)
		{
			finish_status_ = true;
			cout << "Approaching completed." << endl;
		}
		energy_ = energy;
	}
}

double CdVC::Energy(RowVectorXd f)
{
	double energy = 0;
	for (size_t i = 0; i < nV; i++)
	{
		auto vi = heMesh->Vertices().at(i);
		for (auto vj : vi->AdjVertices())
		{
			int j = heMesh->Index(vj);
			energy += CotanWeight(vi, vj) * pow(f[i] - f[j], 2);
		}
	}
	return energy;
}

double CdVC::Energy(MatrixXd f)
{
	double energy = 0;
	for (size_t i = 0; i < f.rows(); i++)
	{
		energy += Energy(RowVectorXd(f.row(i)));
	}
	return energy;
}
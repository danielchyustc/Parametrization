#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/CdVB.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool CdVB::Run() {
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

bool CdVB::Iterate()
{
	for(int i = 0; i < 100; i++)
		UpdateParaMap();
	UpdateTriMesh();
	return true;
}

bool CdVB::IterateTillDone()
{
	while (!finish_status_)
		UpdateParaMap();
	UpdateTriMesh();
	return true;
}

void CdVB::InitSolver()
{
	L.resize(nV, nV);
	double factor = (double)nV / 500;
	cout << "Initializing spqarse matrix..." << endl;
	for (size_t i = 0; i < nV; i++)
	{
		auto vi = heMesh->Vertices().at(i);
		double diagonal = 1;
		for (auto vj : vi->AdjVertices())
		{
			int j = heMesh->Index(vj);
			diagonal += factor;
			L.insert(i, j) = -factor;
		}
		L.insert(i, i) = diagonal;
		if ((i + 1) % 1000 == 0)
		{
			cout << i + 1 << " points initialized" << endl;
		}
	}
	cout << "Prefactoring sparse matrix..." << endl;
	solver.compute(L);
	cout << "Solver initialized." << endl;
}

void CdVB::UpdateParaMap()
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
			para_map_.col(i).normalize();
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

double CdVB::Energy(RowVectorXd f)
{
	double energy = 0;
	for (size_t i = 0; i < nV; i++)
	{
		auto vi = heMesh->Vertices().at(i);
		for (auto vj : vi->AdjVertices())
		{
			int j = heMesh->Index(vj);
			energy += pow(f[i] - f[j], 2);
		}
	}
	return energy;
}

double CdVB::Energy(MatrixXd f)
{
	double energy = 0;
	for (size_t i = 0; i < f.rows(); i++)
	{
		energy += Energy(RowVectorXd(f.row(i)));
	}
	return energy;
}

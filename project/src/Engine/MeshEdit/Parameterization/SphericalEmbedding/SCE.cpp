#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/SCE.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool SCE::Run() {
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::SphericalEmbedding::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	InitParaMap();
	energy_ = Energy(para_map_);
	UpdateTriMesh();
	return true;
}

bool SCE::Iterate()
{
	for(int i = 0; i < 100; i++)
		UpdateParaMap();
	UpdateTriMesh();
	return true;
}

bool SCE::IterateTillDone()
{
	while (!finish_status_)
		UpdateParaMap();
	UpdateTriMesh();
	return true;
}



void SCE::UpdateParaMap()
{
	if (!finish_status_)
	{
		iter_count_++;
		cout << iter_count_ << "th iteration" << endl;
		MatrixXd delta_h(3, nV);
		for (size_t i = 0; i < nV; i++)
		{
			delta_h.col(i) = Projection(para_map_.col(i)) * LaplaceAt(para_map_, i);
		}
		para_map_ -= delta_t * delta_h;
		Vector3d barycenter = Barycenter(para_map_);
		for (size_t i = 0; i < nV; i++)
		{
			para_map_.col(i) -= barycenter;
			para_map_.col(i).normalize();
		}
		float energy = Energy(para_map_);
		if (fabs(energy - energy_) < epsilon)
		{
			finish_status_ = true;
			cout << "Approaching completed." << endl;
		}
		energy_ = energy;
	}
}

double SCE::LaplaceAt(RowVectorXd f, int i)
{
	double value = 0;
	auto vi = heMesh->Vertices().at(i);
	for (auto vj : vi->AdjVertices())
	{
		int j = heMesh->Index(vj);
		value += CotanWeight(i, vi, j, vj) * (f[i] - f[j]);
	}
	return value;
}

VectorXd SCE::LaplaceAt(MatrixXd f, int i)
{
	VectorXd value(f.rows());
	for (size_t k = 0; k < f.rows(); k++)
	{
		value[k] = LaplaceAt(RowVectorXd(f.row(k)), i);
	}
	return value;
}

MatrixXd SCE::Laplace(MatrixXd f)
{
	MatrixXd value(f.rows(), f.cols());
	for (size_t i = 0; i < f.cols(); i++)
	{
		value.col(i) = LaplaceAt(f, i);
	}
	return value;
}

double SCE::Energy(RowVectorXd f)
{
	double energy = 0;
	for (size_t i = 0; i < nV; i++)
	{
		auto vi = heMesh->Vertices().at(i);
		for (auto vj : vi->AdjVertices())
		{
			int j = heMesh->Index(vj);
			energy += CotanWeight(i, vi, j, vj) * pow(f[i] - f[j], 2);
		}
	}
	return energy;
}

double SCE::Energy(MatrixXd f)
{
	double energy = 0;
	for (size_t i = 0; i < f.rows(); i++)
	{
		energy += Energy(RowVectorXd(f.row(i)));
	}
	return energy;
}


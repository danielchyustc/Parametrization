#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/SBE.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool SBE::Run() {
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

bool SBE::Iterate()
{
	for(int i = 0; i < 100; i++)
		UpdateParaMap();
	UpdateTriMesh();
	return true;
}

bool SBE::IterateTillDone()
{
	while (!finish_status_)
		UpdateParaMap();
	UpdateTriMesh();
	return true;
}

void SBE::UpdateParaMap()
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
		for (size_t i = 0; i < nV; i++)
		{
			para_map_.col(i) -= delta_t * delta_h.col(i);
			para_map_.col(i).normalize();
		}
		double energy = Energy(para_map_);
		if (fabs(energy - energy_) < epsilon)
		{
			finish_status_ = true;
			cout << "Approaching completed." << endl;
		}
		energy_ = energy;
	}
}

double SBE::LaplaceAt(RowVectorXd f, int i)
{
	double value = 0;
	auto vi = heMesh->Vertices().at(i);
	for (auto vj : vi->AdjVertices())
	{
		int j = heMesh->Index(vj);
		value += f[i] - f[j];
	}
	return value;
}

VectorXd SBE::LaplaceAt(MatrixXd f, int i)
{
	VectorXd value(f.rows());
	for (size_t k = 0; k < f.rows(); k++)
	{
		value[k] = LaplaceAt(RowVectorXd(f.row(k)), i);
	}
	return value;
}

MatrixXd SBE::Laplace(MatrixXd f)
{
	MatrixXd value(f.rows(), f.cols());
	for (size_t i = 0; i < f.cols(); i++)
	{
		value.col(i) = LaplaceAt(f, i);
	}
	return value;
}

double SBE::Energy(RowVectorXd f)
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

double SBE::Energy(MatrixXd f)
{
	double energy = 0;
	for (size_t i = 0; i < f.rows(); i++)
	{
		energy += Energy(RowVectorXd(f.row(i)));
	}
	return energy;
}

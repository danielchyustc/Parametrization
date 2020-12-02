#include <Engine/MeshEdit/Parameterization/Parameterization.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

Parameterization::Parameterization(Ptr<TriMesh> triMesh) 
	: triMesh(triMesh), heMesh(make_shared<HEMesh<V>>())
{ 
	preproc = Preprocessing::New(triMesh);
	nV = preproc->NumEdges();
	nT = preproc->NumTriangles();
	InitHeMesh();
}

Parameterization::Parameterization(Ptr<Parameterization> para)
	: heMesh(make_shared<HEMesh<V>>())
{ 
	if (para != nullptr)
	{
		preproc = para->preproc;
		triMesh = para->triMesh;
		nV = para->nV;
		nT = para->nT;
		InitHeMesh();
		view_mode_ = para->view_mode_;
		if (para->finish_status_)
		{
			para_map_ = para->para_map_;
			init_status_ = true;
		}
		else
		{
			cout << "Warning: previous parameterization was not finished." << endl;
		}
		if (para->cot_mat_available_)
		{
			CotW = para->CotW;
			cot_mat_available_ = true;
		}
	}
	else
	{
		cout << "Error: null parameterization!" << endl;
	}
}

Parameterization::Parameterization(Ptr<Preprocessing> preproc)
	: preproc(preproc), heMesh(make_shared<HEMesh<V>>())
{
	triMesh = preproc->GetTriMesh();
	nV = preproc->NumVertices();
	nT = preproc->NumTriangles();
	InitHeMesh();
}

bool Parameterization::InitHeMesh()
{
	// init half-edge structure
	vector<vector<size_t>> triangles;
	triangles.reserve(nT);
	for (auto triangle : triMesh->GetTriangles())
	{
		triangles.push_back({ triangle->idx[0], triangle->idx[1], triangle->idx[2] });
	}
	heMesh->Reserve(nV);
	heMesh->Init(triangles);

	// positions of triangle mesh -> positions of half-edge structure
	for (size_t i = 0; i < nV; i++)
	{
		auto v = heMesh->Vertices().at(i);
		v->pos = triMesh->GetPositions()[i].cast_to<vecf3>();
	}
	return true;
}


Vector3d Parameterization::NormalAt(int i)
{
	auto vi = heMesh->Vertices().at(i);
	vecf3 normal = { 0, 0, 0 };
	for (auto vj : vi->AdjVertices())
	{
		if (vi->Degree() <= 2)
		{
			normal += vi->pos - vj->pos;
		}
		else
		{
			auto vk = vi->HalfEdgeTo(vj)->Next()->End();
			normal += (vk->pos - vi->pos).cross(vj->pos - vi->pos);
		}
	}
	normal.normalize_self();
	return Vector3d(normal[0], normal[1], normal[2]);
}

double Parameterization::Cotan(V* vi, V* vj)
{
	auto edge = vi->HalfEdgeTo(vj);
	if (edge->Polygon() == nullptr)
	{
		return 0;
	}
	else
	{
		auto vk = edge->Next()->End();
		vecf3 eki = vi->pos - vk->pos;
		vecf3 ekj = vj->pos - vk->pos;
		return abs(eki.cos_theta(ekj) / eki.sin_theta(ekj));
	}
}

double Parameterization::CotanWeight(int i, V* vi, int j, V* vj)
{
	if (cot_mat_available_)
	{
		return CotW.coeff(i, j);
	}
	else
	{
		return (Cotan(vi, vj) + Cotan(vj, vi)) / 2;
	}
}

double Parameterization::CotanWeight(int i, V* vi, V* vj)
{
	if (cot_mat_available_)
	{
		size_t j = heMesh->Index(vj);
		return CotW.coeff(i, j);
	}
	else
	{
		return (Cotan(vi, vj) + Cotan(vj, vi)) / 2;
	}
}

double Parameterization::CotanWeight(V* vi, V* vj)
{
	if (cot_mat_available_)
	{
		size_t i = heMesh->Index(vi);
		size_t j = heMesh->Index(vj);
		return CotW.coeff(i, j);
	}
	else
	{
		return (Cotan(vi, vj) + Cotan(vj, vi)) / 2;
	}
}

double Parameterization::CotanWeight(int i, int j)
{
	if (cot_mat_available_)
	{
		return CotW.coeff(i, j);
	}
	else
	{
		auto vi = heMesh->Vertices().at(i);
		auto vj = heMesh->Vertices().at(j);
		return (Cotan(vi, vj) + Cotan(vj, vi)) / 2;
	}
}


double Parameterization::SurfaceElement(V* vi)
{
	vecf3 signed_area = { 0, 0, 0 };
	for (auto vj : vi->AdjVertices())
	{
		signed_area += CircumCenter(vi, vj).cross(CircumCenter(vj, vi)) / 2;
	}
	return signed_area.norm();
}

double Parameterization::Barycenter(RowVectorXd f)
{
	double total_area = 0;
	double integral = 0;
	for (size_t i = 0; i < nV; i++)
	{
		double weight = SurfaceElement(i);
		total_area += weight;
		integral += weight * f[i];
	}
	return integral / total_area;
}

VectorXd Parameterization::Barycenter(MatrixXd f)
{
	VectorXd barrycenter(f.rows());
	for (size_t k = 0; k < f.rows(); k++)
	{
		barrycenter[k] = Barycenter(RowVectorXd(f.row(k)));
	}
	return barrycenter;
}

vecf3 Parameterization::CircumCenter(V* vi, V* vj)
{
	auto vk = vi->HalfEdgeTo(vj)->Next()->End();
	vecf3 eij = vj->pos - vi->pos;
	vecf3 eik = vk->pos - vi->pos;
	vecf3 eji = vi->pos - vj->pos;
	vecf3 ejk = vk->pos - vj->pos;
	vecf3 eki = vi->pos - vk->pos;
	vecf3 ekj = vj->pos - vk->pos;
	float weight_i = eij.sin_theta(eik) * eij.cos_theta(eik);
	float weight_j = eji.sin_theta(ejk) * eji.cos_theta(ejk);
	float weight_k = eki.sin_theta(ekj) * eki.cos_theta(ekj);
	return (vi->pos * weight_i + vj->pos * weight_j + vk->pos * weight_k) / (weight_i + weight_j + weight_k);
}

double Parameterization::TriArea(int t)
{
	auto edge = heMesh->Polygons().at(t)->HalfEdge();
	auto vi = edge->Origin();
	auto vj = edge->End();
	auto vk = edge->Next()->End();
	vecf3 vij = vj->pos - vi->pos;
	vecf3 vik = vk->pos - vi->pos;
	return vij.cross(vik).norm() / 2;
}

vecf3 Parameterization::TriNormal(int t)
{
	auto edge = heMesh->Polygons().at(t)->HalfEdge();
	auto vi = edge->Origin();
	auto vj = edge->End();
	auto vk = edge->Next()->End();
	vecf3 vij = vj->pos - vi->pos;
	vecf3 vik = vk->pos - vi->pos;
	return vij.cross(vik).normalize();
}

double Parameterization::VoronoiArea(V* vi)
{
//	return SurfaceElement(vi);
	double area = 0;
	for (auto vj : vi->AdjVertices())
	{
		auto edge = vi->HalfEdgeTo(vj);
		if (edge->Polygon() != nullptr)
		{
			auto vk = edge->Next()->End();
			area += (vj->pos - vi->pos).cross(vk->pos - vi->pos).norm();
		}
	}
	return area;
}

double Parameterization::TotalArea()
{
	double area = 0;
	for (size_t i = 0; i < nV; i++)
	{
		area += VoronoiArea(i);
	}
	return area;
}

double Parameterization::MeanCurvature(int i)
{
	vecf3 laplace = { 0, 0, 0 };
	auto vi = heMesh->Vertices().at(i);
	for (auto vj : vi->AdjVertices())
	{
		laplace += CotanWeight(i, vi, vj) * (vj->pos - vi->pos);
	}
	return 0.5 * laplace.norm() / VoronoiArea(i);
}

double Parameterization::GaussCurvature(int i)
{
	double K = 2 * 3.1415926;
	auto vi = heMesh->Vertices().at(i);
	for (auto vj : vi->AdjVertices())
	{
		auto vk = vi->HalfEdgeTo(vj)->Next()->End();
		double theta = acos((vj->pos - vi->pos).cos_theta(vk->pos - vi->pos));
		K = K - theta;
	}
	return K / VoronoiArea(i);
}

void Parameterization::BuildCotanWeightMatrix()
{
	cout << "Building cotangent weight matrix..." << endl;
	CotW = SparseMatrix<double, RowMajor>(nV, nV);
	for (size_t i = 0; i < nV; i++)
	{
		CotW.insert(i, i) = 0;
		auto vi = heMesh->Vertices().at(i);
		for (auto vj : vi->AdjVertices())
		{
			size_t j = heMesh->Index(vj);
			double weight = (Cotan(vi, vj) + Cotan(vj, vi)) / 2;
			CotW.coeffRef(i, i) += weight;
			CotW.insert(i, j) = -weight;
		}
		if ((i + 1) % 10000 == 0)
		{
			cout << i + 1 << " rows built." << endl;
		}
	}
	CotW.makeCompressed();
	cot_mat_available_ = true;
	cout << "Cotangent weight matrix built." << endl;
}

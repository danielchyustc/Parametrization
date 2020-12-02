#include <Engine/MeshEdit/MinSurf.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

MinSurf::MinSurf(Ptr<TriMesh> triMesh)
	: heMesh(make_shared<HEMesh<V>>()), weight_type_(kEqual)
{
	Init(triMesh);
}

MinSurf::MinSurf(Ptr<TriMesh> triMesh, WeightType weight_type)
	: heMesh(make_shared<HEMesh<V>>()), weight_type_(weight_type)
{
	Init(triMesh);
}

void MinSurf::Clear() {
	heMesh->Clear();
	triMesh = nullptr;
}

bool MinSurf::Init(Ptr<TriMesh> triMesh) {
	Clear();

	if (triMesh == nullptr)
		return true;

	if (triMesh->GetType() == TriMesh::INVALID) {
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is invalid\n");
		return false;
	}

	// init half-edge structure
	size_t nV = triMesh->GetPositions().size();
	vector<vector<size_t>> triangles;
	triangles.reserve(triMesh->GetTriangles().size());
	for (auto triangle : triMesh->GetTriangles())
		triangles.push_back({ triangle->idx[0], triangle->idx[1], triangle->idx[2] });
	heMesh->Reserve(nV);
	heMesh->Init(triangles);

	if (!heMesh->IsTriMesh() || !heMesh->HaveBoundary()) {
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is not a triangle mesh or hasn't a boundaries\n");
		heMesh->Clear();
		return false;
	}

	// triangle mesh's positions ->  half-edge structure's positions
	for (int i = 0; i < nV; i++) {
		auto v = heMesh->Vertices().at(i);
		v->pos = triMesh->GetPositions()[i].cast_to<vecf3>();
	}

	this->triMesh = triMesh;
	return true;
}

bool MinSurf::Run() {
	if (heMesh->IsEmpty() || !triMesh) {
		printf("ERROR::MinSurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	Minimize();

	// half-edge structure -> triangle mesh
	size_t nV = heMesh->NumVertices();
	size_t nF = heMesh->NumPolygons();
	vector<pointf3> positions;
	vector<unsigned> indice;
	vector<normalf> normals = triMesh->GetNormals();
	vector<pointf2> texcoords = triMesh->GetTexcoords();
	positions.reserve(nV);
	indice.reserve(3 * nF);
	for (auto v : heMesh->Vertices())
		positions.push_back(v->pos.cast_to<pointf3>());
	for (auto f : heMesh->Polygons()) { // f is triangle
		for (auto v : f->BoundaryVertice()) // vertices of the triangle
			indice.push_back(static_cast<unsigned>(heMesh->Index(v)));
	}

	triMesh->Init(indice, positions, normals, texcoords);

	return true;
}

void MinSurf::Minimize() {
	// TODO
	size_t nV = heMesh->Vertices().size();
	SparseMatrix<float> A(nV, nV);
	MatrixXf B(nV, 3);
	// construct matrix A, B
	for (size_t i = 0; i < nV; i++)
	{
		auto vi = heMesh->Vertices().at(i);
		if (vi->IsBoundary())
		{
			A.insert(i, i) = 1;
			B(i, 0) = vi->pos.at(0);
			B(i, 1) = vi->pos.at(1);
			B(i, 2) = vi->pos.at(2);
		}
		else
		{
			A.insert(i, i) = 0;
			B.row(i) = Vector3f::Zero();
			for (auto vj : vi->AdjVertices())
			{
				float weight;
				switch (weight_type_)
				{
				case kEqual:
					weight = 1;
					break;
				case kCotangent:
					weight = CotangentWeight(vi, vj);
					break;
				default:
					break;
				}
				A.coeffRef(i, i) += weight;
				A.insert(i, heMesh->Index(vj)) = -weight;
			}
		}
	}

	// solve sparse linear equations
	SparseLU<SparseMatrix<float>> solver;
	solver.compute(A);
	MatrixXf solution = solver.solve(B);
	for (int i = 0; i < nV; i++)
	{
		auto vi = heMesh->Vertices().at(i);
		vi->pos.at(0) = solution(i, 0);
		vi->pos.at(1) = solution(i, 1);
		vi->pos.at(2) = solution(i, 2);
	}

	cout << "Success::MinSurf::Minimize:" << endl
		<< "\t" << "minimal surface successfully constructed" << endl;
}

float MinSurf::CotangentWeight(V* vi, V* vj)
{
	V* v_pre = vi->HalfEdgeTo(vj)->RotatePre()->End();
	V* v_next = vi->HalfEdgeTo(vj)->RotateNext()->End();
	vecf3 v_pi = vi->pos - v_pre->pos;
	vecf3 v_pj = vj->pos - v_pre->pos;
	vecf3 v_ni = vi->pos - v_next->pos;
	vecf3 v_nj = vj->pos - v_next->pos;
	float cot_alpha = v_pi.cos_theta(v_pj) / v_pi.sin_theta(v_pj);
	float cot_beta = v_ni.cos_theta(v_nj) / v_ni.sin_theta(v_nj);
	float weight = cot_alpha + cot_beta;
	return weight > 0 ? weight : 0;
}
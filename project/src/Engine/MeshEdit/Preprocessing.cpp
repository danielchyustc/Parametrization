#include <Engine/MeshEdit/Preprocessing.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

Preprocessing::Preprocessing(Ptr<TriMesh> triMesh)
	: heMesh(make_shared<HEMesh<V>>())
{
	Init(triMesh);
}

void Preprocessing::Clear() {
	heMesh->Clear();
	triMesh = nullptr;
}

bool Preprocessing::Init(Ptr<TriMesh> triMesh) {
	Clear();
	Glue(triMesh);

	if (triMesh == nullptr)
		return true;

	if (triMesh->GetType() == TriMesh::INVALID) {
		printf("ERROR::Preprocessing::Init:\n"
			"\t""trimesh is invalid\n");
		return false;
	}

	nV = triMesh->GetPositions().size();
	nT = triMesh->GetTriangles().size();

	// init half-edge structure
	vector<vector<size_t>> triangles;
	triangles.reserve(nT);
	for (auto triangle : triMesh->GetTriangles())
		triangles.push_back({ triangle->idx[0], triangle->idx[1], triangle->idx[2] });
	heMesh->Reserve(nV);
	heMesh->Init(triangles);

	nE = heMesh->NumEdges();
	nB = heMesh->NumBoundaries();
	euler_char = nV - nE + nT;
	genus = (2 - euler_char - nB) / 2;
	isTriMesh = heMesh->IsTriMesh();
	isClosed = !heMesh->HaveBoundary();

	if (!isTriMesh) {
		printf("ERROR::Preprocessing::Init:\n"
			"\t""trimesh is not a triangle mesh!\n");
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

bool Preprocessing::Glue(Ptr<TriMesh> triMesh) {
	if (triMesh == nullptr || triMesh->GetType() == TriMesh::INVALID) {
		printf("ERROR::Glue::Run:\n"
			"\t""triMesh == nullptr || triMesh->GetType() == TriMesh::INVALID\n");
		return false;
	}

	auto& positions = triMesh->GetPositions();
	vector<unsigned> indice;
	vector<pointf3> uniquePos;
	vector<normalf> normals;
	vector<pointf2> texcoords;
	map<pointf3, unsigned> pos2idx;

	for (auto triangle : triMesh->GetTriangles()) {
		for (int i = 0; i < 3; i++) {
			auto pos = positions[triangle->idx[i]];
			auto target = pos2idx.find(pos);
			if (target == pos2idx.end()) {
				pos2idx[pos] = static_cast<unsigned>(uniquePos.size());
				uniquePos.push_back(pos);
				normals.push_back(triMesh->GetNormals()[triangle->idx[i]]);
				texcoords.push_back(triMesh->GetTexcoords()[triangle->idx[i]]);
				target = pos2idx.find(pos);
			}
			indice.push_back(target->second);
		}
	}
	printf("Glue: %zd -> %zd\n", positions.size(), uniquePos.size());

	triMesh->Init(indice, uniquePos, normals, texcoords);

	return true;
}
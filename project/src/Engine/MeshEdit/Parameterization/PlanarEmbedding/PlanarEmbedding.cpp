#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/PlanarEmbedding.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;


PlanarEmbedding::PlanarEmbedding(Ptr<TriMesh> triMesh)
	: Parameterization(triMesh)
{
	if (preproc->IsClosed())
	{
		cout << "Error: the mesh does not have boundaries!" << endl;
	}
}

PlanarEmbedding::PlanarEmbedding(Ptr<Parameterization> embd)
	: Parameterization(embd)
{
	if (preproc->IsClosed())
	{
		cout << "Error: the mesh does not have boundaries!" << endl;
	}
}

PlanarEmbedding::PlanarEmbedding(Ptr<Preprocessing> preproc)
	: Parameterization(preproc)
{
	if (preproc->IsClosed())
	{
		cout << "Error: the mesh does not have boundaries!" << endl;
	}
}

bool PlanarEmbedding::Run() {
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::PlanarEmbedding::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	UpdateParaMap();
	UpdateTriMesh();
	
	return true;
}

void PlanarEmbedding::UpdateTriMesh()
{
	// half-edge structure -> triangle mesh
	size_t nV = heMesh->NumVertices();
	size_t nF = heMesh->NumPolygons();
	vector<pointf3> positions;
	vector<unsigned> indice;
	vector<normalf> normals = triMesh->GetNormals();
	vector<pointf2> texcoords;
	positions.reserve(nV);
	indice.reserve(3 * nF);
	texcoords.reserve(nV);
	for (size_t i = 0; i < nV; i++)
	{
		positions.push_back(heMesh->Vertices().at(i)->pos.cast_to<pointf3>());
		texcoords.push_back({ para_map_(i, 0), para_map_(i, 1) });
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
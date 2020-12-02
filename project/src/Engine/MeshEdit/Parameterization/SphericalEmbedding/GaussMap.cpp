#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/GaussMap.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool GaussMap::Run() {
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::Minsurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	InitParaMap();

	UpdateTriMesh();
	return true;
}

void GaussMap::InitParaMap()
{
	para_map_.resize(3, nV);
	for (size_t i = 0; i < nV; i++)
	{
		para_map_.col(i) = NormalAt(i);
	}
	iter_count_ = 0;
	finish_status_ = true;
}
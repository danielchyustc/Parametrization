#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/SphericalEmbedding.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;


SphericalEmbedding::SphericalEmbedding(Ptr<TriMesh> triMesh)
	: Parameterization(triMesh)
{
	if (!preproc->IsClosed())
	{
		cout << "Error: the mesh is not closed!" << endl;
	}
}

SphericalEmbedding::SphericalEmbedding(Ptr<Parameterization> embd)
	: Parameterization(embd)
{
	if (!preproc->IsClosed())
	{
		cout << "Error: the mesh is not closed!" << endl;
	}
}

SphericalEmbedding::SphericalEmbedding(Ptr<Preprocessing> preproc)
	: Parameterization(preproc)
{
	if (!preproc->IsClosed())
	{
		cout << "Error: the mesh is not closed!" << endl;
	}
}

bool SphericalEmbedding::Run() {
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::SphericalEmbedding::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	InitParaMap();

	UpdateTriMesh();
	return true;
}

bool SphericalEmbedding::Iterate()
{
	for(int i = 0; i < 100; i++)
		UpdateParaMap();
	UpdateTriMesh();
	return true;
}

bool SphericalEmbedding::IterateTillDone()
{
	while (!finish_status_)
		UpdateParaMap();
	UpdateTriMesh();
	return true;
}

void SphericalEmbedding::InitParaMap()
{
	if (!init_status_)
	{
		para_map_.resize(3, nV);
		for (size_t i = 0; i < nV; i++)
		{
			para_map_.col(i) = NormalAt(i);
		}
		init_status_ = true;
	}
}

void SphericalEmbedding::UpdateTriMesh()
{
	// half-edge structure -> triangle mesh
	size_t nF = heMesh->NumPolygons();
	vector<pointf3> positions;
	vector<unsigned> indice;
	vector<normalf> normals;
	vector<pointf2> texcoords;
	positions.reserve(nV);
	indice.reserve(3 * nF);
	normals.reserve(nV);
	texcoords.reserve(nV);
	for (size_t i = 0; i < nV; i++)
	{
		Vector2d projection = StereoProjection(para_map_.col(i));
		switch (view_mode_)
		{
		case kMesh:
			positions.push_back(heMesh->Vertices().at(i)->pos.cast_to<pointf3>());
			break;
		case kSphere:
			positions.push_back({ para_map_(0, i), para_map_(1, i), para_map_(2, i) });
			break;
		case kPlane:
			positions.push_back({ projection[0], projection[1], 0 });
			break;
		default:
			break;
		}
		float x = atan2(para_map_(1, i), para_map_(0, i)) / 3.1415926;
		float y = 0.5 + asin(para_map_(2, i)) / 3.1415926;
		texcoords.push_back({ 1 - fabs(x), y });
		Vector3d normal = NormalAt(i);
		normals.push_back({ normal[0], normal[1], normal[2] });
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

Matrix3d SphericalEmbedding::Projection(Vector3d v)
{
	return Matrix3d::Identity() - v * v.transpose() / v.squaredNorm();
}

void SphericalEmbedding::MoveBarycenterToOrigin()
{
	Vector3d center = Barycenter(para_map_);
	size_t move_count = 0;
	while (center.norm() > 1e-5 && move_count < 100)
	{
		for (size_t i = 0; i < nV; i++)
		{
			para_map_.col(i) -= center;
			para_map_.col(i).normalize();
		}
		center = Barycenter(para_map_);
		move_count++;
		if (move_count % 10 == 0)
		{
			cout << move_count << " th move to the origin." << endl;
		}
	}

	if (move_count == 0)
	{
		for (size_t i = 0; i < nV; i++)
		{
			para_map_.col(i).normalize();
		}
	}

	cout << "Spherical Embedding: Barycenter moved to the origin." << move_count << " moves in total." << endl;

}

void SphericalEmbedding::AreaCorrection()
{
//	VectorXd mesh_area_ratio(nV);
//	VectorXd para_area_ratio(nV);
	
/*
	double mesh_total_area = 0;
	for (size_t i = 0; i < nV; i++)
	{
		double mesh_area = VoronoiArea(i);
		mesh_area_ratio(i) = mesh_area;
		mesh_total_area += mesh_area;
	}
	mesh_area_ratio /= mesh_total_area;

	size_t correction_count = 0;
//	do{
//	while(correction_count < 10)
//	{
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
*/

	VectorXd para_mesh_ratio(nV);
	for (size_t count = 0; count < 10; count++)
	{
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
		}
//		cout << "Spherical Embedding: para mesh ratio computed." << endl;

		MatrixXd new_map(3, nV);
		new_map.setZero();
		for (size_t i = 0; i < nV; i++)
		{
			auto vi = heMesh->Vertices().at(i);
			for (auto vj : vi->AdjVertices())
			{
				size_t j = heMesh->Index(vj);
				new_map.col(i) += para_mesh_ratio(j) * CotW.coeff(i, j) * para_map_.col(j);
			}
			new_map.col(i).normalize();
		}
		para_map_ = new_map;
	/*
		Vector3d center = Barycenter(para_map_);
		for (size_t i = 0; i < nV; i++)
		{
			para_map_.col(i) -= center;
			para_map_.col(i).normalize();
		}
		*/
		MoveBarycenterToOrigin();
	}

	cout << "10 times." << endl;
/*
		for (size_t i = 0; i < nV; i++)
		{
			double ratio = mesh_area_ratio(i) / para_area_ratio(i);
			ratio = ratio < 10 ? ratio : 10;
			ratio = ratio > 0.1 ? ratio : 0.1;
			para_map_.col(i) *= ratio;
		}
*/
	//	MoveBarycenterToOrigin();

	//	correction_count++;
	//	if (correction_count % 10 == 0)
	//	{
	//		cout << correction_count << " th area correction done." << endl;
	//	}
//	}
//	} while ((para_area_ratio - mesh_area_ratio).norm() > 1e-5);
}
#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/ARAP.h>
#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/ASAP.h>
#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SVD> 

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool ARAP::Run() {
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::Minsurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}
	
	if (!init_status_)
	{
		cout << "Error: need initialization." << endl;
		return false;
	}

	InitPara();
	InitFlatTri();
	InitParaSolver();
	iter_count_ = 0;

	UpdateTriMesh();
	return true;
}

void ARAP::UpdateTriMesh()
{
	// half-edge structure -> triangle mesh
	vector<pointf3> positions;
	vector<unsigned> indice;
	vector<normalf> normals = triMesh->GetNormals();
	vector<pointf2> texcoords;
	positions.reserve(nV);
	indice.reserve(3 * nT);
	texcoords.reserve(nV);

	double min_x = para_map_(0, 0);
	double min_y = para_map_(0, 1);
	double max_x = para_map_(0, 0);
	double max_y = para_map_(1, 0);

	for (size_t i = 0; i < nV; i++)
	{
		double x = para_map_(i, 0);
		double y = para_map_(i, 1);
		if (x < min_x)
			min_x = x;
		else if (x > max_x)
			max_x = x;

		if (y < min_y)
			min_y = y;
		else if (y > max_y)
			max_y = y;
	}

	for (size_t i = 0; i < nV; i++)
	{
		switch (view_mode_)
		{
		case kPlane:
			positions.push_back({ para_map_(i, 0), para_map_(i, 1), 0 });
			break;
		case kMesh:
			positions.push_back(heMesh->Vertices().at(i)->pos.cast_to<pointf3>());
			break;
		default:
			break;
		}
		double x = (para_map_(i, 0) - min_x) / (max_x - min_x);
		double y = (para_map_(i, 1) - min_y) / (max_y - min_y);
		texcoords.push_back({ x, y });
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

void ARAP::InitFlatTri()
{
	flat_tri_ = (Matrix<double, 3, 2>*)malloc(nT * sizeof(Matrix<double, 3, 2>));
	if (flat_tri_ == nullptr)
	{
		cout << "No enough memory!" << endl;
	}

	for (size_t i = 0; i < nT; i++)
	{
		auto ti = heMesh->Polygons().at(i);
		auto vi0 = ti->HalfEdge()->Origin();
		auto vi1 = ti->HalfEdge()->End();
		auto vi2 = ti->HalfEdge()->Next()->End();
		vecf3 ei1 = vi1->pos - vi0->pos;
		vecf3 ei2 = vi2->pos - vi0->pos;
		double cos_theta = ei1.cos_theta(ei2);
		flat_tri_[i] = Matrix<double, 3, 2>();
		flat_tri_[i].row(0) = RowVector2d::Zero();
		flat_tri_[i].row(1) = RowVector2d(ei1.norm(), 0);
		flat_tri_[i].row(2) = RowVector2d(ei2.norm() * cos_theta, ei2.norm() * sqrt(1 - pow(cos_theta, 2)));
	}
}

void ARAP::InitParaSolver()
{
	SparseMatrix<double> A(nV, nV);
	for (size_t i = 0; i < nV; i++)
	{
		if (i == 0)
		{
			A.insert(i, i) = 1;
		}
		else
		{
			auto vi = heMesh->Vertices().at(i);
			A.insert(i, i) = 0;
			for (auto vj : vi->AdjVertices())
			{
				A.insert(i, heMesh->Index(vj)) = 0;
			}
		}
	}

	// add cotangent weight
	for (size_t t = 0; t < nT; t++)
	{
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
		for (int k = 0; k < 3; k++, edge = edge->Next())
		{
			auto vi = edge->Origin();
			auto vj = edge->End();
			int i = heMesh->Index(vi);
			int j = heMesh->Index(vj);
			double cot_theta = Cotan(vi, vj);
			if (i != 0)
			{
				A.coeffRef(i, i) += cot_theta;
				A.coeffRef(i, j) -= cot_theta;
			}
			if (j != 0)
			{
				A.coeffRef(j, i) -= cot_theta;
				A.coeffRef(j, j) += cot_theta;
			}
		}
	}
	para_solver_.compute(A);
}

bool ARAP::Iterate()
{
	iter_count_++;
	cout << iter_count_ << "th iteration" << endl;
	UpdateParaMap();
	UpdateTriMesh();
	return true;
}


void ARAP::UpdateParaMap()
{
	MatrixXd B(nV, 2);
	for (size_t i = 0; i < nV; i++)
	{
		B.row(i) = RowVector2d::Zero();
	}

	for (size_t t = 0; t < nT; t++)
	{
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
		auto vt0 = edge->Origin();
		auto vt1 = edge->End();
		auto vt2 = edge->Next()->End();
		Matrix<double, 3, 2> xt = flat_tri_[t];
		Vector2d ut0 = para_map_.row(heMesh->Index(vt0)).transpose();
		Vector2d ut1 = para_map_.row(heMesh->Index(vt1)).transpose();
		Vector2d ut2 = para_map_.row(heMesh->Index(vt2)).transpose();
		Matrix2d S = Cotan(vt0, vt1) * (ut0 - ut1) * (xt.row(0) - xt.row(1))
			+ Cotan(vt1, vt2) * (ut1 - ut2) * (xt.row(1) - xt.row(2))
			+ Cotan(vt2, vt0) * (ut2 - ut0) * (xt.row(2) - xt.row(0));
		JacobiSVD<Matrix2d> svd(S, ComputeFullU | ComputeFullV);
		Matrix2d U = svd.matrixU();
		Matrix2d V = svd.matrixV();
		Matrix2d L = V * U.transpose();
	//	Matrix2f L;
	//	L(0, 0) = cos(atan2(S(0, 1) - S(1, 0), S.trace()));
	//	L(0, 1) = -sin(atan2(S(0, 1) - S(1, 0), S.trace()));
	//	L(1, 0) = sin(atan2(S(0, 1) - S(1, 0), S.trace()));
	//	L(1, 1) = cos(atan2(S(0, 1) - S(1, 0), S.trace()));

		for (int k = 0; k < 3; k++, edge = edge->Next())
		{
			auto vi = edge->Origin();
			auto vj = edge->End();
			int i = heMesh->Index(vi);
			int j = heMesh->Index(vj);
			double cot_theta = Cotan(vi, vj);
			B.row(i) += (xt.row(k) - xt.row((k + 1) % 3)) * L * cot_theta;
			B.row(j) += (xt.row((k + 1) % 3) - xt.row(k)) * L * cot_theta;
		}
	}

	para_map_ = para_solver_.solve(B);
}
#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/VLG.h>
#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SVD> 

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool VLG::Run() {
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::Minsurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	InitPara();
	iter_count_ = 0;
	init_status_ = true;
	finish_status_ = false;

	UpdateTriMesh();
	return true;
}

void VLG::InitPara()
{
	F0.resize(nT);
	F1.resize(nT);
	F2.resize(nT);
	input_mesh_pos.resize(3, nV);
	position.resize(2 * nV);
	velocity.resize(2 * nV);
	para_map_.resize(nV, 2);

	area_weight.resize(nT);
	double area_sum = 0.0;
	for (size_t t = 0; t < nT; t++)
	{
		double area_t = TriArea(t);
		area_weight[t] = area_t;
		area_sum += area_t;
	}
	area_weight /= area_sum;
	double area_coeff = 1 / sqrt(area_sum);
	for (size_t i = 0; i < nV; i++)
	{
		vecf3 pos = heMesh->Vertices().at(i)->pos;
		input_mesh_pos.col(i) = Vector3d(pos[0], pos[1], pos[2]) * area_coeff;
	}
	velocity.setZero();

	for (size_t t = 0; t < nT; t++)
	{
		auto edge = heMesh->Polygons().at(t)->HalfEdge();
		auto vi = edge->Origin();
		auto vj = edge->End();
		auto vk = edge->Next()->End();
		F0[t] = heMesh->Index(vi);
		F1[t] = heMesh->Index(vj);
		F2[t] = heMesh->Index(vk);
	}
	source_lc.resize(nT, 3);
	source_local_coordinate();
	Tutte();
	target_lc.resize(nT, 3);
	target_local_coordinate();
}

void VLG::source_local_coordinate()
{
	int f0, f1, f2;
	Vector3d e01, e02, n_, x_, y_;
	for (size_t t = 0; t < nT; t++)
	{
		f0 = F0[t];
		f1 = F1[t];
		f2 = F2[t];

		e01 = input_mesh_pos.col(f1) - input_mesh_pos.col(f0);
		e02 = input_mesh_pos.col(f2) - input_mesh_pos.col(f0);
		n_ = e01.cross(e02).normalized();
		x_ = e01.normalized();
		y_ = n_.cross(x_);

		source_lc(t, 0) = e01.norm();
		source_lc(t, 1) = e02.dot(x_);
		source_lc(t, 2) = e02.dot(y_);
	}
}

void VLG::target_local_coordinate()
{
	int f0, f1, f2;
	double x0, x1, x2, y0, y1, y2;
	double x, y, z;
	for (size_t t = 0; t < nT; t++)
	{
		f0 = F0[t];
		f1 = F1[t];
		f2 = F2[t];

		x0 = position[f0];	y0 = position[f0 + nV];
		x1 = position[f1];	y1 = position[f1 + nV];
		x2 = position[f2];	y2 = position[f2 + nV];

		x = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		y = ((x1 - x0) * (x2 - x0) + (y1 - y0) * (y2 - y0)) / x;
		z = ((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)) / x;

		target_lc(t, 0) = x;
		target_lc(t, 1) = y;
		target_lc(t, 2) = z;
	}
}

void VLG::Tutte()
{
	SparseMatrix<double> A(nV, nV);
	MatrixXd B(nV, 2);
	// construct matrix A, B
	for (size_t i = 0; i < nV; i++)
	{
		auto vi = heMesh->Vertices().at(i);
		A.insert(i, i) = 0;
		B.row(i) = Vector2d::Zero();
		for (auto vj : vi->AdjVertices())
		{
			A.coeffRef(i, i) += 1;
			A.insert(i, heMesh->Index(vj)) = -1;
		}
	}

	// set boundary to given value
	size_t nB = heMesh->Boundaries()[0].size();
	double area_coeff = 1 / sqrt(3.1415926);
	for (size_t k = 0; k < nB; k++)
	{
		auto vi = heMesh->Boundaries()[0][k]->Origin();
		size_t i = heMesh->Index(vi);
		for (auto vj : vi->AdjVertices())
		{
			A.coeffRef(i, heMesh->Index(vj)) = 0;
		}
		A.coeffRef(i, i) = 1;
		B(i, 0) = cos(k * 2 * 3.1415926 / nB) * area_coeff;
		B(i, 1) = -sin(k * 2 * 3.1415926 / nB) * area_coeff;
	}

	// solve sparse linear equations
	SparseLU<SparseMatrix<double>> solver;
	solver.compute(A);
	para_map_ = solver.solve(B);

	for (size_t i = 0; i < nV; i++)
	{
		position(i) = para_map_(i, 0);
		position(i + nV) = para_map_(i, 1);
	}
}

bool VLG::Iterate()
{
	if (iter_count_ >= 10000)
	{
		cout << "Have reached maximal iteration number!" << endl;
		return true;
	}
	if (finish_status_)
	{
		cout << "Already finished!" << endl;
		return true;
	}

	solve_velocity();
	position += velocity;
	target_local_coordinate();
	iter_count_++;

	cout << iter_count_ << "th iteration done!" << endl;
	for (size_t i = 0; i < nV; i++)
	{
		para_map_(i, 0) = position(i);
		para_map_(i, 1) = position(i + nV);
	}
	UpdateTriMesh();
	return true;
}

void VLG::solve_velocity()
{
	SparseMatrix<double, RowMajor> A(2 * nV, 2 * nV);
	VectorXd B(2 * nV);
	int f0, f1, f2;
	double x0, x1, x2, y0, y1, y2;
	double x, y, z, delta_x, delta_y, delta_z;
	VectorXd a0(6), a1(6), a2(6), a3(6), u(6), v(6), w(6), b(6);
	MatrixXd h(6, 6);
	VectorXi index(6);
	int t_fix = 0;
	double min_distortion = 100;
	A.setZero();
	B.setZero();
	for (size_t t = 0; t < nT; t++)
	{
		f0 = F0[t];
		f1 = F1[t];
		f2 = F2[t];

		x0 = position[f0];	y0 = position[f0 + nV];
		x1 = position[f1];	y1 = position[f1 + nV];
		x2 = position[f2];	y2 = position[f2 + nV];

		x = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		y = ((x1 - x0) * (x2 - x0) + (y1 - y0) * (y2 - y0)) / x;
		z = ((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)) / x;

		delta_x = source_lc(t, 0) - x;
		delta_y = source_lc(t, 1) - y;
		delta_z = source_lc(t, 2) - z;

		double distortion = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
		if (distortion < min_distortion)
		{
			t_fix = t;
			min_distortion = distortion;
		}

		a0 << x0 - x1, x1 - x0, 0, y0 - y1, y1 - y0, 0;
		a1 << x0 - x1, 0, x1 - x0, y0 - y1, 0, y1 - y0;
		a2 << y1 - y0, 0, y0 - y1, x0 - x1, 0, x1 - x0;
		a3 << y0 - y1, y1 - y0, 0, x1 - x0, x0 - x1, 0;

		u = a0 / x;
		v = (a1 - z * a3 / x) / x;
		w = (a2 + y * a3 / x) / x;

		h = u * u.transpose() + v * v.transpose() + w * w.transpose();
		b = u * delta_x + v * delta_y + w * delta_z;

		index << f0, f1, f2, f0 + nV, f1 + nV, f2 + nV;
		for (size_t i = 0; i < 6; i++)
		{
			B(index[i]) += b[i];
			for (size_t j = 0; j < 6; j++)
				A.coeffRef(index[i], index[j]) += h(i, j);
		}
	}
	f0 = F0[t_fix];
	f1 = F1[t_fix];
	A.coeffRef(f0, f0) += 1;
	A.coeffRef(f0 + nV, f0 + nV) += 1;
	A.coeffRef(f1 + nV, f1 + nV) += 1;
	A.makeCompressed();
	// solve sparse linear equations
	SparseLU<SparseMatrix<double, RowMajor>> solver;
	solver.compute(A);
	velocity = solver.solve(B);
}


void VLG::UpdateTriMesh()
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
	double max_y = para_map_(0, 1);

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

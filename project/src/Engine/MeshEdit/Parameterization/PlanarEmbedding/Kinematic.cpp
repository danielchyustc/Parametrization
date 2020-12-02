#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/Kinematic.h>
#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SVD> 

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool Kinematic::Run() {
	if (heMesh->IsEmpty() || !triMesh)
	{
		printf("ERROR::Minsurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	InitPara();
	InitSolver();
	init_status_ = true;
	finish_status_ = false;

	UpdateTriMesh();
	return true;
}

void Kinematic::InitPara()
{
	F0.resize(nT);
	F1.resize(nT);
	F2.resize(nT);
	input_mesh_pos.resize(3, nV);
	position.resize(2 * nV);
	velocity.resize(2 * nV);
	acceleration.resize(2 * nV);
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
	acceleration.setZero();

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
	Pre_calculate();
	Tutte();
}

void Kinematic::Pre_calculate()
{
	source_p.resize(nT);
	for (size_t t = 0; t < nT; t++)
	{
		Matrix2d p;
		local_coordinate_inverse(t, p);
		source_p[t] = p;
	}
}

void Kinematic::local_coordinate_inverse(int t, Matrix2d& p)
{
	int f0 = F0[t];
	int f1 = F1[t];
	int f2 = F2[t];
	
	Vector3d e01 = input_mesh_pos.col(f1) - input_mesh_pos.col(f0);
	Vector3d e02 = input_mesh_pos.col(f2) - input_mesh_pos.col(f0);
	Vector3d n_(e01.cross(e02).normalized());
	Vector3d x_ = e01.normalized();
	Vector3d y_(n_.cross(x_));

	Matrix2d X;
	X(0, 0) = e01.norm();	X(0, 1) = e02.dot(x_);
	X(1, 0) = 0;			X(1, 1) = e02.dot(y_);
	p = X.inverse();
}

void Kinematic::Tutte()
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

void Kinematic::InitSolver()
{
	energy_area_process.reserve(10000);
	Energysource();

	energy_area_process.push_back(energy_area);
	energy_pre = 0;
	energy_cur = energy_uniform;

	iter_count_ = 0;
	delta_t = 0.001;
	m = 2;
}

void Kinematic::Energysource()
{
	double end_e_one_temp = 0, end_e_area = 0;

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_1, E_2;

	Matrix2d P, Q, J;

	const double* pos = position.data();
	for (int t = 0; t < nT; t++)
	{
		f0 = F0[t];
		f1 = F1[t];
		f2 = F2[t];

		x0 = pos[f0];
		y0 = pos[f0 + nV];

		x1 = pos[f1];
		y1 = pos[f1 + nV];

		x2 = pos[f2];
		y2 = pos[f2 + nV];

		Q(0,0) = x1 - x0;	Q(0,1) = x2 - x0;
		Q(1,0) = y1 - y0;	Q(1,1) = y2 - y0;

		P = source_p[t];
		J = Q * P;

		det = J.determinant();
		E_1 = J.squaredNorm();
		E_2 = E_1 / (det * det);

		end_e_one_temp += E_1 + E_2;
		end_e_area += (E_1 + E_2) * area_weight[t];
	}
	energy_uniform = end_e_one_temp / nT;
	energy_area = end_e_area;
}

bool Kinematic::Iterate()
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

	iter_count_++;
	for (int i = 1; i <= 1000; i++)
	{
		calc_accel();
//		cout << "acceleration:\t" << acceleration.transpose() << endl;
		position += velocity * delta_t + 0.5 * acceleration * delta_t * delta_t;
		velocity += acceleration * delta_t;
		velocity *= 0.99;
	}

	cout << iter_count_ << "th iteration done!" << endl;
	for (size_t i = 0; i < nV; i++)
	{
		para_map_(i, 0) = position(i);
		para_map_(i, 1) = position(i + nV);
	}
	UpdateTriMesh();
	return true;
}

void Kinematic::calc_accel()
{
	double area_now;
	int f0, f1, f2;
	Matrix2d P, Q, J;

	double x0, y0, x1, y1, x2, y2;
	double lim_factor;
	double det, tr, distortion;
	double r0, r1, r2, r3;
	double d00, d01, d02,
		d10, d11, d12;

	cout << "flipped or largely distorted:\t";
	for (int i = 0; i < nT; ++i)
	{
		area_now = area_weight[i];
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = position[f0];
		y0 = position[f0 + nV];

		x1 = position[f1];
		y1 = position[f1 + nV];

		x2 = position[f2];
		y2 = position[f2 + nV];

		Q(0, 0) = x1 - x0;	Q(0, 1) = x2 - x0;
		Q(1, 0) = y1 - y0;	Q(1, 1) = y2 - y0;

		P = source_p[i];
		J = Q * P;

		det = J.determinant();
		tr = J.squaredNorm();
		distortion = abs(tr * (1.0 + 1.0 / (det * det)));
		if (det < 0) cout << i << " f\t";
		if (distortion > 1000) cout << i << " ld\t";

		lim_factor = exp(-distortion / 1000);

		d00 = -P(0, 0) - P(1, 0); d01 = P(0, 0); d02 = P(1, 0);
		d10 = -P(0, 1) - P(1, 1); d11 = P(0, 1); d12 = P(1, 1);

		r0 = 2 * area_now * lim_factor * ((1 + 1 / (det * det)) * J(0, 0) - tr * J(1, 1) / (det * det * det));
		r1 = 2 * area_now * lim_factor * ((1 + 1 / (det * det)) * J(0, 1) + tr * J(1, 0) / (det * det * det));
		r2 = 2 * area_now * lim_factor * ((1 + 1 / (det * det)) * J(1, 0) + tr * J(0, 1) / (det * det * det));
		r3 = 2 * area_now * lim_factor * ((1 + 1 / (det * det)) * J(1, 1) - tr * J(0, 0) / (det * det * det));

		acceleration(f0) -= r0 * d00 + r1 * d10;
		acceleration(f1) -= r0 * d01 + r1 * d11;
		acceleration(f2) -= r0 * d02 + r1 * d12;
		acceleration(f0 + nV) -= r2 * d00 + r3 * d10;
		acceleration(f1 + nV) -= r2 * d01 + r3 * d11;
		acceleration(f2 + nV) -= r2 * d02 + r3 * d12;
	}
	cout << endl;
	acceleration /= m;
}


void Kinematic::UpdateTriMesh()
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

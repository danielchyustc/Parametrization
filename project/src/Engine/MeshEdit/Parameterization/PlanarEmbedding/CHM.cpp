#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/CHM.h>
#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SVD> 

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool CHM::Run() {
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

void CHM::InitPara()
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

	energy = 0;
	neg_gradient.resize(2 * nV);
}

void CHM::source_local_coordinate()
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

void CHM::Tutte()
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

bool CHM::Iterate()
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
	
	for (size_t i = 0; i < 10; i++)
	{
		negative_gradient();
		solve_velocity();
		velocity.normalize();
		double temp_t;
		max_step(position, velocity, temp_t);

		double alpha = min(1.0, 0.8 * temp_t);
		backtracking_line_search(position, velocity, neg_gradient, alpha);
		position += alpha * velocity;
		cout << alpha << endl;
	}

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

void CHM::solve_velocity()
{
	SparseMatrix<double, RowMajor> A(2 * nV, 2 * nV);
	int f0, f1, f2;
	double x0, x1, x2, y0, y1, y2;
	double p0, p1, p2;
	double x, y, z;
	VectorXd a0(6), a1(6), a2(6), a3(6);
	MatrixXd h(6, 6), u(6, 3);
	Matrix3d g1, g2, g;
	VectorXi index(6);
	A.setIdentity();
	A /= nT;
	for (size_t t = 0; t < nT; t++)
	{
		f0 = F0[t];
		f1 = F1[t];
		f2 = F2[t];

		x0 = position[f0];	y0 = position[f0 + nV];
		x1 = position[f1];	y1 = position[f1 + nV];
		x2 = position[f2];	y2 = position[f2 + nV];

		p0 = source_lc(t, 0);
		p1 = source_lc(t, 1);
		p2 = source_lc(t, 2);

		x = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		y = ((x1 - x0) * (x2 - x0) + (y1 - y0) * (y2 - y0)) / x;
		z = ((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)) / x;

		a0 << x0 - x1, x1 - x0, 0, y0 - y1, y1 - y0, 0;
		a1 << x0 - x1, 0, x1 - x0, y0 - y1, 0, y1 - y0;
		a2 << y1 - y0, 0, y0 - y1, x0 - x1, 0, x1 - x0;
		a3 << y0 - y1, y1 - y0, 0, x1 - x0, x0 - x1, 0;

		u.col(0) = a0 / x;
		u.col(1) = (a1 - z * a3 / x) / x;
		u.col(2) = (a2 + y * a3 / x) / x;

		g1 << z/x, 0, 0,
			  0,   0, 0,
			  0,   0, x/z;
		g2 << y*y/(x*z), -y/z, 0,
			  -y/z,      x/z,  0,
			  0,         0,    0;
		g = 4 * g1 + g2;
		h = u * g * u.transpose();

		index << f0, f1, f2, f0 + nV, f1 + nV, f2 + nV;
		for (size_t i = 0; i < 6; i++)
			for (size_t j = 0; j < 6; j++)
				A.coeffRef(index[i], index[j]) += h(i, j);
	}
	/*
	f0 = F0[0];
	f1 = F1[0];
	A.coeffRef(f0, f0) += 1;
	A.coeffRef(f0 + nV, f0 + nV) += 1;
	A.coeffRef(f1 + nV, f1 + nV) += 1;
	*/
	A.makeCompressed();
	// solve sparse linear equations
	SparseLU<SparseMatrix<double, RowMajor>> solver;
	solver.compute(A);
	velocity = solver.solve(neg_gradient);
}

void CHM::sym_dir_energy()
{
	int f0, f1, f2;
	double x0, x1, x2, y0, y1, y2;
	double p0, p1, p2;
	Matrix2d P, Q, J;
	double tr, det;
	energy = 0;
	for (size_t t = 0; t < nT; t++)
	{
		f0 = F0[t];
		f1 = F1[t];
		f2 = F2[t];

		x0 = position(f0);	y0 = position(f0 + nV);
		x1 = position(f1);	y1 = position(f1 + nV);
		x2 = position(f2);	y2 = position(f2 + nV);

		p0 = 1 / source_lc(t, 0);
		p1 = -source_lc(t, 1) / (source_lc(t, 0) * source_lc(t, 2));
		p2 = 1 / source_lc(t, 2);

		P << p0, p1,
			0, p2;
		Q << x1 - x0, x2 - x0,
			y1 - y0, y2 - y0;
		J = Q * P;

		det = J.determinant();
		tr = J.squaredNorm();
		energy += area_weight[t] * (1.0 + 1.0 / (det * det)) * tr;
	}
}

void CHM::sym_dir_energy(VectorXd pos, double& e)
{
	int f0, f1, f2;
	double x0, x1, x2, y0, y1, y2;
	double p0, p1, p2;
	Matrix2d P, Q, J;
	double tr, det;
	e = 0;
	for (size_t t = 0; t < nT; t++)
	{
		f0 = F0[t];
		f1 = F1[t];
		f2 = F2[t];

		x0 = pos(f0);	y0 = pos(f0 + nV);
		x1 = pos(f1);	y1 = pos(f1 + nV);
		x2 = pos(f2);	y2 = pos(f2 + nV);

		p0 = 1 / source_lc(t, 0);
		p1 = -source_lc(t, 1) / (source_lc(t, 0) * source_lc(t, 2));
		p2 = 1 / source_lc(t, 2);

		P << p0, p1,
			0, p2;
		Q << x1 - x0, x2 - x0,
			y1 - y0, y2 - y0;
		J = Q * P;

		det = J.determinant();
		tr = J.squaredNorm();
		e += area_weight[t] * (1.0 + 1.0 / (det * det)) * tr;
	}
}

void CHM::negative_gradient()
{
	int f0, f1, f2;
	double x0, x1, x2, y0, y1, y2;
	double p0, p1, p2;
	Matrix2d P, Q, J;
	double tr, det;
	double K1, K2;
	double d00, d01, d10, d11;
	for (size_t t = 0; t < nT; t++)
	{
		f0 = F0[t];
		f1 = F1[t];
		f2 = F2[t];

		x0 = position(f0);	y0 = position(f0 + nV);
		x1 = position(f1);	y1 = position(f1 + nV);
		x2 = position(f2);	y2 = position(f2 + nV);

		p0 = 1 / source_lc(t, 0);
		p1 = -source_lc(t, 1) / (source_lc(t, 0) * source_lc(t, 2));
		p2 = 1 / source_lc(t, 2);

		P << p0, p1,
			  0, p2;
		Q << x1 - x0, x2 - x0,
			 y1 - y0, y2 - y0;
		J = Q * P;

		det = J.determinant();
		tr = J.squaredNorm();
		K1 = 2 * (1.0 + 1.0 / (det * det)) * area_weight[t];
		K2 = 2 * tr / (det * det * det) * area_weight[t];

		d00 = J(0, 0) * K1 - J(1, 1) * K2;	d01 = J(0, 1) * K1 + J(1, 0) * K2;
		d10 = J(1, 0) * K1 + J(0, 1) * K2;	d11 = J(1, 1) * K1 - J(0, 0) * K2;

		neg_gradient(f0) += d00 * p0 + d01 * (p1 + p2);
		neg_gradient(f1) -= d00 * p0 + d01 * p1;
		neg_gradient(f2) -= d01 * p2;
		neg_gradient(f0 + nV) += d10 * p0 + d11 * (p1 + p2);
		neg_gradient(f1 + nV) -= d10 * p0 + d11 * p1;
		neg_gradient(f2 + nV) -= d11 * p2;
	}
	neg_gradient.normalize();
}

void CHM::max_step(const VectorXd& xx, const VectorXd& dd, double& step)
{
	double temp_t = numeric_limits<double>::infinity();
	int f0, f1, f2;
	double a, b, c, b1, b2, tt, tt1, tt2;
	double x0, x1, x2, x3, x4, x5, d0, d1, d2, d3, d4, d5;
	const double* x = xx.data();
	const double* d = dd.data();
	for (int i = 0; i < nT; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = x[f0]; x1 = x[f1]; x2 = x[f2]; x3 = x[f0 + nV]; x4 = x[f1 + nV]; x5 = x[f2 + nV];
		d0 = d[f0]; d1 = d[f1]; d2 = d[f2]; d3 = d[f0 + nV]; d4 = d[f1 + nV]; d5 = d[f2 + nV];

		a = (d1 - d0) * (d5 - d3) - (d4 - d3) * (d2 - d0);
		b1 = (d1 - d0) * (x5 - x3) + (x1 - x0) * (d5 - d3);
		b2 = (x4 - x3) * (d2 - d0) + (x2 - x0) * (d4 - d3);
		b = b1 - b2;
		c = (x1 - x0) * (x5 - x3) - (x4 - x3) * (x2 - x0);
		tt = numeric_limits<double>::infinity();
		//tt = 10000;
		if (b * b - 4 * a * c >= 0)
		{
			tt1 = 1 / (2 * a) * (-b + sqrt(b * b - 4 * a * c));
			tt2 = 1 / (2 * a) * (-b - sqrt(b * b - 4 * a * c));
			if (tt1 > 0 && tt2 > 0)
			{
				tt = min(tt1, tt2);
			}
			if (tt1 > 0 && tt2 < 0)
			{
				tt = tt1;
			}
			if (tt1 < 0 && tt2 > 0)
			{
				tt = tt2;
			}
		}
		if (temp_t > tt)
		{
			temp_t = tt;
		}

	}
	step = temp_t;
}

void CHM::backtracking_line_search(const VectorXd& x, const VectorXd& d, const VectorXd& negetive_grad, double& alpha)
{
	double h = 0.5;
	double tt = -(negetive_grad.transpose() * d)(0, 0);
	double c = 0.2;
	double ex;
	sym_dir_energy(x, ex);
	double e;
	VectorXd x_new = x + alpha * d;
	sym_dir_energy(x_new, e);
	while (e > ex + alpha * c * tt)
	{
		alpha = h * alpha;
		x_new = x + alpha * d;
		sym_dir_energy(x_new, e);
	}
}

void CHM::UpdateTriMesh()
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

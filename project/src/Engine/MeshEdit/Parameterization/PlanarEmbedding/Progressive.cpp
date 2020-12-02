#include <Engine/MeshEdit/Parameterization/PlanarEmbedding/Progressive.h>
#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SVD> 

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool Progressive::Run() {
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

void Progressive::InitPara()
{
	F0.resize(nT);
	F1.resize(nT);
	F2.resize(nT);
	input_mesh_pos.resize(3, nV);
	position_of_mesh.resize(2 * nV);
	negative_grad_norm.resize(2 * nV);
	para_map_.resize(nV, 2);

	area.resize(nT);
	area_uniform.resize(nT);
	area_src.resize(nT);

	pardiso = NULL;
	convgence_con_rate = 1e-6;
	MAX_ITER_NUM = 5000;
	bound_distortion_K = 250;

	double area_sum = 0.0;
	for (size_t t = 0; t < nT; t++)
	{
		double area_t = TriArea(t);
		area_sum += area_t;
		area_src[t] = area_t;
	}
	area_src /= area_sum;
	originmesh_area_sqrt = sqrt(area_sum);
	area_uniform.setConstant(1.0 / nT);
	area.setConstant(1.0 / nT);

	double area_same_factor = 1.0 / originmesh_area_sqrt;
	for (size_t i = 0; i < nV; i++)
	{
		vecf3 pos = heMesh->Vertices().at(i)->pos;
		input_mesh_pos.col(i) = Vector3d(pos[0], pos[1], pos[2]) * area_same_factor;
//		cout << input_mesh_pos.col(i) << endl;
	}
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

void Progressive::Pre_calculate()
{
	source_p.resize(nT);
	for (size_t t = 0; t < nT; t++)
	{
		Matrix2d p;
		local_coordinate_inverse(t, p);
		source_p[t] = p;
	}
	update_p = source_p;

	pardiso_i.clear(); pardiso_i.reserve(2 * nV + 1);
	pardiso_ia.clear(); pardiso_ia.reserve(2 * nV + 1);
	pardiso_ja.clear(); pardiso_ja.reserve(8 * nV);

	typedef Triplet<int> T;
	std::vector<T> tripletlist;
	for (int i = 0; i < 2 * nV; i++)
	{
		pardiso_ia.push_back(pardiso_ja.size());
		if (i < nV)
		{
			auto vi = heMesh->Vertices().at(i);
			vector<int> row_id;

			row_id.push_back(i);
			row_id.push_back(i + nV);

			for (auto vj : vi->AdjVertices())
			{
				int id_neighbor = heMesh->Index(vj);
				row_id.push_back(id_neighbor);
				row_id.push_back(id_neighbor + nV);
			}
			std::sort(row_id.begin(), row_id.end(), less<int>());
			vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);

			int dd = 0;

			for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
			{
				pardiso_ja.push_back(row_id[k]);
				pardiso_i.push_back(i);
				tripletlist.push_back(T(i, row_id[k], dd));
				++dd;
			}
		}
		else
		{
			auto vi = heMesh->Vertices().at(i - nV);

			vector<int> row_id;

			row_id.push_back(i);

			for (auto vj : vi->AdjVertices())
			{
				int id_neighbor = heMesh->Index(vj) + nV;
				row_id.push_back(id_neighbor);
			}
			std::sort(row_id.begin(), row_id.end(), less<int>());
			vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);

			int dd = 0;

			for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
			{
				pardiso_ja.push_back(row_id[k]);
				pardiso_i.push_back(i);
				tripletlist.push_back(T(i, row_id[k], dd));
				++dd;
			}
		}
	}
	SparseMatrix<int> find_id_in_rows;
	find_id_in_rows.resize(2 * nV, 2 * nV);
	find_id_in_rows.setFromTriplets(tripletlist.begin(), tripletlist.end());

	pardiso_ia.push_back(pardiso_ja.size());

	id_h00.resize(nT); id_h01.resize(nT); id_h02.resize(nT); id_h03.resize(nT); id_h04.resize(nT); id_h05.resize(nT);
	id_h11.resize(nT); id_h12.resize(nT); id_h13.resize(nT); id_h14.resize(nT); id_h15.resize(nT);
	id_h22.resize(nT); id_h23.resize(nT); id_h24.resize(nT); id_h25.resize(nT);
	id_h33.resize(nT); id_h34.resize(nT); id_h35.resize(nT);
	id_h44.resize(nT); id_h45.resize(nT);
	id_h55.resize(nT);

	for (int i = 0; i < nT; i++)
	{
		int f0 = F0[i]; int f1 = F1[i]; int f2 = F2[i]; int f3 = F0[i] + nV; int f4 = F1[i] + nV; int f5 = F2[i] + nV;

		int min01 = min(f0, f1); int max01 = f0 + f1 - min01;
		int min02 = min(f0, f2); int max02 = f0 + f2 - min02;
		int min12 = min(f1, f2); int max12 = f1 + f2 - min12;

		id_h00[i] = pardiso_ia[f0]; id_h01[i] = pardiso_ia[min01] + find_id_in_rows.coeff(min01, max01); id_h02[i] = pardiso_ia[min02] + find_id_in_rows.coeff(min02, max02);
		id_h03[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f3); id_h04[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f4); id_h05[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f5);

		id_h11[i] = pardiso_ia[f1]; id_h12[i] = pardiso_ia[min12] + find_id_in_rows.coeff(min12, max12);
		id_h13[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f3); id_h14[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f4); id_h15[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f5);

		id_h22[i] = pardiso_ia[f2];
		id_h23[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f3); id_h24[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f4); id_h25[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f5);

		id_h33[i] = pardiso_ia[f3]; id_h34[i] = pardiso_ia[min01 + nV] + find_id_in_rows.coeff(min01 + nV, max01 + nV); id_h35[i] = pardiso_ia[min02 + nV] + find_id_in_rows.coeff(min02 + nV, max02 + nV);

		id_h44[i] = pardiso_ia[f4]; id_h45[i] = pardiso_ia[min12 + nV] + find_id_in_rows.coeff(min12 + nV, max12 + nV);

		id_h55[i] = pardiso_ia[f5];

	}
}

void Progressive::local_coordinate_inverse(int t, Matrix2d& p)
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

void Progressive::Tutte()
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
		B(i, 0) = area_coeff * cos(k * 2 * 3.1415926 / nB);
		B(i, 1) = -area_coeff * sin(k * 2 * 3.1415926 / nB);
	}

	// solve sparse linear equations
	SparseLU<SparseMatrix<double>> solver;
	solver.compute(A);
	para_map_ = solver.solve(B);

	for (size_t i = 0; i < nV; i++)
	{
		position_of_mesh(i) = para_map_(i, 0);
		position_of_mesh(i + nV) = para_map_(i, 1);
	}
}

void Progressive::InitSolver()
{
	if (pardiso != NULL)
	{
		delete pardiso;
		pardiso = NULL;
	}
	pardiso = new PardisoSolver();
	pardiso->ia = pardiso_ia;
	pardiso->ja = pardiso_ja;
	pardiso->a.resize(pardiso_ja.size());
	pardiso->nnz = pardiso_ja.size();
	pardiso->num = 2 * nV;

	pardiso->pardiso_init();

	energy_area_process.reserve(MAX_ITER_NUM);
	Energysource();

	energy_area_process.push_back(energy_area);
	energy_pre = 0;
	energy_cur = energy_uniform;

	Intp_T_Min = 0;
	changetocm_flag = 0;
	flag_1 = false;
	flag_2 = false;

	iter_count_ = 0;
	slim_iter_num = 0;
	cm_iter_num = 0;
	sum_iter_num = 0;

	conv_percent = 1;

	g_norm = 1.0;
}

void Progressive::Energysource()
{
	double end_e_one_temp = 0, end_e_area = 0;

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_1, E_2;

	Matrix2d P, Q, J;

	const double* pos = position_of_mesh.data();
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
		end_e_area += (E_1 + E_2) * area_src[t];
	}
	energy_uniform = end_e_one_temp / nT;
	energy_area = end_e_area;
}

bool Progressive::Iterate()
{
	if (iter_count_ >= MAX_ITER_NUM)
	{
		cout << "Have reached maximal iteration number!" << endl;
		delete pardiso;
		pardiso = NULL;
		return true;
	}
	cout << "1111" << endl;
	if (finish_status_)
	{
		cout << "Already finished!" << endl;
		delete pardiso;
		pardiso = NULL;
		return true;
	}
	cout << "2222" << endl;
	if (changetocm_flag < 0.99 && conv_percent > 0.1 && Intp_T_Min < 0.999 && !flag_1)
	{
		cout << "3333" << endl;
		iter_count_++;
		energy_pre = energy_cur;
		Update_source_same_t();
		SLIM();
		energy_area_process.push_back(energy_area);
		slim_iter_num++;
		sum_iter_num++;
		cout << "4444" << endl;
		energy_cur = energy_uniform;
		conv_percent = abs(energy_cur - energy_pre) / energy_pre;
		calc_gradient_norm(position_of_mesh);
		if (conv_percent <= convgence_con_rate || g_norm <= convgence_con_rate)
		{
			flag_1 = true;
		}
		cout << "5555" << endl;
	}
	else if (conv_percent > 0.01 && Intp_T_Min < 0.999)
	{
		cout << "6666" << endl;
		iter_count_++;
		energy_pre = energy_cur;
		Update_source_same_t();
		CM();
		energy_area_process.push_back(energy_area);
		cm_iter_num++;
		sum_iter_num++;
		cout << "7777" << endl;
		energy_cur = energy_uniform;
		conv_percent = abs(energy_cur - energy_pre) / energy_pre;
		calc_gradient_norm(position_of_mesh);
		if (conv_percent <= convgence_con_rate || g_norm <= convgence_con_rate)
		{
			finish_status_ = true;
		}
		cout << "8888" << endl;
	}
	else
	{
		cout << "9999" << endl;
		if (!flag_2)
		{
			recover_to_src();
			energy_cur = energy_area;
			flag_2 = true;
		}

		iter_count_++;
		energy_pre = energy_cur;
		CM();
		energy_area_process.push_back(energy_area);
		sum_iter_num++;

		energy_cur = energy_area;
		conv_percent = abs(energy_cur - energy_pre) / energy_pre;
		calc_gradient_norm(position_of_mesh);
		if (conv_percent <= convgence_con_rate || g_norm <= convgence_con_rate)
		{
			finish_status_ = true;
		}
	}

	cout << iter_count_ << "th iteration done!" << endl;
	for (size_t i = 0; i < nV; i++)
	{
		para_map_(i, 0) = position_of_mesh(i);
		para_map_(i, 1) = position_of_mesh(i + nV);
		//		cout << para_map_(i, 0) << "\t" << para_map_(i, 1) << endl;
	}
	UpdateTriMesh();
	return true;
}

void Progressive::Update_source_same_t()
{
	double t_min = 1;
	int geqK = -1;

	vector<double> all_s0(nT), all_s1(nT);

	vector<Matrix2d> all_w(nT);

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det;
	double E_d;
	double tt;
	double new_sig0, new_sig1;
	Matrix2d P, Q, J;

	double* position = position_of_mesh.data();
	//	for (size_t i = 0; i < nV; i++)
	//		cout << "position" << position[i] << "\t" << position[i + nV] << endl;

	for (int i = 0; i < nT; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = position[f0];
		y0 = position[f0 + nT];

		x1 = position[f1];
		y1 = position[f1 + nT];

		x2 = position[f2];
		y2 = position[f2 + nT];

		Q(0,0) = x1 - x0;	Q(0,1) = x2 - x0;
		Q(1,0) = y1 - y0;	Q(1,1) = y2 - y0;

		P = source_p[i];
		J = Q * P;

		det = J.determinant();
		E_d = (1 + 1 / (det * det)) * J.squaredNorm();

		JacobiSVD<Matrix2d> svd(J, ComputeFullU | ComputeFullV);
		Vector2d Sig = svd.singularValues();
		double sig0 = Sig(0);
		double sig1 = Sig(1);
		all_s0[i] = sig0;
		all_s1[i] = sig1;
//		cout << "sig" << sig0 << "\t" << sig1 << endl;

		double temp = 1 / (sig1 * sig1 - sig0 * sig0);
		all_w[i] = temp * J.transpose() * J - 0.5 * temp * (sig0 * sig0 + sig1 * sig1) * Matrix2d::Identity();

		if (E_d <= bound_distortion_K) geqK++;
		else
		{
			tt = newton_equation(sig0, sig1, bound_distortion_K);
			if (tt < t_min) t_min = tt;
		}
	}

	changetocm_flag = (geqK + 1.0) / nT;

	for (int i = 0; i < nT; ++i)
	{
		double sig0 = all_s0[i];
		double sig1 = all_s1[i];

		new_sig0 = pow(sig0, t_min - 1);
		new_sig1 = pow(sig1, t_min - 1);

		double delta_new = new_sig1 - new_sig0;
		double plus_new = 0.5 * (new_sig1 + new_sig0);

		Matrix2d W = delta_new * all_w[i] + Matrix2d::Identity() * plus_new;
		P = source_p[i];
		update_p[i] = P * W;
//		cout << "pii" << update_p[i] << endl;
	}

	Intp_T_Min = t_min;
}

void Progressive::SLIM()
{
	double area_now;
	int f0, f1, f2;
	Matrix2d P, Q, J;

	double x0, y0, x1, y1, x2, y2;

	double sig0, sig1;

	double det, tr;
	double r0, r1, r2, r3;
	double d00, d01, d02,
		d10, d11, d12;

	double new_sig0, new_sig1;
	double temp;
	Matrix2d W;
	double p1, p2, p3, w1, w2, w3;

	double h00, h01, h02, h03, h04, h05,
		h11, h12, h13, h14, h15,
		h22, h23, h24, h25,
		h33, h34, h35,
		h44, h45,
		h55;
	double* position = position_of_mesh.data();

	//	for (size_t i = 0; i < nV; i++)
	//		cout << "position" << position[i] << "\t" << position[i + nV] << endl;

	int nnz = pardiso_ja.size();
	pardiso_a.clear(); pardiso_b.clear();
	pardiso_a.resize(nnz, 0.0);
	pardiso_b.resize(2 * nV, 0.0);

	for (int i = 0; i < nT; ++i)
	{
		area_now = area[i];
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
		
		JacobiSVD<Matrix2d> svd(J, ComputeFullU | ComputeFullV);
		Vector2d Sig = svd.singularValues();
		sig0 = Sig(0);
		sig1 = Sig(1);

		new_sig0 = sqrt(1 + 1 / sig0 + 1 / (sig0 * sig0) + 1 / (sig0 * sig0 * sig0)); new_sig1 = sqrt(1 + 1 / sig1 + 1 / (sig1 * sig1) + 1 / (sig1 * sig1 * sig1));

		temp = (new_sig1 - new_sig0) / (sig1 * sig1 - sig0 * sig0);
		W = temp * J.transpose() * J + 0.5 * ((new_sig0 + new_sig1) - (sig0 * sig0 + sig1 * sig1) * temp) * Matrix2d::Identity();

		p1 = P.row(0).squaredNorm(); p2 = P.row(0).dot(P.row(1)); p3 = P.row(1).squaredNorm();
		w1 = W.col(0).squaredNorm(); w2 = W.col(0).dot(W.col(1)); w3 = W.col(1).squaredNorm();

		area_now *= 2;

		h00 = area_now * (p1 + p2 + p2 + p3) * w1; h01 = -area_now * (p1 + p2) * w1; h02 = -area_now * (p2 + p3) * w1; h03 = area_now * (p1 + p2 + p2 + p3) * w2; h04 = -area_now * (p1 + p2) * w2; h05 = -area_now * (p2 + p3) * w2;
		h11 = area_now * p1 * w1;                  h12 = area_now * p2 * w1;    	 h13 = -area_now * (p1 + p2) * w2; h14 = area_now * p1 * w2;                  h15 = area_now * p2 * w2;
		h22 = area_now * p3 * w1;                  h23 = -area_now * (p2 + p3) * w2; h24 = area_now * p2 * w2;         h25 = area_now * p3 * w2;
		h33 = area_now * (p1 + p2 + p2 + p3) * w3; h34 = -area_now * (p1 + p2) * w3; h35 = -area_now * (p2 + p3) * w3;
		h44 = area_now * p1 * w3;                  h45 = area_now * p2 * w3;
		h55 = area_now * p3 * w3;

		det = J.determinant();
		tr = J.squaredNorm();

		d00 = -P(0,0) - P(1,0); d01 = P(0,0); d02 = P(1,0);
		d10 = -P(0,1) - P(1,1); d11 = P(0,1); d12 = P(1,1);

		r0 = area_now * ((1 + 1 / (det * det)) * J(0,0) - tr * J(1,1) / (det * det * det));
		r1 = area_now * ((1 + 1 / (det * det)) * J(0,1) + tr * J(1,0) / (det * det * det));
		r2 = area_now * ((1 + 1 / (det * det)) * J(1,0) + tr * J(0,1) / (det * det * det));
		r3 = area_now * ((1 + 1 / (det * det)) * J(1,1) - tr * J(0,0) / (det * det * det));


		pardiso_b[f0] -= r0 * d00 + r1 * d10;
		pardiso_b[f1] -= r0 * d01 + r1 * d11;
		pardiso_b[f2] -= r0 * d02 + r1 * d12;
		pardiso_b[f0 + nV] -= r2 * d00 + r3 * d10;
		pardiso_b[f1 + nV] -= r2 * d01 + r3 * d11;
		pardiso_b[f2 + nV] -= r2 * d02 + r3 * d12;

		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;

	}

	pardiso->a = pardiso_a;
	pardiso->rhs = pardiso_b;

	pardiso->factorize();
	pardiso->pardiso_solver();

	//	for (size_t i = 0; i < nV; i++)
	//		cout << "position" << position_of_mesh(i) << "\t" << position_of_mesh(i + nV) << endl;

	vector<double> result_d = pardiso->result;

//	for (size_t i = 0; i < nV; i++)
//		cout << "d" << result_d[i] << "\t" << result_d[i + nV] << endl;

	VectorXd negative_grad(2 * nV), d(2 * nV);
	for (int i = 0; i < 2 * nV; i++)
	{
		negative_grad(i) = pardiso_b[i];
		d(i) = result_d[i];
	}

	double temp_t;
	max_step(position_of_mesh, d, temp_t);

	double alpha = min(1.0, 0.8 * temp_t);
	backtracking_line_search(position_of_mesh, d, negative_grad, alpha);
	position_of_mesh += alpha * d;

	Energysource();
}

void Progressive::recover_to_src()
{
	area = area_src;
	update_p = source_p;
}

void Progressive::CM()
{
	double area_now;
	int f0, f1, f2;
	Matrix2d P, Q, J;

	double x0, y0, x1, y1, x2, y2;

	double hi_0, hi_1;

	double alpha_0, alpha_1, beta_0, beta_1;

	double s1, s2, sig0, sig1;

	double alpha_norm, beta_norm;
	double h_u, h_v, walpha, wbeta;

	double a1x0, a1x1, a1x2, a1x3, a1x4, a1x5,
		a2x0, a2x1, a2x2, a2x3, a2x4, a2x5;

	double aa, bb;
	double uu, vv, uv;
	double u, v;

	double h00, h01, h02, h03, h04, h05,
		h11, h12, h13, h14, h15,
		h22, h23, h24, h25,
		h33, h34, h35,
		h44, h45,
		h55;
	double* position = position_of_mesh.data();
	int nnz = pardiso_ja.size();
	pardiso_a.clear(); pardiso_b.clear();
	pardiso_a.resize(nnz, 0.0);
	pardiso_b.resize(2 * nV, 0.0);

	for (int i = 0; i < nT; ++i)
	{
		area_now = area[i];
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

		alpha_0 = J(0,0) + J(1,1); alpha_1 = J(1,0) - J(0,1);
		beta_0 = J(0,0) - J(1,1);  beta_1 = J(1,0) + J(0,1);

		alpha_norm = 0.5 * sqrt(alpha_0 * alpha_0 + alpha_1 * alpha_1);
		beta_norm = 0.5 * sqrt(beta_0 * beta_0 + beta_1 * beta_1);

		s1 = P(0,0) * (P(0,0) + P(1,0)) + P(0,1) * (P(0,1) + P(1,1));
		s2 = P(1,0) * (P(0,0) + P(1,0)) + P(1,1) * (P(0,1) + P(1,1));

		double h1 = P.row(0).squaredNorm();
		double h2 = P.row(0).dot(P.row(1));
		double h3 = P.row(1).squaredNorm();
		double h4 = P.determinant();

		double p00 = P(0, 0), p01 = P(0, 1), p10 = P(1, 0), p11 = P(1,1);

		a1x0 = alpha_0 * (-p00 - p10) + alpha_1 * (p01 + p11);  a1x1 = alpha_0 * p00 - alpha_1 * p01; a1x2 = alpha_0 * p10 - alpha_1 * p11;
		a1x3 = alpha_0 * (-p01 - p11) + alpha_1 * (-p00 - p10); a1x4 = alpha_0 * p01 + alpha_1 * p00; a1x5 = alpha_0 * p11 + alpha_1 * p10;

		a2x0 = beta_0 * (-p00 - p10) + beta_1 * (-p01 - p11);   a2x1 = beta_0 * p00 + beta_1 * p01;   a2x2 = beta_0 * p10 + beta_1 * p11;
		a2x3 = beta_0 * (p01 + p11) + beta_1 * (-p00 - p10);    a2x4 = -beta_0 * p01 + beta_1 * p00;  a2x5 = -beta_0 * p11 + beta_1 * p10;

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		hi_0 = 2 + 6 * 1 / (sig0 * sig0 * sig0 * sig0); hi_1 = 2 + 6 * 1 / (sig1 * sig1 * sig1 * sig1);

		aa = 0.25 / alpha_norm; bb = 0.25 / beta_norm;

		uu = aa * aa * (area_now * hi_0 + area_now * hi_1);
		vv = bb * bb * (area_now * hi_0 + area_now * hi_1);
		uv = aa * bb * (area_now * hi_0 - area_now * hi_1);

		h_u = area_now * (2 * sig0 - 2 * 1 / (sig0 * sig0 * sig0));
		h_v = area_now * (2 * sig1 - 2 * 1 / (sig1 * sig1 * sig1));

		walpha = h_u + h_v;
		wbeta = h_u - h_v;

		double hwa1 = (walpha * 0.25 / alpha_norm); double hwa2 = -(walpha * 0.25 * 0.25 / (alpha_norm * alpha_norm * alpha_norm));
		double hwb1 = (wbeta * 0.25 / beta_norm); double hwb2 = -(wbeta * 0.25 * 0.25 / (beta_norm * beta_norm * beta_norm));


		h00 = uu * a1x0 * a1x0 + vv * a2x0 * a2x0 + uv * a1x0 * a2x0 + uv * a2x0 * a1x0; h01 = uu * a1x0 * a1x1 + vv * a2x0 * a2x1 + uv * a1x0 * a2x1 + uv * a2x0 * a1x1; h02 = uu * a1x0 * a1x2 + vv * a2x0 * a2x2 + uv * a1x0 * a2x2 + uv * a2x0 * a1x2; h03 = uu * a1x0 * a1x3 + vv * a2x0 * a2x3 + uv * a1x0 * a2x3 + uv * a2x0 * a1x3; h04 = uu * a1x0 * a1x4 + vv * a2x0 * a2x4 + uv * a1x0 * a2x4 + uv * a2x0 * a1x4; h05 = uu * a1x0 * a1x5 + vv * a2x0 * a2x5 + uv * a1x0 * a2x5 + uv * a2x0 * a1x5;

		h11 = uu * a1x1 * a1x1 + vv * a2x1 * a2x1 + uv * a1x1 * a2x1 + uv * a2x1 * a1x1; h12 = uu * a1x1 * a1x2 + vv * a2x1 * a2x2 + uv * a1x1 * a2x2 + uv * a2x1 * a1x2; h13 = uu * a1x1 * a1x3 + vv * a2x1 * a2x3 + uv * a1x1 * a2x3 + uv * a2x1 * a1x3; h14 = uu * a1x1 * a1x4 + vv * a2x1 * a2x4 + uv * a1x1 * a2x4 + uv * a2x1 * a1x4; h15 = uu * a1x1 * a1x5 + vv * a2x1 * a2x5 + uv * a1x1 * a2x5 + uv * a2x1 * a1x5;

		h22 = uu * a1x2 * a1x2 + vv * a2x2 * a2x2 + uv * a1x2 * a2x2 + uv * a2x2 * a1x2; h23 = uu * a1x2 * a1x3 + vv * a2x2 * a2x3 + uv * a1x2 * a2x3 + uv * a2x2 * a1x3; h24 = uu * a1x2 * a1x4 + vv * a2x2 * a2x4 + uv * a1x2 * a2x4 + uv * a2x2 * a1x4; h25 = uu * a1x2 * a1x5 + vv * a2x2 * a2x5 + uv * a1x2 * a2x5 + uv * a2x2 * a1x5;

		h33 = uu * a1x3 * a1x3 + vv * a2x3 * a2x3 + uv * a1x3 * a2x3 + uv * a2x3 * a1x3; h34 = uu * a1x3 * a1x4 + vv * a2x3 * a2x4 + uv * a1x3 * a2x4 + uv * a2x3 * a1x4; h35 = uu * a1x3 * a1x5 + vv * a2x3 * a2x5 + uv * a1x3 * a2x5 + uv * a2x3 * a1x5;

		h44 = uu * a1x4 * a1x4 + vv * a2x4 * a2x4 + uv * a1x4 * a2x4 + uv * a2x4 * a1x4; h45 = uu * a1x4 * a1x5 + vv * a2x4 * a2x5 + uv * a1x4 * a2x5 + uv * a2x4 * a1x5;

		h55 = uu * a1x5 * a1x5 + vv * a2x5 * a2x5 + uv * a1x5 * a2x5 + uv * a2x5 * a1x5;

		if (walpha >= 0)
		{
			h00 += hwa1 * (s1 + s2) + hwa2 * a1x0 * a1x0; h01 += hwa1 * (-s1) + hwa2 * a1x0 * a1x1; h02 += hwa1 * (-s2) + hwa2 * a1x0 * a1x2; h03 += hwa2 * a1x0 * a1x3; h04 += hwa1 * (h4)+hwa2 * a1x0 * a1x4; h05 += hwa1 * (-h4) + hwa2 * a1x0 * a1x5;
			h11 += hwa1 * (h1)+hwa2 * a1x1 * a1x1;        h12 += hwa1 * (h2)+hwa2 * a1x1 * a1x2;    h13 += hwa1 * (-h4) + hwa2 * a1x1 * a1x3; h14 += hwa2 * a1x1 * a1x4; h15 += hwa1 * (h4)+hwa2 * a1x1 * a1x5;
			h22 += hwa1 * (h3)+hwa2 * a1x2 * a1x2;        h23 += hwa1 * (h4)+hwa2 * a1x2 * a1x3;    h24 += hwa1 * (-h4) + hwa2 * a1x2 * a1x4; h25 += hwa2 * a1x2 * a1x5;
			h33 += hwa1 * (s1 + s2) + hwa2 * a1x3 * a1x3; h34 += hwa1 * (-s1) + hwa2 * a1x3 * a1x4; h35 += hwa1 * (-s2) + hwa2 * a1x3 * a1x5;
			h44 += hwa1 * (h1)+hwa2 * a1x4 * a1x4;        h45 += hwa1 * (h2)+hwa2 * a1x4 * a1x5;
			h55 += hwa1 * (h3)+hwa2 * a1x5 * a1x5;

		}
		h00 += hwb1 * (s1 + s2) + hwb2 * a2x0 * a2x0; h01 += hwb1 * (-s1) + hwb2 * a2x0 * a2x1; h02 += hwb1 * (-s2) + hwb2 * a2x0 * a2x2; h03 += hwb2 * a2x0 * a2x3; h04 += hwb1 * (-h4) + hwb2 * a2x0 * a2x4; h05 += hwb1 * (h4)+hwb2 * a2x0 * a2x5;
		h11 += hwb1 * (h1)+hwb2 * a2x1 * a2x1;        h12 += hwb1 * (h2)+hwb2 * a2x1 * a2x2;    h13 += hwb1 * (h4)+hwb2 * a2x1 * a2x3;    h14 += hwb2 * a2x1 * a2x4; h15 += hwb1 * (-h4) + hwb2 * a2x1 * a2x5;
		h22 += hwb1 * (h3)+hwb2 * a2x2 * a2x2;        h23 += hwb1 * (-h4) + hwb2 * a2x2 * a2x3; h24 += hwb1 * (h4)+hwb2 * a2x2 * a2x4;    h25 += hwb2 * a2x2 * a2x5;
		h33 += hwb1 * (s1 + s2) + hwb2 * a2x3 * a2x3; h34 += hwb1 * (-s1) + hwb2 * a2x3 * a2x4; h35 += hwb1 * (-s2) + hwb2 * a2x3 * a2x5;
		h44 += hwb1 * (h1)+hwb2 * a2x4 * a2x4;        h45 += hwb1 * (h2)+hwb2 * a2x4 * a2x5;
		h55 += hwb1 * (h3)+hwb2 * a2x5 * a2x5;

		u = aa * walpha; v = bb * wbeta;

		pardiso_b[f0] -= (u * a1x0 + v * a2x0);
		pardiso_b[f1] -= (u * a1x1 + v * a2x1);
		pardiso_b[f2] -= (u * a1x2 + v * a2x2);
		pardiso_b[f0 + nV] -= (u * a1x3 + v * a2x3);
		pardiso_b[f1 + nV] -= (u * a1x4 + v * a2x4);
		pardiso_b[f2 + nV] -= (u * a1x5 + v * a2x5);

		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;
	}

	pardiso->a = pardiso_a;
	pardiso->rhs = pardiso_b;

	pardiso->factorize();
	pardiso->pardiso_solver();

	vector<double> result_d = pardiso->result;

	VectorXd negative_grad(2 * nV), d(2 * nV);
	for (int i = 0; i < 2 * nV; i++)
	{
		negative_grad(i) = pardiso_b[i];
		d(i) = result_d[i];
	}

	pardiso->free_numerical_factorization_memory();

	double temp_t;
	max_step(position_of_mesh, d, temp_t);

	double alpha = 0.95 * temp_t;

	backtracking_line_search(position_of_mesh, d, negative_grad, alpha);

	position_of_mesh += alpha * d;
	Energysource();
}

void Progressive::max_step(const VectorXd& xx, const VectorXd& dd, double& step)
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

void Progressive::calc_gradient_norm(const VectorXd& x)
{
	double area_now;
	int f0, f1, f2;
	Matrix2d P, Q, J;

	double x0, y0, x1, y1, x2, y2;

	double det, tr;
	double r0, r1, r2, r3;
	double d00, d01, d02,
		d10, d11, d12;

	negative_grad_norm.setZero();

	const double* position = x.data();

	for (int i = 0; i < nT; ++i)
	{
		area_now = area[i];
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = position[f0];	y0 = position[f0 + nV];
		x1 = position[f1];	y1 = position[f1 + nV];
		x2 = position[f2];	y2 = position[f2 + nV];

		Q(0, 0) = x1 - x0;	Q(0, 1) = x2 - x0;
		Q(1, 0) = y1 - y0;	Q(1, 1) = y2 - y0;

		P = source_p[i];
		J = Q * P;

		det = J.determinant();
		tr = J.squaredNorm();

		d00 = -P(0, 0) - P(1, 0); d01 = P(0, 0); d02 = P(1, 0);
		d10 = -P(0, 1) - P(1, 1); d11 = P(0, 1); d12 = P(1, 1);

		r0 = 2 * area_now * ((1 + 1 / (det * det)) * J(0, 0) - tr * J(1, 1) / (det * det * det));
		r1 = 2 * area_now * ((1 + 1 / (det * det)) * J(0, 1) + tr * J(1, 0) / (det * det * det));
		r2 = 2 * area_now * ((1 + 1 / (det * det)) * J(1, 0) + tr * J(0, 1) / (det * det * det));
		r3 = 2 * area_now * ((1 + 1 / (det * det)) * J(1, 1) - tr * J(0, 0) / (det * det * det));

		negative_grad_norm(f0) -= r0 * d00 + r1 * d10;
		negative_grad_norm(f1) -= r0 * d01 + r1 * d11;
		negative_grad_norm(f2) -= r0 * d02 + r1 * d12;
		negative_grad_norm(f0 + nV) -= r2 * d00 + r3 * d10;
		negative_grad_norm(f1 + nV) -= r2 * d01 + r3 * d11;
		negative_grad_norm(f2 + nV) -= r2 * d02 + r3 * d12;
	}

	g_norm = negative_grad_norm.norm();
}

void Progressive::backtracking_line_search(const VectorXd& x, const VectorXd& d, const VectorXd& negetive_grad, double& alpha)
{
	double h = 0.5;
	double tt = -(negetive_grad.transpose() * d)(0, 0);
	double c = 0.2;
	double ex;
	Energy(x, ex);
	double e;
	VectorXd x_new = x + alpha * d;
	Energy(x_new, e);
	while (e > ex + alpha * c * tt)
	{
		alpha = h * alpha; 
		x_new = x + alpha * d;
		Energy(x_new, e);
	}
}

void Progressive::Energy(const VectorXd& position, double& energyupdate)
{
	double energy = 0;


	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_d;
	Matrix2d P, Q, J;
	const double* pos = position.data();
	for (int i = 0; i < nT; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + nV];

		x1 = pos[f1];
		y1 = pos[f1 + nV];

		x2 = pos[f2];
		y2 = pos[f2 + nV];

		Q(0, 0) = x1 - x0;	Q(0, 1) = x2 - x0;
		Q(1, 0) = y1 - y0;	Q(1, 1) = y2 - y0;

		P = update_p[i];
		J = Q * P;

		det = J.determinant();
		E_d = (1.0 + 1.0 / (det * det)) * J.squaredNorm();

		energy += area[i] * E_d;
	}
	energyupdate = energy;
}

double Progressive::newton_equation(const double& a, const double& b, const double& K)
{
	double tt = 1;
	double E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K;
	while (abs(E_d) > 1e-5)
	{
		tt = tt - 1 / (2 * log(a) * pow(a, 2 * tt) + 2 * log(b) * pow(b, 2 * tt) + 2 * log(1 / a) * pow(1 / a, 2 * tt) + 2 * log(1 / b) * pow(1 / b, 2 * tt)) * (pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K);
		E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K;
	}
	return tt;
}


void Progressive::UpdateTriMesh()
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

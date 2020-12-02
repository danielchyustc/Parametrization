#include <Engine/MeshEdit/Parameterization/SphericalEmbedding/ProjASAP.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>

using namespace Ubpa;

using namespace std;
using namespace Eigen;

bool ProjASAP::Run() {
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

void ProjASAP::InitParaMap()
{
	// find flat point
	double min_curv_met = 100.0;
	size_t i_flat = 0;
	for (size_t i = 0; i < nV; i++)
	{
		double K = GaussCurvature(i);
		double H = MeanCurvature(i);
		double curv_met = 2 * H * H - K;
		if (curv_met < min_curv_met)
		{
			min_curv_met = curv_met;
			i_flat = i;
		}
	}
	auto vi_flat = heMesh->Vertices().at(i_flat);
	size_t nFN = vi_flat->Degree();
	size_t nV1 = nV - 1;
	size_t nT1 = nT - nFN;

	cout << "Projected ASAP: Flattest point found." << endl;

	// compute area
	double area_i = SurfaceElement(i_flat);
	double total_area = 0;
	for (size_t i = 0; i < nV; i++)
	{
		total_area += VoronoiArea(i);
	}
	double R = sqrt(total_area / (4 * area_i));

	cout << "Projected ASAP: Area computed." << endl;



	// init flat triangles
	flat_tri_ = (Matrix<double, 3, 2>*)malloc(nT1 * sizeof(Matrix<double, 3, 2>));
	if (flat_tri_ == nullptr)
	{
		cout << "No enough memory!" << endl;
	}
	size_t flat_tri_index = 0;
	for (size_t i = 0; i < nT; i++)
	{
		auto ti = heMesh->Polygons().at(i);
		auto vi0 = ti->HalfEdge()->Origin();
		auto vi1 = ti->HalfEdge()->End();
		auto vi2 = ti->HalfEdge()->Next()->End();
		if (vi0 == vi_flat || vi1 == vi_flat || vi2 == vi_flat) continue;
		vecf3 ei1 = vi1->pos - vi0->pos;
		vecf3 ei2 = vi2->pos - vi0->pos;
		double cos_theta = ei1.cos_theta(ei2);
		flat_tri_[flat_tri_index] = Matrix<double, 3, 2>();
		flat_tri_[flat_tri_index].row(0) = RowVector2d::Zero();
		flat_tri_[flat_tri_index].row(1) = RowVector2d(ei1.norm(), 0);
		flat_tri_[flat_tri_index].row(2) = RowVector2d(ei2.norm() * cos_theta, ei2.norm() * sqrt(1 - pow(cos_theta, 2)));
		flat_tri_index++;
	}

	// build sparse system
	SparseMatrix<double, RowMajor> A(2 * nV1 + 2 * nT1, 2 * nV1 + 2 * nT1);
	VectorXd B(2 * nV1 + 2 * nT1);
	A.setZero();
	B.setZero();

	auto v1 = vi_flat->AdjVertices().at(0);
	auto v2 = vi_flat->AdjVertices().at(nFN / 2);
	int i1 = heMesh->Index(v1);
	int i2 = heMesh->Index(v2);
	int insert_i1 = (i1 < i_flat) ? i1 : i1 - 1;
	int insert_i2 = (i2 < i_flat) ? i2 : i2 - 1;

	A.insert(insert_i1, insert_i1) = 1;
	A.insert(nV1 + insert_i1, nV1 + insert_i1) = 1;
	A.insert(insert_i2, insert_i2) = 1;
	A.insert(nV1 + insert_i2, nV1 + insert_i2) = 1;
	B(insert_i1) = 0;
	B(nV1 + insert_i1) = 0;
	B(insert_i2) = R;
	B(nV1 + insert_i2) = R;

	// init A and B

	size_t insert_index = 0;
	
	flat_tri_index = 0;
	for (size_t t = 0; t < nT; t++)
	{
		auto triangle = heMesh->Polygons().at(t);
		auto edge = triangle->HalfEdge();
		auto vi = edge->Origin();
		auto vj = edge->End();
		auto vk = edge->Next()->End();
		if (vi == vi_flat || vj == vi_flat || vk == vi_flat) continue;
		int i = heMesh->Index(vi);
		int j = heMesh->Index(vj);
		int k = heMesh->Index(vk);
		int ix = (i < i_flat) ? i : i - 1;
		int jx = (j < i_flat) ? j : j - 1;
		int kx = (k < i_flat) ? k : k - 1;
		int iy = ix + nV1;
		int jy = jx + nV1;
		int ky = kx + nV1;
		int tx = 2 * nV1 + flat_tri_index;
		int ty = 2 * nV1 + nT1 + flat_tri_index;
		RowVector2d xi = flat_tri_[flat_tri_index].row(0);
		RowVector2d xj = flat_tri_[flat_tri_index].row(1);
		RowVector2d xk = flat_tri_[flat_tri_index].row(2);
		vecf3 eij = vj->pos - vi->pos;
		vecf3 eik = vk->pos - vi->pos;
		vecf3 eji = vi->pos - vj->pos;
		vecf3 ejk = vk->pos - vj->pos;
		vecf3 eki = vi->pos - vk->pos;
		vecf3 ekj = vj->pos - vk->pos;

		double cot_ij = eki.cos_theta(ekj) / eki.sin_theta(ekj);
		double cot_jk = eij.cos_theta(eik) / eij.sin_theta(eik);
		double cot_ki = ejk.cos_theta(eji) / ejk.sin_theta(eji);
		double cot_ijk = cot_ij + cot_jk;
		double cot_jki = cot_jk + cot_ki;
		double cot_kij = cot_ki + cot_ij;

		double cot_ij_x = cot_ij * (xi(0) - xj(0));
		double cot_jk_x = cot_jk * (xj(0) - xk(0));
		double cot_ki_x = cot_ki * (xk(0) - xi(0));
		double cot_ijk_x = cot_ij_x - cot_jk_x;
		double cot_jki_x = cot_jk_x - cot_ki_x;
		double cot_kij_x = cot_ki_x - cot_ij_x;

		double cot_ij_y = cot_ij * (xi(1) - xj(1));
		double cot_jk_y = cot_jk * (xj(1) - xk(1));
		double cot_ki_y = cot_ki * (xk(1) - xi(1));
		double cot_ijk_y = cot_ij_y - cot_jk_y;
		double cot_jki_y = cot_jk_y - cot_ki_y;
		double cot_kij_y = cot_ki_y - cot_ij_y;

		double tn = cot_ij * (xi - xj).squaredNorm() + cot_jk * (xj - xk).squaredNorm() + cot_ki * (xk - xi).squaredNorm();

		A.coeffRef(ix, ix) += cot_kij;
		A.coeffRef(ix, jx) -= cot_ij;
		A.coeffRef(ix, kx) -= cot_ki;
		A.coeffRef(ix, tx) += cot_kij_x;
		A.coeffRef(ix, ty) += cot_kij_y;

		A.coeffRef(iy, iy) += cot_kij;
		A.coeffRef(iy, jy) -= cot_ij;
		A.coeffRef(iy, ky) -= cot_ki;
		A.coeffRef(iy, tx) += cot_kij_y;
		A.coeffRef(iy, ty) -= cot_kij_x;

		A.coeffRef(jx, jx) += cot_ijk;
		A.coeffRef(jx, kx) -= cot_jk;
		A.coeffRef(jx, ix) -= cot_ij;
		A.coeffRef(jx, tx) += cot_ijk_x;
		A.coeffRef(jx, ty) += cot_ijk_y;

		A.coeffRef(jy, jy) += cot_ijk;
		A.coeffRef(jy, ky) -= cot_jk;
		A.coeffRef(jy, iy) -= cot_ij;
		A.coeffRef(jy, tx) += cot_ijk_y;
		A.coeffRef(jy, ty) -= cot_ijk_x;

		A.coeffRef(kx, kx) += cot_jki;
		A.coeffRef(kx, ix) -= cot_ki;
		A.coeffRef(kx, jx) -= cot_jk;
		A.coeffRef(kx, tx) += cot_jki_x;
		A.coeffRef(kx, ty) += cot_jki_y;

		A.coeffRef(ky, ky) += cot_jki;
		A.coeffRef(ky, iy) -= cot_ki;
		A.coeffRef(ky, jy) -= cot_jk;
		A.coeffRef(ky, tx) += cot_jki_y;
		A.coeffRef(ky, ty) -= cot_jki_x;
		
		A.coeffRef(tx, ix) -= cot_kij_x;
		A.coeffRef(tx, jx) -= cot_ijk_x;
		A.coeffRef(tx, kx) -= cot_jki_x;
		A.coeffRef(tx, iy) -= cot_kij_y;
		A.coeffRef(tx, jy) -= cot_ijk_y;
		A.coeffRef(tx, ky) -= cot_jki_y;
		A.coeffRef(tx, tx) -= tn;

		A.coeffRef(ty, ix) -= cot_kij_y;
		A.coeffRef(ty, jx) -= cot_ijk_y;
		A.coeffRef(ty, kx) -= cot_jki_y;
		A.coeffRef(ty, iy) += cot_kij_x;
		A.coeffRef(ty, jy) += cot_ijk_x;
		A.coeffRef(ty, ky) += cot_jki_x;
		A.coeffRef(ty, ty) -= tn;

		flat_tri_index++;
		if (flat_tri_index % 1000 == 0)
			cout << flat_tri_index << " blocks built." << endl;
	}
	A.makeCompressed();

	cout << "Projected ASAP: Sparse linear system built." << endl;

	// solve sparse linear system
	SparseLU<SparseMatrix<double>> solver;
	solver.compute(A);
	VectorXd solution = solver.solve(B);

	cout << "Projected ASAP: Sparse linear system solved." << endl;

	// process planar paramap
	MatrixXd planar_para(2, nV1);
	for (size_t i = 0; i < nV1; i++)
	{
		planar_para(0, i) = solution(i);
		planar_para(1, i) = solution(nV1 + i);
	}

	Vector2d planar_center = Barycenter(planar_para);
	double planar_radius = 0;
	for (size_t i = 0; i < nV1; i++)
	{
		double r = (planar_para.col(i) - planar_center).norm();
		if (r > planar_radius) planar_radius = r;
	}
	for (size_t i = 0; i < nV1; i++)
	{
		planar_para.col(i) = (planar_para.col(i) - planar_center) * R / planar_radius;
	}

	cout << "Projected ASAP: Planar paramap processed." << endl;

	// update paramap
	para_map_.resize(3, nV);
	for (size_t i = 0; i < nV; i++)
	{
		if (i == i_flat)
		{
			para_map_.col(i) = Vector3d(0, 0, 1);
		}
		else
		{
			insert_index = (i < i_flat) ? i : i - 1;
			double x = planar_para(0, insert_index);
			double y = planar_para(1, insert_index);
			double r = sqrt(x * x + y * y);
			//float theta = 4 * (1 - 1 / R) * atan(pow(r/R, 1));
			//para_map_.col(i) = Vector3f(sin(theta) * x / r, sin(theta) * y / r, -cos(theta));
			para_map_.col(i) = Vector3d(2 * x, 2 * y, r * r - 1) / (r * r + 1);
		}
	}

	cout << "Projected ASAP: Parameterization map updated." << endl;

//	MoveBarycenterToOrigin();
	init_status_ = true;
	finish_status_ = true;
}